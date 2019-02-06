// us
#include "ctransport.hh"
#include "thread_pool.hpp"

// standard includes
#include <iostream>
#include <string>
#include <future>
#include <functional>
#include <algorithm>

// ROOT includes
#include "TMath.h"
#include "TFile.h"
#include "TNtupleD.h"


//*******
// Collection transport
//*******
Ctransport::Ctransport(std::string fname, int seed) {
  charges.clear();
  times.clear();
  places.clear();
  density = 0.1664; // [kg/m^3] fix NTP (295K) helium gas density
  rnd = new TRandom3(seed);
  readCS(fname); // fixed CS file name
}


Ctransport::~Ctransport() {
  delete rnd;
}


// calculate a signal on electrode for any charges in region of interest
void Ctransport::ctransport(Electrode* electrode, std::list<charge_t> q) {

  if (q.size()<1) {
    std::cout << "Error: container of charges is empty" << std::endl;
    return;
  }
  // got all charges as initial input
  charges.clear(); // copy to data member
  for (charge_t cc : q) {
    charges.push_front(cc); // insert from front
  }
  
  run(electrode);
  return;
}


double Ctransport::pick_target(std::vector<double> weight, int& which) {
  double target_mass;
  double pick = rnd->Rndm();
  if (pick <= weight[0]) { // Helium
    target_mass = 4.0026 * 0.93149; // [GeV/c^2]
    which = 0;
  }
  else if (pick > weight[0] && pick <= weight[0]+weight[1]) { // Ethanol
    target_mass = 46.069 * 0.93149; // [GeV/c^2]
    which = 1;
  }
  else { // Argon
    target_mass = 39.948 * 0.93149; // [GeV/c^2]
    which = 2;
  }
  return target_mass;
}


bool Ctransport::taskfunction(Electrode* electrode, charge_t q) {
  // have a charge and info about all fields for each thread

  //Init
  TVector3 speed;
  double energy;
  double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
  double time_sum = 0.0;
  double target_mass;
  int inel_flag = 0;
  int which = 0; // Helium
  
  TVector3 distance_sum, distance_step;

  distance_step.SetXYZ(0.0,0.0,0.0);
  speed.SetXYZ(0.0,0.0,0.0);
  
  double e_mass = 0.511e-3; // [GeV/c^2]

  std::vector<double> weight; // [%] gas composition volume ratios
  weight.push_back(0.95);
  weight.push_back(0.04);
  weight.push_back(0.01);

  target_mass = 4.0026 * 0.93149; // Helium default first target [GeV/c^2]
  double mumass_eV = 1.0e9 * e_mass * target_mass / (e_mass + target_mass);

  double prob;
  double localdensity = density * 6.023e26 / target_mass;// convert to number density [m^-3]
  // for E=2.12e8, gives E/N = 10Td = 1.e-16 Vcm^2
  
  double time_step, running_time;
  double kmax = 2.e-12; // constant for null coll. method
  double kv;
  double tau = 1/(localdensity * kmax);
  
  double init_energy;
  double speed_start, tangle;

  // speed vector init
  speed.SetXYZ(-1.0,0.0,0.0);
  init_energy = 1.e-9 * 0.025;// thermal start energy [GeV] 
  speed_start = TMath::Sqrt(2.0*init_energy / e_mass * c2);
  speed.SetMag(speed_start);
  tangle = TMath::Pi()*rnd->Rndm();// isotropic
  speed.SetTheta(tangle);
  
  time_sum = running_time = 0.0;

  double x, y, z;
  double xe, ye, ze;
  int elcharge;
  Point3 exyz; // Drift field
  Point3 point;
  Point3 previous;
  bool analytic = false; // default False

  point = q.location; // start location, Point3 object; [cm] from ROOT
  previous = point;

  // starting TVector3 from point
  distance_sum.SetXYZ(point.xc()*0.01,point.yc()*0.01,point.zc()*0.01); // [cm]->[m]

  elcharge = q.charge; // -1: e-
  charge_t cc;

  exyz = electrode->getFieldValue(analytic,point); // [V/m]

  // debug
  //  int nsteps = 0;

  // transport loop
  while (!analytic) { 
	
    // prepare and update
    time_step = time_update(tau);
    running_time += time_step;
    // keep track of total time
    time_sum += time_step;

    // vector addition stepwise turns velocity vector
    speed += speed_update(elcharge,exyz,time_step);
    
    // CMS system energy
    energy = 0.5*mumass_eV*speed.Mag2()/c2; // non-rel. energy in [eV]
    // artificially raise the cross section 
    kv = speed.Mag() * cross_section(energy, which, inel_flag);

    if (inel_flag>0) { // was ionization
      speed.SetXYZ(0.0,0.0,0.0); // inelastic takes energy off e-
      kv = 0.0;
      cc.location = point; // last known collision location
      cc.chargeID = 1; // was an electron
      cc.charge = -1; //
      book_charge(cc); // store in object container
    }
    
    if (kv>=kmax) {
      std::cout << "kmax too small" << std::endl;
      break;
    }
    
    // random number collision decision
    prob = rnd->Rndm();
	    
    // collision decision
    if (prob <= (kv/kmax)) {
      //      nsteps += 1; // collision occurred
      //      if (!(nsteps % 100000)) std::cout << "collision " << nsteps << " : x,y coordinates " << point.xc() << " " << point.yc() << std::endl;
      //      if (nsteps>=5000) analytic = true; // stop after n steps
      // book position of collision
      distance_step = d_update(speed,running_time);
      distance_sum += distance_step; // in [m]
      point.Set(distance_sum.X()*100.0,distance_sum.Y()*100.0,distance_sum.Z()*100.0); // [cm]

      // new speed from elastic collision kinematics
      speed = kin_factor2(speed, target_mass);
      
      // check geometry and fields
      exyz = electrode->getFieldValue(analytic,point);
      // std::cout << "in transport: field values " << exyz.xc() << " " << exyz.yc() << std::endl;
      // std::cout << "in transport: x,y coordinates " << point.xc() << " " << point.yc() << std::endl;
      // std::cout << "collision at energy " << energy << std::endl;
      // std::cout << "speed X: " << speed.X() << " Y: " << speed.Y() << std::endl;
      // std::cout << "time between coll " << running_time << std::endl;

      if (analytic) {
	book_time(time_sum); // e- stopping, record time
	book_place(previous); // stop location
	//	std::cout << "stop collision " << nsteps << " : current point " << point.xc() << " " << point.yc() << " " << point.zc() << std::endl;
      }
      // reset system, continue
      running_time = 0.0;
      previous = point;
      target_mass = pick_target(weight, which);
      mumass_eV = 1.0e9 * e_mass * target_mass / (e_mass + target_mass); // kinematics only
    }
    if (time_sum>=1.0e-5) { // 10 mus, particle got stuck, roughly 10^7 collisions
      analytic = true; // Stop
      std::cout << "STUCK: time = " << time_sum << std::endl;
      std::cout << "STUCK: place= " << point.xc() << " " << point.yc() << std::endl;
    }

  }
  // one charge done
  return false;
}


bool Ctransport::run(Electrode* electrode) {
  bool flag = false;
  // First, prepare electrode object for transport
  electrode->initfields(); // ready to transport

  unsigned int nthreads = std::thread::hardware_concurrency();
  if (nthreads>4) nthreads = 4; // limit max CPU number
  std::vector<std::future<bool> > results; 
  thread_pool* pool = new thread_pool(nthreads); // task pool

  charge_t q;
  int counter = 0;

  while (!charges.empty()) { // stop when refilling stopped
    //    std::cout << "from threads, charge basket size = " << charges.size() << std::endl;

    // empty charges and store tasks in blocks of nthreads
    for (int n=0;n<nthreads && !charges.empty();n++) { // drain charges basket
      q = charges.front(); // get front element of std::list
      results.push_back(pool->async(std::function<bool(Electrode*, charge_t)>(std::bind(&Ctransport::taskfunction, this, std::placeholders::_1, std::placeholders::_2)), electrode, q)); // tasks

      charges.pop_front(); // remove first charge from list
      counter++; // counts tasks/electrons launched
    }

    // drain task pool
    for (std::future<bool>& status : results){ 
      if (status.get()) // wait for completion before the next round
	flag = true;
      else
	continue;
    }
    // all tasks from pool finished - clear it. Next batch of charges in pool.
    results.clear();
  }

  //  std::cout << "from threads, total task counter = " << counter << std::endl;
  // charge loop finished
  delete pool;
  return flag;
}


void Ctransport::book_charge(charge_t q) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  // avalanche limit - hard cut on number of charges
  if (charges.size() > 1000) return;
  charges.push_back(q); // total charge list to be filled/drained in threads
  // avalanche feedback
  //  if (!(charges.size() % 500)) std::cout << "charges booked " << charges.size() << std::endl;
  return;
}


void Ctransport::book_time(double tt) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  times.push_back(tt); // time sum recorded
  return;
}


void Ctransport::book_place(Point3 loc) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  places.push_back(loc); // time sum recorded
  return;
}


int Ctransport::findBin(double en) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  // caution on energy
  if (en<0.0) en = 0.0;
  if (en>=40.0) return (int)energybins.size()-1; // final entry

  std::vector<double>::iterator low;
  low = std::lower_bound(energybins.begin(), energybins.end(), en); // find index with lower bound
  int bin = (low - energybins.begin());
  return bin;
}


double Ctransport::time_update(double tau)
{
    return -tau*TMath::Log(rnd->Rndm());
}

TVector3 Ctransport::speed_update(int charge, Point3 dfield, double time)
{
  TVector3 Efield(charge*dfield.xc(), charge*dfield.yc(), charge*dfield.zc()); // in [V/m]
  double eoverm = 1.759e11; // Coulomb / kg
  TVector3 v = eoverm * Efield * time;
  return v;
}

TVector3 Ctransport::d_update(TVector3 v0, double time)
{
    TVector3 dstep = v0*time; // acceleration done already in speed_update
    return dstep;
}

TVector3 Ctransport::kin_factor2(TVector3 v0, double target_mass)
{
    TVector3 vel = v0;
    //    double theta0 = v0.Theta();
    double phi0 = v0.Phi();
    double energy, transfer;
    double c2 = 2.99792458e8*2.99792458e8; // c^2 [m/s]^2
    double e_mass = 0.511e-3; // [GeV/c^2]
    double reduced_mass = (4.0*target_mass*e_mass)/((target_mass + e_mass)*(target_mass + e_mass));
    double mumass_eV = 1.0e9 * e_mass * target_mass / (e_mass + target_mass); // [eV]

    energy = 0.5 * mumass_eV * v0.Mag2() / c2; // non-rel. energy in [eV]
    double azimuth = angle_function2(energy);
    // TVector3 has theta defined relative to +z axis
    double theta = TMath::Pi()/2.0;// in plane theta

    // double phi = 0.5*TMath::ASin((target_mass + e_mass)/target_mass*TMath::Sin(theta));
    // transfer = TMath::Sqrt((1.0 - reduced_mass * TMath::Cos(phi)*TMath::Cos(phi)));

    vel.SetTheta(theta);
    vel.SetPhi(azimuth+phi0); // relative to previous
    return vel;
}

double Ctransport::angle_function2(double energy)
{
  // needs elastic scattering angular distribution in x,y plane
  return TMath::TwoPi()*rnd->Rndm();
}


// read and prepare the cross sections from file
void Ctransport::readCS(std::string csname) {
  TFile ff(csname.data(),"read");
  TNtupleD* nt = (TNtupleD*)ff.Get("cs");
  int nentries = nt->Draw("csel:csinel:energy","weight>0.5","goff"); // for helium cs histogram

  // set up x-axis for histograms
  double* en = nt->GetV3();

  // set up helium cs
  double* csel = nt->GetV1();
  double* csinel = nt->GetV2();
  for (int i=0;i<nentries;i++) {
    energybins.push_back(en[i]); // sorted already
    HeCSel.push_back(csel[i] * 1.e-4); // convert to SI [m^2]
    HeCSinel.push_back(csinel[i] * 1.e-4);
  }
  std::cout << "In CTransport: helium cross sections from file: " << nentries << std::endl;

  // set up ethanol cs
  nentries = nt->Draw("csel:csinel","weight>0.02 && weight<0.5","goff"); // for ethanol cs histogram
  csel = nt->GetV1();
  csinel = nt->GetV2();
  for (int i=0;i<nentries;i++) {
    EthCSel.push_back(csel[i] * 1.e-4);
    EthCSinel.push_back(csinel[i] * 1.e-4);
  }
  std::cout << "In CTransport: ethanol cross sections from file: " << nentries << std::endl;

  // set up argon cs
  nentries = nt->Draw("csel:csinel","weight<0.02","goff"); // for argon cs histogram
  csel = nt->GetV1();
  csinel = nt->GetV2();
  for (int i=0;i<nentries;i++) {
    ArCSel.push_back(csel[i] * 1.e-4);
    ArCSinel.push_back(csinel[i] * 1.e-4);
  }
  ff.Close();
  std::cout << "In CTransport: argon cross sections from file: " << nentries << std::endl;
}


// cross sections retrieve from MagBoltz output file.
double Ctransport::cross_section(double energy, int which, int &inel_flag)
{
  // retrieve from MagBoltz output file.
  double cs_el, cs_inel,ratio;
  int bin = findBin(energy);
  
  inel_flag = 0; // default case

  // shouldn't happen but does, it seems
  if (bin>=(int)energybins.size()) 
    bin = (int)energybins.size()-1; // final entry
  if (energy >= 40.0) energy = 39.9; // max range of energybins from MagBoltz
 
  // debug
  // if (energy>10.0) 
  //   std::cout << "Inside cross section with en " << energy << " gas which " << which << " and energy bin " << bin << std::endl;

  switch(which) {
  case 0: { // Helium
    cs_el = HeCSel.at(bin);
    cs_inel = HeCSinel.at(bin);
    if (cs_inel>0.0) { // decision elastic to ionization
      ratio = cs_inel / (cs_el+cs_inel);
      if (rnd->Rndm() < ratio) {
	inel_flag = 1; // ionization electron
	//	std::cout << "helium inelastic collision at energy " << energy << std::endl;
	return 0.0; // stops further transport after inel
      }
      else
	return cs_el;
    }
    else 
      return cs_el;    
  }
  case 1: { // Ethanol
    cs_el = EthCSel.at(bin);
    cs_inel = EthCSinel.at(bin);
    if (cs_inel>0.0) { // decision elastic to ionization
      ratio = cs_inel / (cs_el+cs_inel);
      if (rnd->Rndm() < ratio) {
	inel_flag = 1; // ionization electron
	//	std::cout << "ethanol inelastic collision at energy " << energy << std::endl;
	return 0.0; // stops further transport after inel
      }
      else
	return cs_el;
    }
    else
      return cs_el;
  }
  case 2: { // Argon
    cs_el = ArCSel.at(bin);
    cs_inel = ArCSinel.at(bin);
    if (cs_inel>0.0) { // decision elastic to ionization
      ratio = cs_inel / (cs_el+cs_inel);
      if (rnd->Rndm() < ratio) {
	inel_flag = 1; // ionization electron
	//	std::cout << "argon inelastic collision at energy " << energy << std::endl;
	return 0.0; // stops further transport after inel
      }
      else
	return cs_el;
    }
    else
      return cs_el;
  }
  default: return 0.0; // not needed but stops further transport after inel
  }
}
