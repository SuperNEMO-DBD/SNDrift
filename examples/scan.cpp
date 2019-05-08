// *********************************
// SNDrift: scan drift time from point
//**********************************

#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

// us
#include "ctransport.hh"
#include "electrode.hh"
#include "fields.hh"
#include "geomodel.hh"
#include "getopt_pp.h"
#include "utils.hh"

// ROOT
#include "TFile.h"
#include "TNtupleD.h"
#include "TParameter.h"

void showHelp() {
  std::cout << "collection scan command line option(s) help" << std::endl;
  std::cout << "\t -x , --xstart <x-coordinate start [cm]>" << std::endl;
  std::cout << "\t -y , --ystart <y-coordinate start [cm]>" << std::endl;
  std::cout << "\t -b , --bias <Anode bias in Volt>" << std::endl;
  std::cout << "\t -p , --pressure <tracker gas pressure [mbar]>" << std::endl;
  std::cout << "\t -s , --seed <random number seed offset>" << std::endl;
  std::cout << "\t -o , --outputFile <FULL PATH ROOT FILENAME>" << std::endl;
}



int main(int argc, char** argv) {
  // function declare
  void signal_calculation(int seed, double bias, double xstart, double ystart, double pr, std::string fname);

  int seed;
  double bias, xs, ys, pressure;
  std::string outputFileName;
  GetOpt::GetOpt_pp ops(argc, argv);

  // Check for help request
  if (ops >> GetOpt::OptionPresent('h', "help")){
    showHelp();
    return 0;
  }
  
  ops >> GetOpt::Option('x', "xstart", xs, 3.5);
  ops >> GetOpt::Option('y', "ystart", ys, -2.9);
  ops >> GetOpt::Option('b', "bias", bias, 1000.0);
  ops >> GetOpt::Option('p', "pressure", pressure, 1013.25);
  ops >> GetOpt::Option('s', "seed", seed, 0);
  ops >> GetOpt::Option('o', outputFileName, "");

  if (outputFileName=="")
    outputFileName = "avalanche.root";

  //run the code
  signal_calculation(seed, bias, xs, ys, pr, outputFileName);
  
  return 0;
}



void signal_calculation(int seed, double bias, double xstart, double ystart, double pr, std::string fname) {

  charge_t hit;
  Point3 loc(xstart, ystart, 0.0); // [cm] unit from root geometry
  int qq = -1; // [e]
  hit.location = loc;
  hit.charge = qq;
  hit.chargeID = 0;
  std::list<charge_t> hits;
  hits.push_back(hit); // let's have the one only

  //----------------------------------------------------------
  // Geometry
  const char* gfname = "data/trackergeom.gdml";
  GeometryModel* gmodel = new GeometryModel(gfname);

  //----------------------------------------------------------
  // FEM fields from file
  // hard-coded field map files
  ComsolFields* fem = new ComsolFields("data/sntracker_driftField.root");
  fem->setBias(bias);

  //----------------------------------------------------------
  // Transport
  std::string fn = "data/trackergasCS.root";
  Ctransport* ctr = new Ctransport(fn, seed);
  ctr->setDensity(0.1664 * pr / 1013.25); // [kg/m^3]  pr [mbar] / NTP (295K) helium gas density
  // setting up

  //----------------------------------------------------------
  // transport start
  //----------------------------------------------------------
  fem->read_fields();
  Electrode* anode = new Electrode(fem, gmodel);

  ctr->ctransport(anode, hits);
  // for (double tt : ctr->getDriftTimes())
  //   std::cout << "Drift time: " << tt << std::endl;
  // for (Point3 loc : ctr->getLocations())
  //   std::cout << "Places: " << loc.xc() << " " << loc.yc() << " " << loc.zc() << std::endl;

  //----------------------------------------------------------
  // to storage
  //----------------------------------------------------------
  // metainfo
  TParameter<double> xpar("xstart",xstart);
  TParameter<double> ypar("ystart",ystart);
  TParameter<double> bpar("bias",bias);

  // file
  TFile ff(fname.c_str(),"RECREATE");
  TNtupleD* ntstore = new TNtupleD("drift_results","Stopping locations and times","dtime:sx:sy:sz");
  std::vector<double> dts = ctr->getDriftTimes();
  std::vector<Point3> sloc = ctr->getLocations();
  std::cout << "charges written: " << sloc.size() << std::endl;
  std::cout << "Drift time = " << *std::max_element(dts.begin(), dts.end()) << std::endl;
  for (unsigned int i=0;i<dts.size();i++)
    ntstore->Fill(dts.at(i), sloc.at(i).xc(), sloc.at(i).yc(), sloc.at(i).zc());

  // store metainfo in both ntuples
  ntstore->GetUserInfo()->Add(&xpar);
  ntstore->GetUserInfo()->Add(&ypar);
  ntstore->GetUserInfo()->Add(&bpar);

  ntstore->Write();

  ff.Close();

  delete anode;
  delete ctr;
  delete fem;
  delete gmodel;

  return;
}


