// *********************************
// SNDrift: Monte-Carlo drift time from point
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
  std::cout << "Monte-Carlo scan command line option(s) help" << std::endl;
  std::cout << "\t -x , --xstart <x-coordinate start [cm]>" << std::endl;
  std::cout << "\t -y , --ystart <y-coordinate start [cm]>" << std::endl;
  std::cout << "\t -b , --bias <Anode bias in Volt>" << std::endl;
  std::cout << "\t -s , --seed <random number seed offset>" << std::endl;
  std::cout << "\t -n , --nsim <number of Monte Carlo simulations>" << std::endl;
  std::cout << "\t -d , --dataDir <FULL PATH Directory to data file>" << std::endl;
  std::cout << "\t -o , --outputFile <FULL PATH ROOT FILENAME>" << std::endl;
}



int main(int argc, char** argv) {
  // function declare
  void signal_calculation(int seed, int nsim, double bias, double xstart, double ystart, std::string fname);

  int seed, nsim;
  double bias, xs, ys;
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
  ops >> GetOpt::Option('s', "seed", seed, 0);
  ops >> GetOpt::Option('n', "nsim", nsim, 10);
  ops >> GetOpt::Option('d', dataDirName, "data/");
  ops >> GetOpt::Option('o', outputFileName, "");

  if (outputFileName=="")
    outputFileName = "drifttimes.root";

  //run the code
  signal_calculation(seed, nsim, bias, xs, ys, outputFileName);
  
  return 0;
}



void signal_calculation(int seed, int nsim, double bias, double xstart, double ystart, std::string fname) {

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
  std::string gfname = dataDirName+"trackergeom.gdml";
  GeometryModel* gmodel = new GeometryModel(gfname);

  //----------------------------------------------------------
  // FEM fields from file
  // hard-coded field map files
  std::string femname = dataDirName+"sntracker_driftField.root";
  ComsolFields* fem = new ComsolFields(femname.data());
  fem->setBias(bias);

  //----------------------------------------------------------
  // Transport
  std::string fn = dataDirName+"trackergasCS.root";
  Ctransport* ctr = new Ctransport(fn, seed);
  // setting up

  //----------------------------------------------------------
  // transport start
  //----------------------------------------------------------
  fem->read_fields();
  Electrode* anode = new Electrode(fem, gmodel);

  std::vector<double> timestore;
  for (int nn=0; nn<nsim; nn++) { // Monte Carlo loop
    ctr->ctransport(anode, hits);
    std::vector<double> dts = ctr->getDriftTimes();
    timestore.push_back(*std::max_element(dts.begin(), dts.end()));
    // some feedback
    std::cout << "drift time: " << *std::max_element(dts.begin(), dts.end()) << std::endl;
  }

  //----------------------------------------------------------
  // to storage
  //----------------------------------------------------------
  // metainfo
  TParameter<double> xpar("xstart",xstart);
  TParameter<double> ypar("ystart",ystart);
  TParameter<double> bpar("bias",bias);

  // file
  TFile ff(fname.c_str(),"RECREATE");
  TNtupleD* ntstore = new TNtupleD("drift_results","Stopping times","sx:sy:dtime");
  for (double tt : timestore)
    ntstore->Fill(xstart, ystart, tt);

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


