#ifndef SNDRIFT_CTRANSPORT_HH
#define SNDRIFT_CTRANSPORT_HH

#include <list>
#include <vector>
#include <string>
#include <mutex>

// ROOT
#include "TRandom3.h"
#include "TVector3.h"

//local
#include "utils.hh"
#include "electrode.hh"

//***********************************
// Charge signal class
// to be used as an interface
// to the algorithm.
//***********************************
class Ctransport {
 private:
  double density;
  std::list<charge_t> charges;
  std::vector<double> times;
  std::vector<Point3> places;
  std::mutex mtx;
  TRandom3* rnd;
  std::vector<double> energybins;
  std::vector<double> HeCSel; // three gas cross section containers
  std::vector<double> EthCSel;
  std::vector<double> ArCSel;
  std::vector<double> HeCSinel; // three gas cross section containers
  std::vector<double> EthCSinel;
  std::vector<double> ArCSinel;

  // used by task function
  void book_charge(charge_t q);
  void book_time(double tt);
  void book_place(Point3 loc);
  void readCS(std::string csname);
  int  findBin(double en);
  double time_update(double tau);
  double angle_function2(double energy);
  double cross_section(double energy, int which, int &inel_flag);
  double pick_target(std::vector<double> weight, int& which);
  TVector3 speed_update(int charge, Point3 dfield, double time);
  TVector3 d_update(TVector3 v0, double time);
  TVector3 kin_factor2(TVector3 v0, double tm);


 protected:
  bool run(Electrode* electrode);
  bool taskfunction(Electrode* electrode, charge_t q);

 public:
  // Constructor
  Ctransport(std::string fname, int seed);
  
  // Default destructor
  ~Ctransport();

  // Methods
  // preparation, required input from main()
  // otherwise no transport possible
  // work on this electrode id with charges
  void ctransport(Electrode* electrode, std::list<charge_t> q); 
  std::vector<double> getDriftTimes() {return times;}
  std::vector<Point3> getLocations() {return places;}
  double getDensity() {return density;}
  void setDensity(double d) {density = d;};
};
#endif
