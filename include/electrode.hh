#ifndef SNDRIFT_WIRE_HH
#define SNDRIFT_WIRE_HH

#include <mutex>

//local
#include "utils.hh"
#include "geomodel.hh"
#include "fields.hh"


//***********************************
// Charge signal class
// to be used as an interface
// to the algorithm.
//***********************************
class Electrode {
 private:
  // pointer to geometry for constructing fields
  GeometryModel* gm;
  // pointer to comsol fields for constructing fields
  ComsolFields* femfields;

  Fields* field; // specific for each electrode, constructed at creation
  std::mutex mtx;

 protected:

 public:
  // Constructor
  // using geometry data
  Electrode(ComsolFields* fem, GeometryModel* g);
  
  // Default destructor
  ~Electrode();

  // access
  void initfields(); // out of constructor - takes time.

  Point3 getFieldValue(bool& analytic, Point3 p);
};
#endif
