#ifndef SNDRIFT_FIELDS_HH
#define SNDRIFT_FIELDS_HH

#include <vector>

// ROOT includes
#include "TKDTree.h"
#include "TString.h"

//local
#include "utils.hh"
#include "geomodel.hh"


//***********************************
// Field map classes
//***********************************
class ComsolFields {
 private:
  // scaling the weighting field for two electrodes
  double bias;
  // which ROOT file to read the weighting field
  TString fname;
  std::vector<Point3> coords;
  std::vector<Point3> dmap;

 protected:

 public:
  // Constructor
  ComsolFields(TString fname); // from file
  
  // Default destructor
  ~ComsolFields() {;}

  // Methods
  void read_fields();
  void setBias(double b) {bias = b;};
  std::vector<Point3> positions() {return coords;}
  std::vector<Point3> driftmap() {return dmap;}
};


class Fields {
 private:
  // pointer to geometry for asking
  GeometryModel* gm;
  // container for field coordinates here
  TKDTreeID* coordinates;
  // storage container
  double* allx;
  double* ally;
  double* alldx;
  double* alldy;

 protected:
  void prepare_fields(ComsolFields* fem);
  Point3 getFieldValue(Point3 p, bool& analytic);  

 public:
  // Constructor
  Fields(ComsolFields* fem, GeometryModel* gm); // from file
  
  // Default destructor
  ~Fields();

  // Methods
  // return field values in [V/m]
  Point3 getDriftField(Point3 p, bool& analytic);
};
#endif
