// us
#include "electrode.hh"

//****************
// Electrode Model
//****************
// Constructors
Electrode::Electrode(ComsolFields* fem, GeometryModel* g) {
  gm = g;
  femfields = fem;
  active = false;
  field = 0;
}


// default Destructor
Electrode::~Electrode() {
  if (field) delete field;
}


void Electrode::initfields() {
  active = true;
  field = new Fields(femfields, gm); // create from file + geometry info
}


Point3 Electrode::getFieldValue(bool& analytic, Point3 p) {
  std::lock_guard<std::mutex> lck (mtx); // protect thread access
  Point3 triplet;
  
  triplet = field->getDriftField(p, analytic);

  return triplet;
}


