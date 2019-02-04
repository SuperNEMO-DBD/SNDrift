#include "catch.hpp"

// us
#include "ctransport.hh"
#include "fields.hh"
#include "geomodel.hh"


int check_geometry(){
  // reach from testing directory
  const char* gfname = "../data/trackergeom.gdml";
  GeometryModel* gmodel = new GeometryModel(gfname);
  if (!gmodel) return 0;
  else return 1; // fine
}


int check_wire(){
  const char* gfname = "../data/trackergeom.gdml";
  GeometryModel* gmodel = new GeometryModel(gfname);
  int value =  gmodel->whereami(3.6, -2.9, 0.0); // should hit a wire
  return value; // should be -1 for wire
}


unsigned int check_fields(){
  // reach from testing directory
  ComsolFields* fem = new ComsolFields("../data/sntracker_driftField.root");
  fem->read_fields();
  return fem->positions().size(); // should be currently 258462
}


double check_readcs(){
  // reach from testing directory
  std::string fn = "../data/trackergasCS.root";
  Ctransport* ctr = new Ctransport(fn, 0);
  return ctr->getDensity(); // should be helium 0.1664
}


TEST_CASE( "Geometry in", "[sndrift][geo_in]" ) {
  REQUIRE( check_geometry() == 1 );
}

TEST_CASE( "Geometry wire", "[sndrift][wiretest]" ) {
  REQUIRE( check_wire() == -1 );
}

TEST_CASE( "Fields in", "[sndrift][fieldtest]" ) {
  REQUIRE( check_fields() == 258462 );
}

TEST_CASE( "CS in", "[sndrift][cstest]" ) {
  REQUIRE( check_readcs() == 0.1664 );
}
