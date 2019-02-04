#ifndef SNDRIFT_UTILS_HH
#define SNDRIFT_UTILS_HH

#include <vector>

// ROOT includes
#include "TVector3.h"
#include "TH1F.h"

//*****************
// Geometry helpers
//*****************
class Point3 {
 private:
  double xcoord;
  double ycoord;
  double zcoord;
  
 public:
  Point3();
  Point3(double x, double y, double z);
  ~Point3() {;}
  
  void Set(double x, double y, double z);
  double xc() {return xcoord;}
  double yc() {return ycoord;}
  double zc() {return zcoord;}
};


class LineSegment3 {
 private:
  TVector3 vec;
  Point3 point_A;
  Point3 point_B;
    
 public:
  LineSegment3();
  LineSegment3(Point3 a, Point3 b);
  ~LineSegment3() {;}
  
  TVector3 vector() {return vec;}
  Point3 start_point() {return point_A;}
  Point3 end_point() {return point_B;}
  
  double length() {return vec.Mag();}
  LineSegment3 connect_point(Point3 p); // find connection
  LineSegment3 connect_line(LineSegment3 ll); // find connection
  void Set(Point3 a, Point3 b);
};


class BBox {
 private:
  Point3 point_A;
  TVector3 xaxis;
  TVector3 yaxis;
  TVector3 zaxis;
    
protected:

 public:
  BBox();
  BBox(Point3 a, Point3 b, Point3 c, Point3 d);
  ~BBox() {;}
  
  void Set(Point3 a, Point3 b, Point3 c, Point3 d);
  bool IsInside(Point3 p);
};


// Helper structures
struct charge_t {
  Point3 location;
  int charge;
  int chargeID; // distinguish e- (1) and gamma (0)
};


typedef std::vector<Point3> path_t;

#endif
