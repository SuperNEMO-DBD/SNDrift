// us
#include "utils.hh"

// general
#include "math.h"

//************************
// Euclid member functions
//************************

//*******
// Point3
//*******
Point3::Point3() {
  xcoord = 0.0;
  ycoord = 0.0;
  zcoord = 0.0;
}

Point3::Point3(double x, double y, double z) {
  xcoord = x;
  ycoord = y;
  zcoord = z;
}

void Point3::Set(double x, double y, double z) {
  xcoord = x;
  ycoord = y;
  zcoord = z;
}




//**************
// Line Segment3
//**************
LineSegment3::LineSegment3() {
  point_A = Point3();
  point_B = Point3();
  vec = TVector3(0.0,0.0,0.0);
}

LineSegment3::LineSegment3(Point3 a, Point3 b) {
  point_A = a;
  point_B = b;
  vec = TVector3(b.xc() - a.xc(), b.yc() - a.yc(), b.zc() - a.zc());
}


void LineSegment3::Set(Point3 a, Point3 b) {
  point_A = a;
  point_B = b;
  vec = TVector3(b.xc() - a.xc(), b.yc() - a.yc(), b.zc() - a.zc());
}


LineSegment3 LineSegment3::connect_point(Point3 p) {
  double d = vec.Mag2();
  if (d==0.0) {
    return LineSegment3();
  }
  double u = ((p.xc() - point_A.xc())*vec.X() + (p.yc() - point_A.yc())*vec.Y() + (p.zc() - point_A.zc())*vec.Z()) / d;
  Point3 pline(point_A.xc()+u*vec.X(),point_A.yc()+u*vec.Y(),point_A.zc()+u*vec.Z());
  LineSegment3 ls = LineSegment3(p,pline);
  return ls; // point to line point
}


LineSegment3 LineSegment3::connect_line(LineSegment3 ll) {
  TVector3 vll = ll.vector();
  Point3 pll = ll.start_point();
  TVector3 p13(pll.xc()-point_A.xc(),pll.yc()-point_A.yc(),pll.zc()-point_A.zc());
  double d1343 = p13.Dot(vec);
  double d4321 = vec.Dot(vll);
  double d1321 = p13.Dot(vll);
  double d4343 = vec.Mag2();
  double denom = vll.Mag2() * d4343 - d4321*d4321;
  if (denom==0.0) {
    Point3 pzero(0.0,0.0,0.0);
    return LineSegment3(pzero,pzero);
  }
  double ua = (d1343 * d4321 - d1321 * d4343) / denom;
  double ub = (d1343 + d4321 * ua) / d4343;
  Point3 p1(ll.start_point().xc() + ua * ll.vector().X(),ll.start_point().yc() + ua * ll.vector().Y(),ll.start_point().zc() + ua * ll.vector().Z());
  Point3 p2(point_A.xc() + ub * vec.X(),point_A.yc() + ub * vec.Y(),point_A.zc() + ub * vec.Z());
  return LineSegment3(p1,p2);
}


//**************
// Bounding Box
//**************
BBox::BBox() {
  point_A = Point3();
  xaxis = TVector3(0.0,0.0,0.0);
  yaxis = TVector3(0.0,0.0,0.0);
  zaxis = TVector3(0.0,0.0,0.0);
}


BBox::BBox(Point3 a, Point3 b, Point3 c, Point3 d) {
  point_A = a;
  TVector3 o(point_A.xc(),point_A.yc(),point_A.zc());
  TVector3 tox(b.xc(),b.yc(),b.zc());
  TVector3 toy(c.xc(),c.yc(),c.zc());
  TVector3 toz(d.xc(),d.yc(),d.zc());
  xaxis = tox - o;
  yaxis = toy - o;
  zaxis = toz - o;
}


void BBox::Set(Point3 a, Point3 b, Point3 c, Point3 d) {
  point_A = a;
  TVector3 o(point_A.xc(),point_A.yc(),point_A.zc());
  TVector3 tox(b.xc(),b.yc(),b.zc());
  TVector3 toy(c.xc(),c.yc(),c.zc());
  TVector3 toz(d.xc(),d.yc(),d.zc());
  xaxis = tox - o;
  yaxis = toy - o;
  zaxis = toz - o;
}


bool BBox::IsInside(Point3 p) {
  bool inside = 0;
  TVector3 o(point_A.xc(),point_A.yc(),point_A.zc());
  TVector3 check(p.xc(),p.yc(),p.zc());
  TVector3 topoint = check - o;
  double cpx = topoint.Dot(xaxis); // dot products
  double cpy = topoint.Dot(yaxis);
  double cpz = topoint.Dot(zaxis);

  if (cpx>0.0 && cpx<=xaxis.Mag2()) {
    if (cpy>0.0 && cpy<=yaxis.Mag2()) {
      if (cpz>0.0 && cpz<=zaxis.Mag2()) {
	inside = 1;
      }
    }
  }
  return inside;
}


