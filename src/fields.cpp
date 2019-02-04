#include <iostream>

// us
#include "fields.hh"

// ROOT includes
#include "TFile.h"
#include "TNtupleD.h"
#include "TMath.h"

ComsolFields::ComsolFields(TString fn) {
  bias = 1.0;
  fname = fn;
  coords.clear();
  dmap.clear();
}


// come now as 2D data in x,y from comsol
void ComsolFields::read_fields() {
  Point3 p;
  Point3 pe;

  TFile* ffd = new TFile(fname.Data(),"read");
  TNtupleD* ntd = (TNtupleD*)ffd->Get("drift");
  int entries = ntd->GetEntries();

  double x,y;
  double wx,wy;

  ntd->SetBranchAddress("x",&x);
  ntd->SetBranchAddress("y",&y);

  ntd->SetBranchAddress("ex",&wx);
  ntd->SetBranchAddress("ey",&wy);


  // comsol y-coord becomes z-coord in geometry
  for (int i=0;i<entries;i++){
    ntd->GetEntry(i);
    p.Set(x*1.0e2, y*1.0e2, 0.0); // [m]->[cm]
    coords.push_back(p);
    pe.Set(bias*wx, bias*wy, 0.0); // electrode weighting to drift
    dmap.push_back(pe);
  }
  std::cout << "in Comsol Fields: read field entries " << entries << std::endl;
  // all done and in memory
  ffd->Close();
}


Fields::Fields(ComsolFields* fem, GeometryModel* g) {
  gm = g; // have access to geometry model

  allx = 0; // null ptr
  ally = 0;
  alldx = 0;
  alldy = 0;
  
  prepare_fields(fem);
}

Fields::~Fields() {
  if (allx) { // all set together
    delete [] allx ;
    delete [] ally ;
    delete [] alldx;
    delete [] alldy;
  }
  if (coordinates) delete coordinates;
}

void Fields::prepare_fields(ComsolFields* fem) {
  std::vector<Point3> cdata = fem->positions();
  int nentries = cdata.size();
  //  std::cout << "in Fields::prepare fields." << std::endl;

  coordinates = new TKDTreeID(nentries,2,1);
  
  allx = new double [nentries];
  ally = new double [nentries];

  for (int i=0;i<nentries;i++){
    // no transf needed, requests come as vectors
    allx[i] = cdata[i].xc();
    ally[i] = cdata[i].yc();
  }
  coordinates->SetData(0,allx);
  coordinates->SetData(1,ally);
  coordinates->Build();

  alldx = new double [nentries];
  alldy = new double [nentries];

  std::vector<Point3> ddata = fem->driftmap();
  for (int i=0;i<nentries;i++){
    alldx[i] = ddata[i].xc();
    alldy[i] = ddata[i].yc();
  }  
  // all done and in memory

}




Point3 Fields::getDriftField(Point3 p, bool& analytic) {
  Point3 triplet = getFieldValue(p, analytic);
  return triplet;
}



Point3 Fields::getFieldValue(Point3 p, bool& analytic) {

  // common routine to ask for field value

  // ask the geometry
  double xv = p.xc();
  double yv = p.yc();
  double zv = p.zc();

  int value = gm->whereami(xv,yv,zv);
  //  std::cout << "in Fields::answer to whereami: " << value << std::endl;
  Point3 triplet;
  
  if (value==1) { // comsol region
    double point[2];
    double dist[8]; // check on nearest 8 neighbours in grid
    int indx[8];
    
    TVector3 fieldvec;
    TVector3 sumvec;
    std::vector<TVector3> nnvec;
    
    double dsum = 0.0;
    point[0] = xv;  // relative to origin x
    point[1] = yv;  // relative to origin y

    //    std::cout << "in Fields: point coordinates " << xv << " " << yv << " " << zv << std::endl;
        
    coordinates->FindNearestNeighbors(point,8,indx,dist);
    for (int j=0;j<8;j++) {
      fieldvec.SetXYZ(alldx[indx[j]], alldy[indx[j]], 0.0);
      //      std::cout << "in Fields: nearest coords: " << allx[indx[j]] << " " << ally[indx[j]] << std::endl;
      //      std::cout << "in Fields: Drift field value: " << alldx[indx[j]] << " " << alldy[indx[j]] << std::endl;
      nnvec.push_back(fieldvec);
      dsum += dist[j];
    }

    double denom = 0.0;
    for (int j=0;j<8;j++) denom += (1.0-dist[j]/dsum);
    
    sumvec.SetXYZ(0.,0.,0.);
    for (int j=0;j<8;j++) {
      fieldvec = nnvec.at(j)*((1.0-dist[j]/dsum)/denom);
      sumvec += fieldvec;
    }
    // getting the proportions right between x and y field components 
    triplet.Set(sumvec.X(),sumvec.Y(),sumvec.Z());
    //    std::cout << "in Fields: point coordinates " << xv << " " << yv << " " << zv << std::endl;
    //    std::cout << "in Fields: average field value " << sumvec.X() << " " << sumvec.Y() << " " << sumvec.Z() << std::endl;
    
  }
  else { // any other region than Comsol like a wire or world.
    // outside anything relevant, stop transport.
    triplet.Set(-1.0,0.0,0.0);
    analytic = true; // trigger to stop
  }
  return triplet;
}
