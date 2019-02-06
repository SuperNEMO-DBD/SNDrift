#include <iostream>
// ROOT items
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoTorus.h"
#include "TGeoVolume.h"
#include "TGeoCompositeShape.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMath.h"

int main() {
  void trackervol();

  // simply call geometry function
  trackervol();
  return 1;
}


void trackervol()
{
  TGeoManager* geom = new TGeoManager("sndriftvol","tracker geometry");

  TGeoMaterial* mat = new TGeoMaterial("Vacuum", 0, 0, 1.e-20);
  TGeoMaterial* tgas = new TGeoMaterial("He", 4, 2, 1.e-20);
  TGeoMaterial *matMet = new TGeoMaterial("Al", 26, 13, 2.7);

  TGeoMedium* med = new TGeoMedium("Vacuum",1,mat);
  TGeoMedium* wiremed = new TGeoMedium("Al",2,matMet);
  TGeoMedium* tgasmed = new TGeoMedium("He",3,tgas);

  // units [cm]
  // note hard-coded names for volumes is required
  double xhalf = 8.0; // tracker box from Comsol
  double yhalf = 21.7; // tracker box from Comsol
  double halfheight = 2.0; // z is not relevant here but make some space

  double wirerad = 0.0025; // 25 mum wire radius as in Comsol
  double anoderad = 0.002; // 20 mum wire radius as in Comsol

  // geo translation constants, wire centre to centre
  double fw2fw = 0.911; // field wire spacing 9.11 mm
  double fwcorner = 1.289; // field wire corner 12.89 mm
  double a2fw  = 2.2; // anode spacing to field wire 22 mm

  // make shape components
  TGeoTube *fwtub  = new TGeoTube("C",0,wirerad,halfheight);
  TGeoTube *atub  = new TGeoTube("A",0,anoderad,halfheight);

  TGeoVolume* world = geom->MakeBox("World",med,4*xhalf+0.1,4*yhalf+0.1,halfheight+0.1); // larger than rest
  TGeoVolume* comsol = geom->MakeBox("Comsol",tgasmed,xhalf,yhalf,halfheight); // field volume

  TGeoVolume *fwire = new TGeoVolume("FWire", fwtub, wiremed);
  TGeoVolume *awire = new TGeoVolume("AWire", atub, wiremed);

  geom->SetTopVolume(world);

  // place all wires into tracker chamber
  // unit cell first, open right side and bottom
  double basey = yhalf - 0.7; // negative y like in comsol [cm]
  double basexmid = -xhalf + 2.689; // left top corner, build direction down
  double basexleft = -xhalf + 1.4; // left mid corner, build direction down

  for (int row=0; row<9; row++) { // 9 rows total
    for (int i=0; i<3; i++) { // 3 columns total
      for (int j=0; j<3; j++) { // layer 1
	comsol->AddNode(fwire, j + i*6 + row*21, new TGeoTranslation(basexmid + j*fw2fw + i*2*a2fw, basey - row*2*a2fw, 0.0));
      }
      // left column
      comsol->AddNode(fwire, 3 + i*6 + row*21, new TGeoTranslation(basexleft + i*2*a2fw, basey - fwcorner - row*2*a2fw, 0.0));
      comsol->AddNode(fwire, 4 + i*6 + row*21, new TGeoTranslation(basexleft + i*2*a2fw, basey - fwcorner - fw2fw - row*2*a2fw, 0.0));
      // anode
      comsol->AddNode(awire, i + row*3, new TGeoTranslation(basexleft + a2fw + i*2*a2fw, basey - fwcorner - fw2fw - row*2*a2fw, 0.0));
      // left column
      comsol->AddNode(fwire, 5 + i*6 + row*21, new TGeoTranslation(basexleft + i*2*a2fw, basey - fwcorner - 2*fw2fw - row*2*a2fw, 0.0));
    }
  }

  // right and bottom perimeter wires
  for (int i=0; i<3; i++) { // 3 columns total
    for (int j=0; j<3; j++) { // final layer, bottom
      comsol->AddNode(fwire, j + 194 + i*3, new TGeoTranslation(basexmid + j*fw2fw + i*2*a2fw, basey - 2*fwcorner - 2*fw2fw - 16*a2fw, 0.0));
    }
  }
  for (int row=0; row<9; row++) { // 9 rows total
    for (int j=0; j<3; j++) { // final layer, right
      comsol->AddNode(fwire, 18 + j + row*21, new TGeoTranslation(basexleft + 6*a2fw, basey - fwcorner - j*fw2fw - row*2*a2fw, 0.0));
    }
  }

  // tracker gas volume placed in world as in Comsol, bottom right quarter
  world->AddNode(comsol,1, new TGeoTranslation(xhalf, -yhalf, 0.0));

  world->SetLineColor(1);
  comsol->SetLineColor(2);
  comsol->SetVisibility(kTRUE);
  fwire->SetLineColor(2);
  awire->SetLineColor(2);
  geom->SetTopVisible();
  geom->SetVisOption(1);
  geom->CloseGeometry();
  // *** GEOMETRY closed ***

  geom->Export("trackergeom.gdml");

  // numerical volume position check
  double xcurrent, ycurrent;
  double xv =  3.4; // xstart
  double yv = -2.9; // row
  // double xv =  3.56486; // xstart
  // double yv = -2.89117; // row
  double zv = 0.0;
  //  step through geom
  for (int n=0;n<210;n++) {
    xcurrent = xv + n*0.001;
    ycurrent = yv;
    geom->SetCurrentPoint(xcurrent, ycurrent, zv);
    TGeoNode* nd = geom->FindNode();
    TGeoVolume* vol = nd->GetVolume();
    
    TString region(vol->GetName()); // changed from GetNumber()
    std::cout << "Geometry model: region = " << region << std::endl;
    std::cout << "Geometry model: coords: " << xcurrent << " " << ycurrent << " " << zv << std::endl;
  }
  //  world->Draw();

}


