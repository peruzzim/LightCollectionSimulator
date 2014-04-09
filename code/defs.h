#ifndef __DEFS_H__
#define __DEFS_H__

#include <TROOT.h>
#include "TStyle.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TRandom3.h"
#include <iostream>
#include <assert.h>

using namespace TMath;
using namespace std;

typedef TVector3 CaloPoint;

typedef enum {
  kNE=0,
  kSE,
  kSW,
  kNW
} XYDirection;

typedef enum {
  kNS=0,
  kEW,
} AxisDirection;


// PARAMETERS (all dimensions in millimeters) ////////////////////////////////

class SimulationParameters{

 public:

  Double_t xtalsize = 24;
  Double_t xtalheight = 10;
  Double_t chamfersize = 2.5;
  Double_t eff_reflection_total = 1;
  Double_t eff_reflection_paper = 0.95;
  Double_t eff_lightcoll = 1;
  Double_t light_att_length = 631; // remember: in mm // = transmission 80% after 141 mm
  Double_t limit_arrival_time = 5e3; // in ps
  Double_t material_refr_index = 1.62;
  Double_t photons_per_mev = 1000;
  Double_t scintillation_typ_timescale = 10; // in ps
  Double_t efficiency_fiber_collection = 1;
  Double_t efficiency_WLS_mechanism = 0.01;
  Double_t PMT_qe = 0.25;
  Double_t PMT_gain = 1e7;
  Double_t electron_charge = 1.6e-7; // in pC
  Double_t conversion_chamferphotons_to_charge = efficiency_fiber_collection*efficiency_fiber_collection*PMT_qe*PMT_gain*electron_charge;
  Double_t speed_of_light = 0.3/material_refr_index; // in mm/ps
  Int_t    max_distance_unit_propagation = 1000;

};



class CalculateDirection{
 public:
  static bool isNbound(XYDirection a) {return (a==kNE || a==kNW);}
  static bool isSbound(XYDirection a) {return (a==kSE || a==kSW);}
  static bool isWbound(XYDirection a) {return (a==kNW || a==kSW);}
  static bool isEbound(XYDirection a) {return (a==kNE || a==kSE);}

};




class LightRay{
public:
  CaloPoint origin;
  TVector3 dirvect;
  XYDirection xydir;

  LightRay(CaloPoint origin_, TVector3 dirvect_): origin(origin_), dirvect(dirvect_){
    Double_t phi = dirvect.Phi();
    if (phi>=0 && phi<Pi()/2) xydir=kNE;
    else if (phi>=Pi()/2) xydir=kNW;
    else if (phi<0 && phi>=-Pi()/2) xydir=kSE;
    else xydir=kSW;
  }

  void Print(){
    cout << "LightRay" << endl;
    origin.Print();
    dirvect.Print();
  };


};

class MatchObject{
public:
  Int_t index;
  AxisDirection axdir;
  XYDirection incomingdir;
  Int_t othercoordinate;
  Double_t position;
};

class GoodMatchObject{
public:
  Int_t x;
  Int_t y;
  Double_t positionx;
  Double_t positiony;
  void Set(Int_t x_, Int_t y_){ x=x_; y=y_;}
};

#endif
