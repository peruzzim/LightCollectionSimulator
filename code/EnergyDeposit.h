#ifndef __ENERGYDEPOSIT_H__
#define __ENERGYDEPOSIT_H__

#include "defs.h"
using namespace TMath;
using namespace std;

class EnergyDeposit{
 public:
  CaloPoint position;
  Double_t energy;
  Double_t time;
  Int_t totalrays;

  SimulationParameters pars;

 EnergyDeposit(CaloPoint pos = TVector3(0,0,0), Double_t E=0, Double_t t=0) : position(pos),energy(E),time(t){
    totalrays = energy*pars.photons_per_mev;
  };
  
  CaloPoint getpositionlabframe();
  void setpositionlabframe(TVector3 labframeposition);

};
#endif
