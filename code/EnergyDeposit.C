#ifndef __ENERGYDEPOSIT_C__
#define __ENERGYDEPOSIT_C__

#include "EnergyDeposit.h"
using namespace TMath;
using namespace std;

CaloPoint EnergyDeposit::getpositionlabframe(){
    TRotation r;
    r.RotateZ(Pi()/4);
    TVector3 a=position;
    a = r*a;
    return a;
  };
void EnergyDeposit::setpositionlabframe(TVector3 labframeposition){
  TRotation r;
  r.RotateZ(-Pi()/4);
  position = r*labframeposition;
};

#endif
