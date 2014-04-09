#ifndef __LIGHTSIMULATOR_H__
#define __LIGHTSIMULATOR_H__
#include "defs.h"
#include "EnergyDeposit.h"
#include "SimulationOutput.h"

using namespace TMath;
using namespace std;

class LightSimulator{
  
public:

  LightSimulator(EnergyDeposit deposit_, Int_t ntoys_);
  ~LightSimulator(){
    if (output) delete output;
  };

  EnergyDeposit deposit;
  Int_t nrays;
  Int_t ntoys;

  SimulationParameters pars;

  Double_t limit_angle = ASin(1./pars.material_refr_index);

  bool debug = false;
  void SetDebug(bool debug_);

  bool isinsideboundaries();

  inline void swap(Double_t &x, Double_t &y);
  inline Double_t round(Double_t a);
  
  bool propray(const LightRay &lr, AxisDirection axdir, Int_t index, Int_t &other_coordinate, MatchObject &m);
  Int_t findchannel(Int_t x, Int_t y);
  
  SimulationOutput *output = NULL;

  SimulationOutput* Run();
  void Normalize(TH1F *h);

  vector<GoodMatchObject>  convert_matches(const vector<MatchObject> &a);

};


#endif
