#ifndef __SIMULATIONOUPUT_H__
#define __SIMULATIONOUPUT_H__

#include "defs.h"
#include "EnergyDeposit.h"
#include "TParameter.h"
using namespace TMath;
using namespace std;

class SimulationOutput{

 public:
  TH1F* chamfer_photons[4];
  TH1F* chamfer_pulseshape[4];
  TH1F* reflections_all;
  TH1F* reflections_total;
  TH1F* reflections_paper;
  TH1F* optical_path;
  TH1F* time_of_arrival;
  EnergyDeposit deposit;

  SimulationOutput(){}; 
  SimulationOutput(EnergyDeposit deposit_);
  ~SimulationOutput(){
    for (int i=0; i<4; i++) if (chamfer_photons[i]) delete chamfer_photons[i];
    for (int i=0; i<4; i++) if (chamfer_pulseshape[i]) delete chamfer_pulseshape[i];
    if (reflections_all) delete reflections_all;
    if (reflections_total) delete reflections_total;
    if (reflections_paper) delete reflections_paper;
    if (optical_path) delete optical_path;
    if (time_of_arrival) delete time_of_arrival;    
};
  void Write();

  SimulationParameters pars;

};

#endif
