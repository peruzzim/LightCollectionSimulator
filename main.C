#include "code/defs.h"
#include "code/EnergyDeposit.h"
#include "code/SimulationOutput.h"
#include "code/LightSimulator.h"
#include <iostream>
#include <unistd.h>

using namespace std;

int main(int argc, char* argv[]) {

  TFile *f = new TFile ("scan.root","recreate");
  f->cd();

  SimulationParameters pars;
  float x=0;
  float y=0;
  float z=0;
  float energy=1;
  float t=0;
  int toys=1;

  char ch;
  while ((ch = getopt(argc, argv, "E:x:y:z:t:n:")) != -1 ) {
    switch (ch) {
    case 'E': energy = atof(optarg); break;
    case 'x': x = atof(optarg); break;
    case 'y': y = atof(optarg); break;
    case 'z': z = atof(optarg); break;
    case 't': t = atof(optarg); break;
    case 'n': toys = atoi(optarg); break;
    default:
      cout << "Problem in option parsing" << endl;
      cout << "Usage: ./simulation.run -x xcoord -y ycoord -z zcoord (all w.r.t. center of the crystal) -E energy(MeV) -t time(ps) -n toys" << endl;
      return -1;
    }
  }

  argc -= optind;
  argv += optind;


  EnergyDeposit myEnergyDeposit(TVector3(0,0,0),energy,t);
  myEnergyDeposit.setpositionlabframe(TVector3(x,y,z));
  cout << Form("Runtime parameters (x,y,z)=(%f,%f,%f), E=%f, t=%f, n=%d",x,y,z,energy,t,toys) << endl;
  LightSimulator *a = new LightSimulator(myEnergyDeposit,toys);
  SimulationOutput *output = a->Run();
  output->Write();

  f->Close();

  system("root -l scan.root plot.C");

  return 0;
  
};


