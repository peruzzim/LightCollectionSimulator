#ifndef __LIGHTSIMULATOR_C__
#define __LIGHTSIMULATOR_C__
#include "LightSimulator.h"

using namespace TMath;
using namespace std;

LightSimulator::LightSimulator(EnergyDeposit deposit_, Int_t ntoys_):deposit(deposit_),ntoys(ntoys_){
  nrays = deposit.totalrays;
};

void LightSimulator::SetDebug(bool debug_){
  debug = debug_;
};

bool LightSimulator::isinsideboundaries(){
  float x = deposit.getpositionlabframe().x();
  float y = deposit.getpositionlabframe().y();
  if (fabs(x)>pars.xtalsize/2) return false;
  if (fabs(y)>pars.xtalsize/2) return false;
  if (y-pars.xtalsize/2>-(x-pars.xtalsize/2)-pars.chamfersize/sqrt(2)) return false;  
  if (y+pars.xtalsize/2<-(x+pars.xtalsize/2)+pars.chamfersize/sqrt(2)) return false;   
  if (y-pars.xtalsize/2>+(x+pars.xtalsize/2)-pars.chamfersize/sqrt(2)) return false;  
  if (y+pars.xtalsize/2<+(x-pars.xtalsize/2)+pars.chamfersize/sqrt(2)) return false; 
  return true;
};

inline void LightSimulator::swap(Double_t &x, Double_t &y){
  Double_t t=x;
  x=y;
  y=t;
};

inline Double_t LightSimulator::round(Double_t a){
  return floor(a+0.5);
};
  
bool LightSimulator::propray(const LightRay &lr, AxisDirection axdir, Int_t index, Int_t &other_coordinate, MatchObject &m){

  Double_t lineat = +pars.xtalsize/sqrt(2)-pars.chamfersize/2+index*pars.xtalsize*sqrt(2);    
  Double_t val1=lr.origin.y();
  Double_t val2=lr.origin.x();
  Double_t valb1=lr.dirvect.y();
  Double_t valb2=lr.dirvect.x();
    
  if (axdir==kNS){
    if (CalculateDirection::isSbound(lr.xydir)){
      lineat=-lineat;
    }
  }
  else {
    if (CalculateDirection::isWbound(lr.xydir)){
      lineat=-lineat;
    }
    swap(val1,val2);
    swap(valb1,valb2);
  }
    
  Double_t deltaA = lineat-val1;
  Double_t deltaB = deltaA/valb1*valb2;
  Double_t intersect = val2+deltaB;
  Double_t sigmas = intersect/(pars.xtalsize*sqrt(2));
  Double_t position = sigmas-Int_t(sigmas);
    
  if (fabs(position)<pars.chamfersize/2/(pars.xtalsize*sqrt(2))) {
    other_coordinate = round(sigmas);
    m.index=index;
    m.axdir=axdir;
    m.incomingdir=lr.xydir;
    m.othercoordinate=other_coordinate;
    m.position=position;
    if (debug) cout << "Match found " << "incdir/axdir " <<  m.incomingdir << m.axdir << " index " << m.index << " other" << m.othercoordinate << endl;
    return true;
  }
  else return false;
    
};

Int_t LightSimulator::findchannel(Int_t x, Int_t y){
  assert ((x+y)%2==0);
  Int_t rx = x%4;
  Int_t ry=y%4;
  if (rx<0) rx+=4;
  if (ry<0) ry+=4;
  if (rx==0) return (ry==0) ? 3 : 1;
  else if (rx==2) return (ry==0) ? 1 : 3;
  else if (rx==1) return (ry==1) ? 2 : 4;
  else if (rx==3) return (ry==1) ? 4 : 2;
  assert(false);
};
  

SimulationOutput* LightSimulator::Run(){
    
  TRandom3 randgen = TRandom3(0);

  if (debug) cout << ntoys*nrays << " light propagations to simulate" << endl;

  long counter_ray=0;

  if (output) delete output;
  output = new SimulationOutput(deposit);
    
  if (!isinsideboundaries()) return output;

  for (int nrun = 0; nrun<ntoys; nrun++){
            
    vector<Int_t> myphotons(4,0);

    for (int nray = 0; nray<nrays; nray++){

      counter_ray++;
      //	if (counter_ray%(ntoys*nrays/10)==0) cout << "Done " << counter_ray << " rays" << endl;

      Double_t phi = randgen.Uniform(-Pi(),Pi());
      Double_t costheta = randgen.Uniform(-1,1);
      TVector3 dir; dir.SetMagThetaPhi(1,ACos(costheta),phi);
      LightRay lr(deposit.position,dir);

      bool matched = false;
      Int_t matchx = 999;
      Int_t matchy = 999;
      Int_t firstmatchx = 999;
      Int_t firstmatchy = 999;

      Int_t index=0;
      MatchObject m;
      vector<MatchObject> matches;
      while (index>=0 && index<pars.max_distance_unit_propagation && index<Max(firstmatchx,firstmatchy)){
	if (propray(lr,kNS,index,matchx,m)){
	  if (!matched) firstmatchx=matchx;
	  matched=true;
	  matches.push_back(m);
	}
	if (propray(lr,kEW,index,matchy,m)){
	  if (!matched) firstmatchy=matchy;
	  matched=true;
	  matches.push_back(m);
	}
	index++;
      }

      vector<GoodMatchObject> goodmatches=convert_matches(matches);
	
      Double_t pathxy = pars.max_distance_unit_propagation*10*pars.xtalsize;
      GoodMatchObject finalmatch;
      for (size_t i=0; i<goodmatches.size(); i++){
	Double_t path = sqrt(pow(goodmatches[i].positionx-lr.origin.x(),2)+pow(goodmatches[i].positiony-lr.origin.y(),2));
	if (path<pathxy) {pathxy=path; finalmatch=goodmatches[i];}
      }
      if (pathxy>1e4){
	if (debug) cout << "LOOPER" << endl;
	continue;
      }

      Int_t channel = findchannel(finalmatch.x,finalmatch.y)-1;
     
      if (debug) cout << finalmatch.x << " " << finalmatch.y << " " << pathxy << " " << channel+1 << endl;

      Double_t path3d = pathxy/Sin(lr.dirvect.Theta());
     
      TVector3 rotated_origin_xy = lr.origin; rotated_origin_xy.SetZ(0);
      TVector3 rotated_impact = TVector3(finalmatch.positionx,finalmatch.positiony,0);
      TRotation r;
      r.RotateZ(-Pi()/4);
      rotated_origin_xy = r*rotated_origin_xy;
      rotated_impact = r*rotated_impact;

      Int_t nx = 0;
      Int_t ny = 0;
      Int_t nz = 0;
      {
	Int_t sx = finalmatch.x+finalmatch.y;
	if (sx==0) nx=0;
	else if (sx>0) nx = sx/2-1;
	else nx = TMath::Abs(sx)/2;
	Int_t dy = finalmatch.y-finalmatch.x;
	if (dy==0) ny=0;
	else if (dy>0) ny = dy/2-1;
	else ny = TMath::Abs(dy)/2;
	Float_t finalz = path3d*Cos(lr.dirvect.Theta())+deposit.position.z()-pars.xtalheight/2;
	nz = Int_t((fabs(finalz)+pars.xtalheight/2)/pars.xtalheight);
      }
      Int_t all_crossings = nx+ny+nz;

      Int_t my_paper_refl = 0;

      TVector2 difference = TVector2(fabs(rotated_origin_xy.x()-rotated_impact.x()),fabs(rotated_origin_xy.y()-rotated_impact.y()));
      if (difference.Phi()<limit_angle) my_paper_refl+=nx;
      if (Pi()/2-difference.Phi()<limit_angle) my_paper_refl+=ny;
      if (lr.dirvect.Theta()<limit_angle) my_paper_refl+=nz;

      Double_t att_absorption_point = randgen.Exp(pars.light_att_length);
      if (path3d>att_absorption_point) continue;

      Double_t prob_loss_in_reflections = pow(pars.eff_reflection_paper,my_paper_refl)*pow(pars.eff_reflection_total,all_crossings-my_paper_refl);
      if (randgen.Uniform()>prob_loss_in_reflections*pars.eff_lightcoll) continue;
	
      Double_t scintillation_delay = randgen.Exp(pars.scintillation_typ_timescale);

      Double_t arrival_time = deposit.time+path3d/pars.speed_of_light+scintillation_delay;
      if (arrival_time>pars.limit_arrival_time) continue;

      myphotons.at(channel)++;
      output->chamfer_pulseshape[channel]->Fill(arrival_time);
      output->reflections_paper->Fill(my_paper_refl);
      output->reflections_total->Fill(all_crossings-my_paper_refl);
      output->reflections_all->Fill(all_crossings);
      output->optical_path->Fill(path3d);
	
    }

    for (int i=0; i<4; i++) output->chamfer_photons[i]->Fill(myphotons[i]);

  }

  for (int i=0; i<4; i++) Normalize(output->chamfer_pulseshape[i]);
  Normalize(output->reflections_paper);
  Normalize(output->reflections_total);
  Normalize(output->reflections_all);
  Normalize(output->optical_path);

  return output;

};

void LightSimulator::Normalize(TH1F *h){
  h->Scale(1./h->Integral());
}

vector<GoodMatchObject>  LightSimulator::convert_matches(const vector<MatchObject> &a){
  vector<GoodMatchObject> out;
  GoodMatchObject ma;
  for (size_t i=0; i<a.size(); i++){
    if (a[i].axdir==kNS){
      if (CalculateDirection::isNbound(a[i].incomingdir)) ma.Set(2*(a[i].othercoordinate),2*(a[i].index+1));
      else ma.Set(2*(a[i].othercoordinate),2*(-a[i].index));
    }
    else{
      if (CalculateDirection::isEbound(a[i].incomingdir)) ma.Set(-1+2*(a[i].index+1),1+2*(a[i].othercoordinate));
      else ma.Set(-1+2*(-a[i].index),1+2*(a[i].othercoordinate));
    }
    assert ((ma.x+ma.y)%2==0);

    ma.positionx=ma.x*pars.xtalsize/sqrt(2); 
    ma.positiony=-pars.xtalsize/sqrt(2)+ma.y*pars.xtalsize/sqrt(2);

    if (a[i].axdir==kNS) ma.positionx+=a[i].position*(pars.xtalsize*sqrt(2));
    else ma.positiony-=a[i].position*(pars.xtalsize*sqrt(2));

    out.push_back(ma);
  }
  return out;
};
  


#endif
