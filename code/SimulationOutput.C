#ifndef __SIMULATIONOUTPUT_C__
#define __SIMULATIONOUTPUT_C__

#include "SimulationOutput.h"

using namespace TMath;
using namespace std;

SimulationOutput::SimulationOutput(EnergyDeposit deposit_):deposit(deposit_){
    for (int i=0; i<4; i++) {
      chamfer_photons[i] = new TH1F(Form("photons_ch%d",i+1),Form("photons_ch%d",i+1),int(deposit.totalrays)/10,0,deposit.totalrays);
      chamfer_photons[i]->SetTitle(Form("Photons chamfer ch%d",i+1));
      chamfer_photons[i]->GetXaxis()->SetTitle("Photons");
      chamfer_pulseshape[i] = new TH1F(Form("pulse_ch%d",i+1),Form("pulse_ch%d",i+1),1000,0,pars.limit_arrival_time);
      chamfer_pulseshape[i]->SetTitle(Form("Pulseshape ch%d",i+1));
      chamfer_pulseshape[i]->GetXaxis()->SetTitle("t (ps)");
    }
    reflections_all = new TH1F("reflections_all","reflections_all",60,0,60);
    reflections_total = new TH1F("reflections_total","reflections_total",60,0,60);
    reflections_total->SetLineColor(kGreen);
    reflections_paper = new TH1F("reflections_paper","reflections_paper",60,0,60);
    reflections_paper->SetLineColor(kRed);
    optical_path = new TH1F("optical_path","optical_path",200*10,0,200);
    optical_path->GetXaxis()->SetTitle("mm");
    reflections_all->Sumw2();
    reflections_total->Sumw2();
    reflections_paper->Sumw2();
    optical_path->Sumw2();
  };

void SimulationOutput::Write(){

  TCanvas *c1 = new TCanvas("c1","Photons",1200,600);
  c1->Divide(4,2);
  for (int i=0; i<4; i++) {c1->cd(i+1); chamfer_photons[i]->Draw();}
  for (int i=0; i<4; i++) {c1->cd(4+i+1); chamfer_pulseshape[i]->Draw();}
  c1->Write();

  TCanvas *c2 = new TCanvas("c2","Reflections",800,400);
  c2->Divide(2);
  c2->cd(1);
  reflections_all->Draw("E1");
  reflections_all->GetYaxis()->SetRangeUser(0,0.5);
  reflections_total->Draw("sameE1");
  reflections_paper->Draw("sameE1");
  c2->cd(2);
  optical_path->Draw();
  c2->Write();

  TCanvas *c3 = new TCanvas("c3","Sharing",800,400);
  c3->Divide(2);
  c3->cd(1);

  TH2F *h2d = new TH2F("h2d","Photon distribution in channels",2,-pars.xtalsize/2,pars.xtalsize/2,2,-pars.xtalsize/2,pars.xtalsize/2);
  h2d->SetBinContent(1,1,chamfer_photons[3]->GetMean());
  h2d->SetBinContent(2,1,chamfer_photons[2]->GetMean());
  h2d->SetBinContent(2,2,chamfer_photons[1]->GetMean());
  h2d->SetBinContent(1,2,chamfer_photons[0]->GetMean());
  h2d->SetBinError(1,1,sqrt(chamfer_photons[3]->GetMean()));
  h2d->SetBinError(2,1,sqrt(chamfer_photons[2]->GetMean()));
  h2d->SetBinError(2,2,sqrt(chamfer_photons[1]->GetMean()));
  h2d->SetBinError(1,2,sqrt(chamfer_photons[0]->GetMean()));
  h2d->SetStats(kFALSE);
  h2d->Draw("text E1");


  TMarker *mark = new TMarker(deposit.getpositionlabframe().x(),deposit.getpositionlabframe().y(),20);
  mark->SetMarkerColor(kRed);

  TPaveText *total1 = new TPaveText(-pars.xtalsize/8,-pars.xtalsize/2,pars.xtalsize/8,-pars.xtalsize/2+pars.xtalsize/6,"");
  total1->AddText("Total");
  total1->AddText(Form("%.1f #gamma",h2d->Integral()));

  total1->Draw();
  mark->Draw();
  

  c3->cd(2);
  TH2F *h2r = (TH2F*)(h2d->Clone("h2r"));
  h2r->SetTitle("Photon sharing among channels");
  h2r->Scale(1./h2r->Integral());
  h2r->GetZaxis()->SetRangeUser(0,1);
  h2r->SetContour(50);
  h2r->SetStats(kFALSE);
  h2r->Draw("text colz E1");

  TPaveText *total2 = new TPaveText(-pars.xtalsize/8,-pars.xtalsize/2,pars.xtalsize/8,-pars.xtalsize/2+pars.xtalsize/6,"");
  total2->AddText("Total");
  total2->AddText(Form("%.3f",h2d->Integral()/deposit.totalrays));

  total2->Draw();
  mark->Draw();

  c3->Write();

  deposit.position.Write("deposit_position");
  TParameter<float> *myE = new TParameter<float>("deposit_energy",deposit.energy);
  TParameter<float> *myT = new TParameter<float>("deposit_time",deposit.time);
  TParameter<int> *myN = new TParameter<int>("deposit_totalrays",deposit.totalrays);
  myE->Write();
  myT->Write();
  myN->Write();

  };



#endif
