//Plotting Class includes
#ifndef TRUTHPLOTS_CXX
#define TRUTHPLOTS_CXX

#include "TruthPlots.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>


TruthPlots::TruthPlots()
{
}

TruthPlots::~TruthPlots(){}

void TruthPlots::initialize(Config config)
{ 
  
  if(config.DataYears=="2015+2016"){ dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017"){ dirYear="2017/"; year="17";}
  if(config.DataYears=="2018"){ dirYear="2018/"; year="18";}
  if(config.DataYears=="Full"){ dirYear="Full/"; year="Full";}

  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}

  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
  systematic="_nominal";
 
  //zmumu
  if(config.WZSelection=="zmumu"){
    wzchannel="z"; nameOfSample="zmumu";
    fmc    = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive + nameOfSample + "_" + wzchannel + total+".root").c_str());
    fcw    = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
  }  
  //wplus
  if(config.WZSelection=="wplus" || config.WZSelection=="combined"){
    wzchannel="wplus"; nameOfSample="wplusmunu";
    fmc   = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive + nameOfSample + "_" + wzchannel + total+".root").c_str());
    fcw   = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
  }
  //wminus
  if(config.WZSelection=="wminus" || config.WZSelection=="combined"){
    wzchannel="wminus"; nameOfSample="wminmunu";
    fmc    = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive + nameOfSample + "_" + wzchannel + total+".root").c_str());
    fcw    = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
    fmcp   = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive + nameOfSample + "_" + wzchannel + total+".root").c_str());
    fcwp   = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
  }                

  if(!fmc->IsOpen() || !fcw->IsOpen() || !fmcp->IsOpen() || !fcwp->IsOpen()){
    std::cout<<"Truth plots on, but not enough files to produce them, bye!"<<std::endl;
    exit(1);
  }

  return;
}

void TruthPlots::execute(Config config){

  std::vector<std::string> kine_simple;
  kine_simple.push_back("mw");
  kine_simple.push_back("mwtw");
  kine_simple.push_back("WZ_m");
  //kine_simple.push_back("WZ_mt");
  for(int i=0; i<kine_simple.size();i++){
    //SinglePlot(config, kine_simple[i],wzchannel);
    MatchPlot(config, kine_simple[i] ,wzchannel);
  }

  std::vector<std::string> kine;
  kine.push_back("m"  );
  kine.push_back("pt" ); 
  kine.push_back("eta"); 
  kine.push_back("phi");  
  for(int j=0; j<kine.size();j++)
    LevelsPlot(config,kine[j],wzchannel);

  std::vector<std::string> cw_factor;
  //cw_factor.push_back("hCwz");
  cw_factor.push_back("hCwzEta");

  //cw_factor.push_back("hPtywz");
  cw_factor.push_back("hPtywzEta");

  //  cw_factor.push_back("hStabwz");
  cw_factor.push_back("hStabwzEta");
  cw_factor.push_back("hEta");

  for(int k=0; k<cw_factor.size(); k++)
    CwPlots(config, cw_factor[k]);

  return;
}


void TruthPlots::SinglePlot(Config config, std::string kine, std::string boson){

  TH1D *hPlotTemp, *hTemp1, *hPlot;
  double GeV=1./1000.;
  if(boson=="zmumu")  hPlotTemp = (TH1D*)fmc ->GetObjectUnchecked(("h_"+kine).c_str());
  if(boson=="wplus" || config.WZSelection=="combined")  hPlotTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_"+kine).c_str());
  if(boson=="wminus" && config.WZSelection!="combined") hPlotTemp  = (TH1D*)fmcp ->GetObjectUnchecked(("h_"+kine).c_str());
  if(boson=="wminus" && config.WZSelection=="combined"){
    hTemp1  = (TH1D*)fmcp ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotTemp->Add(hTemp1);
  }

  hPlot = scaleXaxis(hPlotTemp, GeV, kine);
  if(hPlot==NULL) exit(10);

  std::string xlabel,ylabel; double xlow=0., xhigh=0.;
  TString leglabel;

  if(boson=="zmumu") {xlabel="m_{Z} [GeV]"; xlow=60000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="wminus"|| boson=="wplus"){xlabel="m_{W} [GeV]"; xlow=40000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="zmumu")  leglabel="MC Z#rightarrow#mu#mu";
  if(boson=="wplus")  leglabel="MC W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus") leglabel="MC W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="combined") leglabel="MC W#rightarrow#mu#nu";
  ylabel="Events";

  plotAxisLine(hPlot,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0.);

  TString lumi="";
  TCanvas *c = new TCanvas("canvas","canvas",550,450);
  TLegend* leg = new TLegend(0.6,0.85,0.8,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry(hPlot,  "#bf{#scale[0.85]{"+leglabel+", Born boson}}","L");  

  TPaveText *box;
  box = new TPaveText(0.625,0.865,0.8,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi+" Simulation");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  pad->SetLogy();
  pad->SetLogx();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.1125);

  hPlot->SetTitle("");
  hPlot->GetXaxis()->SetTitleOffset(1.45);
  hPlot->DrawCopy("hist");  

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  if(config.WZSelection=="combined") boson = "combined";
  c->Print((config.OutputFileDir + "Plots/"+ dirYear + "Truth/"+"simple_truth_" +boson+"_"+kine+".pdf").c_str());  
  delete hPlot;

  return;
}


void TruthPlots::MatchPlot(Config config, std::string kine, std::string boson){

  TFile *ftemp, *ftemp2;
  std::vector<TH1D*> mass_slices; mass_slices.clear();
  TH1D *hTemp, *hTemp2, *hTempFinal, *hTempFinal2;
  double GeV=1./1000.;
  std::vector<std::string> slices; slices.clear();
  std::vector<int> color; color.clear();

  slices.push_back("");      color.push_back(1);
  slices.push_back("_120");  color.push_back(2);   
  slices.push_back("_180");  color.push_back(3);  
  slices.push_back("_250");  color.push_back(4);  
  slices.push_back("_400");  color.push_back(5);      
  slices.push_back("_600");  color.push_back(7);      
  slices.push_back("_800");  color.push_back(8);  
  slices.push_back("_1000"); color.push_back(9); 
  slices.push_back("_1250"); color.push_back(1);
  slices.push_back("_1500"); color.push_back(2);
  slices.push_back("_1750"); color.push_back(3);
  slices.push_back("_2000"); color.push_back(4); 
  slices.push_back("_2250"); color.push_back(5);
  slices.push_back("_2500"); color.push_back(6);
  slices.push_back("_2750"); color.push_back(7);
  slices.push_back("_3000"); color.push_back(8);
  slices.push_back("_3500"); color.push_back(9); 
  slices.push_back("_4000"); color.push_back(1); 
  slices.push_back("_4500"); color.push_back(2); 
  slices.push_back("_5000"); color.push_back(3);

  for(int k=0;k<(int)slices.size();k++){
    ftemp  = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/"+ nameOfSample+slices[k] + "_" + wzchannel+".root").c_str());
    hTemp = (TH1D*)ftemp ->GetObjectUnchecked(("h_"+kine).c_str());
    hTempFinal = scaleXaxis(hTemp, GeV, kine);
    if(config.WZSelection=="combined"){
      ftemp2 = new TFile((config.OutputFileDir +"Files/"+ dirYear +"Truth/wplusmunu"+slices[k] + "_wplus.root").c_str());
      hTemp2 = (TH1D*)ftemp2 ->GetObjectUnchecked(("h_"+kine).c_str());
      hTempFinal2 = scaleXaxis(hTemp2, GeV, kine);
      hTempFinal->Add(hTempFinal2);
    }      
    mass_slices.push_back(hTempFinal);  
  }

  std::string xlabel,ylabel; double xlow=0., xhigh=0.;
  TString leglabel;

  if(boson=="zmumu") {xlabel="m_{Z} [GeV]"; xlow=60000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="wminus"|| boson=="wplus"){xlabel="m_{W} [GeV]"; xlow=40000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="zmumu")  leglabel="MC Z#rightarrow#mu#mu";
  if(boson=="wplus")  leglabel="MC W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus") leglabel="MC W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="combined") leglabel="MC W#rightarrow#mu#nu";
  ylabel="Events";

  if(kine=="m") xlabel="m_{W} [GeV]";
  if(kine=="mwtw") xlabel="m_{T,W} [GeV]";

  for(int k=1;k<(int)mass_slices.size();k++)
    plotAxisLine(mass_slices[k],color[k],color[k],23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0.);
  plotAxisLine(mass_slices[0],color[0],color[0],23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,10e-8,10e8);

  TString lumi="";
  TCanvas *c = new TCanvas("canvas","canvas",550,450);
  TLegend* leg = new TLegend(0.6,0.85,0.8,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry((TH1D*)mass_slices[0],  "#bf{#scale[0.85]{"+leglabel+", Born boson}}","L");  

  TPaveText *box;
  box = new TPaveText(0.625,0.865,0.8,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi+" Simulation");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  pad->SetLogy();
  pad->SetLogx();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.1125);

  for(int k=0;k<(int)mass_slices.size();k++){
    mass_slices[k]->GetXaxis()->SetTitleOffset(1.45);
    if(k==0) 
      mass_slices[k]->DrawCopy("hist");
    else
      mass_slices[k]->DrawCopy("sames hist");
  }
  //hPlot->DrawCopy("p sames hist");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  if(config.WZSelection=="combined") boson = "combined";
  c->Print((config.OutputFileDir + "Plots/"+ dirYear + "Truth/"+"slices_truth_" +boson+"_"+kine+".pdf").c_str());  

  return;
}


void TruthPlots::LevelsPlot(Config config, std::string kine, std::string boson)
{
  TH1D *hBoson, *hDilBorn, *hDilBare, *hDilDres;
  TH1D *hBosonP, *hDilBornP, *hDilBareP, *hDilDresP;
  TH1D *hBosonTemp, *hDilBornTemp, *hDilBareTemp, *hDilDresTemp;
  TH1D *hBosonTempP, *hDilBornTempP, *hDilBareTempP, *hDilDresTempP;
  double GeV=1./1000.;

  if(boson=="zmumu"){  
    hBosonTemp    = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_"+kine).c_str());
    hDilBornTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_born").c_str());
    hDilBareTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_bare").c_str());
    hDilDresTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_dres").c_str());
  }
  if(boson=="wplus" || config.WZSelection=="combined"){
    hBosonTemp    = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_"+kine).c_str());
    hDilBornTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_born").c_str());
    hDilBareTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_bare").c_str());
    hDilDresTemp  = (TH1D*)fmc  ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_dres").c_str());
  }
  if(boson=="wminus" && config.WZSelection!="combined"){
    hBosonTemp    = (TH1D*)fmc ->GetObjectUnchecked(("h_WZ_"+kine).c_str());
    hDilBornTemp  = (TH1D*)fmc ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_born").c_str());
    hDilBareTemp  = (TH1D*)fmc ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_bare").c_str());
    hDilDresTemp  = (TH1D*)fmc ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_dres").c_str());
  }
  if(boson=="wminus" && config.WZSelection=="combined"){
    hBosonTempP    = (TH1D*)fmcp ->GetObjectUnchecked(("h_WZ_"+kine).c_str());
    hDilBornTempP  = (TH1D*)fmcp ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_born").c_str());
    hDilBareTempP  = (TH1D*)fmcp ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_bare").c_str());
    hDilDresTempP  = (TH1D*)fmcp ->GetObjectUnchecked(("h_WZ_dilep_"+kine+"_dres").c_str());
    hBosonTemp->Add(hBosonTempP); hDilBornTemp->Add(hDilBornTempP); hDilBareTemp->Add(hDilBareTempP); hDilDresTemp->Add(hDilDresTempP);
  }

  if(hBosonTemp==NULL) exit(1);
  hBoson    = scaleXaxis(hBosonTemp,   GeV, kine);
  hDilBorn  = scaleXaxis(hDilBornTemp, GeV, kine);
  hDilBare  = scaleXaxis(hDilBareTemp, GeV, kine);
  hDilDres  = scaleXaxis(hDilDresTemp, GeV, kine);

  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0.;
  TString leglabel;

  if(kine=="m" && boson=="zmumu") {xlabel="m_{Z} [GeV]"; xlow=60000.*GeV; xhigh=7000000.*GeV;}
  if(kine=="m" && (boson=="wminus"||boson=="wplus")){
    xlabel="m_{W} [GeV]";
    xlow = hBoson->GetXaxis()->GetBinLowEdge(1); 
    xhigh= hBoson->GetXaxis()->GetBinUpEdge(hBoson->GetNbinsX()); 
  }
  if(kine=="pt"){xlabel="p_{T} [GeV]"; xlow=0.; xhigh=200*1000.*GeV;}
  if(kine=="eta"){xlabel="#eta"; xlow=-6.; xhigh = 6.;}
  if(kine=="phi"){xlabel="#phi [rad]"; xlow=-5.; xhigh=5.;}
  if(boson=="zmumu")  leglabel="MC Z#rightarrow#mu#mu";
  if(boson=="wplus")  leglabel="MC W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus") leglabel="MC W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="combined") leglabel="MC W#rightarrow#mu#nu";

  ratiolabel="boson/dilepton"; 

  ylabel="Events";
    
  plotAxisLine(hBoson,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,4.);
  plotAxisLine(hDilBorn,kMagenta+2,kMagenta+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(hDilBare,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(hDilDres,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="44.3";
  if(config.DataYears=="2018")      lumi="59.9";
  if(config.DataYears=="2015+2016") lumi="36.1";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.55,0.85,0.75,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry(hBoson,  "#bf{#scale[0.85]{"+leglabel+", Born boson}}","L");  
  leg->AddEntry(hDilBorn,"#bf{#scale[0.85]{"+leglabel+", Dilepton Born}}","PL");  
  leg->AddEntry(hDilBare,"#bf{#scale[0.85]{"+leglabel+", Dilepton bare}}","PL");
  leg->AddEntry(hDilDres,"#bf{#scale[0.85]{"+leglabel+", Dilepton dress}}","PL");  

  TPaveText *box;
  box = new TPaveText(0.575,0.865,0.75,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi+" fb^{-1}");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  if(kine=="m"){
    pad->SetLogy();
    pad->SetLogx();}
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);

  hDilBare->SetTitle("");
  hDilBare->GetXaxis()->SetTitleOffset(1.45);
  hBoson->DrawCopy("hist");  
  hDilBare->DrawCopy("p sames hist");
  hDilBorn->DrawCopy("hist p sames");
  hDilDres->DrawCopy("p sames hist");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  if(kine=="m")padr->SetLogx();
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = hBoson->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[hBoson->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hBoson->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hBoson->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hBoson->GetXaxis()->GetNbins()]=hBoson->GetXaxis()->GetBinUpEdge(hBoson->GetXaxis()->GetNbins());

  TH1D* ratiodataMCor   = new TH1D("ratiodataMCor","",nx,xbins);
  TH1D* ratiodataMC   = new TH1D("ratiodataMC","",nx,xbins);
  TH1D* ratiodataMCmu = new TH1D("ratiodataMCmu","",nx,xbins);
  TH1D* ratioMCmu = new TH1D("ratioMCmu","",nx,xbins);

  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){

    double data=hBoson->GetBinContent(b);
    //if(data==0) continue;
    double mcor=hDilBorn->GetBinContent(b);
    if(mcor==0) continue;
    double mc=hDilBare->GetBinContent(b);
    if(mc==0) continue;
    double mcmu=hDilDres->GetBinContent(b);
    if(mcmu==0) continue;

    double dataratioor=data/mcor;
    double dataratio=data/mc;
    double dataratiocalib=data/mcmu;

    double err=hBoson->GetBinError(b)/mc;
    double error=hBoson->GetBinError(b)/mc;
    double errmu=hBoson->GetBinError(b)/mcmu;
    
    ratiodataMCor->SetBinContent(b,dataratioor);
    ratiodataMCor->SetBinError(b,err);

    ratiodataMC->SetBinContent(b,dataratio);
    ratiodataMC->SetBinError(b,err);

    ratiodataMCmu->SetBinContent(b,dataratiocalib);
    ratiodataMCmu->SetBinError(b,errmu);
  }

  plotAxisLine(ratiodataMCor,kCyan,kCyan,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratioMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  
  //ratioSettings(ratiodataMCor  ,0.79,1.21,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kMagenta+2,0.5,10); 
  ratioSettings(ratiodataMCor  ,0.49,1.51,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kMagenta+2,0.5,10); 
  ratioSettings(ratiodataMC    ,0.49,1.51,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kCyan+3,0.5,10);
  ratioSettings(ratiodataMCmu  ,0.49,1.51,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 
  ratioSettings(ratioMCmu      ,0.49,1.51,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 

  ratiodataMC->DrawCopy("p");
  ratiodataMCor->DrawCopy("p same");
  ratiodataMCmu->DrawCopy("p same");
  //ratioMCmu->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,1.,xhigh,1.);
  //line = new TLine(40000.,1.,150000.,1.);//remove
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  if(config.WZSelection=="combined") boson="combined";
  c->Print((config.OutputFileDir + "Plots/"+dirYear+ "Truth/truth_"+boson+"_"+kine+".pdf").c_str());
  
  delete hBoson; delete hDilBorn; delete hDilBare; delete hDilDres;
  return;
}


void TruthPlots::CwPlots(Config config, std::string factor){

  TH1D *hCw, *hCwEta, *hPurity, *hPurityEta, *hStab, *hStabEta, *hEta;  
  //hCw        = (TH1D*)fcw ->GetObjectUnchecked("hCwz");	  
  hCwEta     = (TH1D*)fcw ->GetObjectUnchecked("hCwzEta");	  
  //hPurity    = (TH1D*)fcw ->GetObjectUnchecked("hPtywz");	  
  hPurityEta = (TH1D*)fcw ->GetObjectUnchecked("hPtywzEta"); 
  //hStab      = (TH1D*)fcw ->GetObjectUnchecked("hStabwz");	  
  hStabEta   = (TH1D*)fcw ->GetObjectUnchecked("hStabwzEta");
  hEta       = (TH1D*)fcw ->GetObjectUnchecked("hEtaCgen");

  TH1D *hPlot;
  if(factor=="hCwzEta") hPlot=(TH1D*)hCwEta->Clone("hCwzEta_cloned");
  if(factor=="hPtywzEta") hPlot=(TH1D*)hPurityEta->Clone("hPurityEta_cloned");
  if(factor=="hStabwzEta") hPlot=(TH1D*)hStabEta->Clone("hStabEta_cloned");
  if(factor=="hEta") hPlot=(TH1D*)hEta->Clone("hEta_cloned");

  if(factor=="hEta")
    for(int k=1;k<hPlot->GetXaxis()->GetNbins()+1;k++)
      hPlot->SetBinContent(k,(hPlot   ->GetBinContent(k)/hPlot   ->GetBinWidth(k)));

  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0.;
  TString leglabel;
  if(config.WZSelection=="wminus") leglabel = "W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="wplus")  leglabel = "W^{+}#rightarrow#mu^{+}#nu";


  if(factor=="hCwz" )     {xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=200000.; ylabel="C_{w}";}
  if(factor=="hCwzEta")   {xlabel="#||{#eta^{#mu}}"; xlow=-2.5; xhigh=2.5; ylabel="C_{w}";}
  if(factor=="hPtywz")    {xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=200000.; ylabel="Purity";}
  if(factor=="hPtywzEta") {xlabel="#||{#eta^{#mu}}"; xlow=-2.5; xhigh=2.5;  ylabel="Purity";}
  if(factor=="hStabwz")   {xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=200000.;  ylabel="Stability";}
  if(factor=="hStabwzEta"){xlabel="#||{#eta^{#mu}}"; xlow=-2.5; xhigh=2.5;  ylabel="Stability";}
  if(factor=="hEta"){xlabel="#||{#eta^{#mu}}"; xlow=-2.5; xhigh=2.5;  ylabel="Events";}

  plotAxisLine(hPlot,kCyan+3,kCyan+3,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,0,2.4,false,0.,1.1);

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2018")      lumi="59.9";
  if(config.DataYears=="2015+2016") lumi="36.1";

  TCanvas *c = new TCanvas("canvas","canvas",550,550);
  TLegend* leg = new TLegend(0.6,0.775,0.8,0.725,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  //  leg->AddEntry(hPlot, "#bf{#scale[0.85]{"+leglabel+", Born level}}","L");  
  leg->AddEntry(hPlot, "#bf{#scale[0.85]{"+leglabel+", Born level}}","L");

  TPaveText *box;
  box = new TPaveText(0.625,0.815,0.8,0.915,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, Simulation");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.05);

  hPlot->SetTitle("");
  hPlot->GetXaxis()->SetTitleOffset(1.45);

  if(factor=="hPtywzEta")
    hPlot->GetYaxis()->SetRangeUser(0,1.5);
  else
    hPlot->GetYaxis()->SetRangeUser(0,1.25);

  if(factor=="hEta")
    hPlot->GetYaxis()->SetRangeUser(0,1.2e8);

  hPlot->SetLineWidth(3);
  hPlot->DrawCopy("hist");  

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  c->cd();
  c->Print((config.OutputFileDir +"Plots/"+dirYear+"Truth/unfolding_"+factor+"_"+config.WZSelection+".pdf").c_str());
  
  delete hPlot;
  return;
}



TH1D* TruthPlots::scaleXaxis(TH1D* hist, double scale, std::string kine)
{
  if(scale == 1) return hist;
  double binMin = hist->GetXaxis()->GetXmin();
  double binMax = hist->GetXaxis()->GetXmax();
  double nBins  = hist->GetXaxis()->GetNbins();
  double binMinScaled = binMin*scale;
  double binMaxScaled = binMax*scale;
  TH1D* scaled = NULL;

  if(kine=="m" || kine=="mw" || kine=="mwtw" || kine=="WZ_m"){ 
    binMinScaled = log10(binMinScaled); binMaxScaled = log10(binMaxScaled);
    double scaledBins[(int)nBins+1];
    double binWidth = (binMaxScaled - binMinScaled)/(nBins);
    for (int i = 0; i < nBins+1; i++){
      scaledBins[i  ] = binMinScaled+(i  )*binWidth;
      scaledBins[i  ] = pow(10,scaledBins[i  ]);
    }
    scaled = new TH1D(hist->GetName(),hist->GetTitle(),nBins,scaledBins);
  }else {
    scaled = new TH1D(hist->GetName(),hist->GetTitle(),nBins,binMinScaled,binMaxScaled);
  }
  
  for(unsigned int b = 1; b< nBins+1; b++) // look at bin numbering scheme documented in TH1::GetBin()
    {
      double val = hist->GetBinContent(b);
      double err = hist->GetBinError(b);
      scaled->SetBinContent(b,val);
      scaled->SetBinError(b,err);
    }
  return scaled;
}


void TruthPlots::plotAxisLine(TH1D* hist, int lineColor, int markerColor,
			      int markerStyle, double markerSize,
			      TString title, TString xlabel, TString ylabel, bool xRange,
			      double xlow, double xhigh, bool yRange, double ylow, double yhigh)
{
  hist->SetLineColor(lineColor);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerColor(markerColor);
  hist->SetMarkerSize(markerSize);
  hist->SetTitle(title);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);

  if(xRange==true){
    hist->GetXaxis()->SetRangeUser(xlow,xhigh);   /*std::cout<<xlow<<"    Parrito"<<xhigh<<std::endl;*/}
  if(yRange==true)
    hist->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}


void TruthPlots::profAxisLine(TProfile* prof, int lineColor, int markerColor,
			      int markerStyle, double markerSize,
			      TString title, TString xlabel, TString ylabel, bool xRange,
			      double xlow, double xhigh, bool yRange, double ylow, double yhigh)
{
  prof->SetLineColor(lineColor);
  prof->SetMarkerStyle(markerStyle);
  prof->SetMarkerColor(markerColor);
  prof->SetMarkerSize(markerSize);
  prof->SetTitle(title);
  prof->GetXaxis()->SetTitle(xlabel);
  prof->GetYaxis()->SetTitle(ylabel);
  if(xRange==true)
    prof->GetXaxis()->SetRangeUser(xlow,xhigh);
  if(yRange==true)
    prof->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}


void TruthPlots::ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
			       double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
			       int markerstyle, int color, double markersize, int divisions)
{
  hist->SetMinimum(min);
  hist->SetMaximum(max);
  hist->GetXaxis()->SetTitle(xlabel);

  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->SetLabelSize(xlabelsize);
  hist->GetYaxis()->SetLabelSize(ylabelsize);
  hist->GetYaxis()->SetTitleOffset(ytitleof);
  hist->GetYaxis()->SetTitleSize(ytitlesize);
  hist->SetLineColor(color);

  hist->Sumw2();
  hist->SetMarkerStyle(markerstyle);
  hist->SetMarkerSize(markersize);
  hist->SetMarkerColor(color);

  hist->SetStats(0);
  hist->GetYaxis()->SetNdivisions(divisions);

  return;
}


void TruthPlots::profRatioSettings(TProfile* prof, double min, double max, TString xlabel, TString ylabel,
			       double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
			       int markerstyle, int color, double markersize, int divisions)
{
  prof->SetMinimum(min);
  prof->SetMaximum(max);
  prof->GetXaxis()->SetTitle(xlabel);
  prof->GetYaxis()->SetTitle(ylabel);
  prof->GetXaxis()->SetLabelSize(xlabelsize);
  prof->GetYaxis()->SetLabelSize(ylabelsize);
  prof->GetYaxis()->SetTitleOffset(ytitleof);
  prof->GetYaxis()->SetTitleSize(ytitlesize);
  prof->SetLineColor(color);

  prof->Sumw2();
  prof->SetMarkerStyle(markerstyle);
  prof->SetMarkerSize(markersize);
  prof->SetMarkerColor(color);

  prof->SetStats(0);
  prof->GetYaxis()->SetNdivisions(divisions);

  return;
}

void TruthPlots::setstyle(){

  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(000000);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelFont(22,"x");
  gStyle->SetLabelFont(22,"y");

  gStyle->SetTitleSize(0.035,"x");
  gStyle->SetTitleSize(0.035,"y");
  gStyle->SetTitleOffset(1.35,"x");
  gStyle->SetTitleOffset(1.35,"y");

  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.18);

  int font = 42;
  gStyle->SetTextFont(font);

  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");

  return;
}

#endif
