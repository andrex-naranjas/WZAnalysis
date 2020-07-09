//Plotting Class includes
#ifndef SFPLOTS_CXX
#define SFPLOTS_CXX

#include "SFPlots.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>


SFPlots::SFPlots()
{
}

SFPlots::~SFPlots(){}

void SFPlots::initialize(Config config)
{ 
  
  if(config.DataYears=="2015+2016"){ dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017"){ dirYear="2017/"; year="17";}
  if(config.DataYears=="2018"){ dirYear="2018/"; year="18";}
  if(config.DataYears=="Full"){ dirYear="Full/"; year="Full";}

  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}
  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
 
  //zmumu
  if(config.WZSelection=="zmumu"){ wzchannel="z"; nameOfSample="zmumu";}  
  //wplus
  if(config.WZSelection=="wplus" || config.WZSelection=="combined"){ wzchannel="wplus"; nameOfSample="wplusmunu";}
  //wminus
  if(config.WZSelection=="wminus" || config.WZSelection=="combined"){ wzchannel="wminus"; nameOfSample="wminmunu";}


  fnom  = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + nameOfSample + "_" + wzchannel + "_nominal" + puname +total+".root").c_str());  

  return;
}

void SFPlots::execute(Config config){
  
  std::vector<std::string> kine_simple;
  kine_simple.push_back("mwt");
  kine_simple.push_back("pt");
  kine_simple.push_back("eta");
  kine_simple.push_back("m");

  std::vector<std::string> SF;
  SF.push_back("Iso");
  SF.push_back("ID");
  SF.push_back("Tri");
  SF.push_back("ttva");
  SF.push_back("PU");
  SF.push_back("jvt_PU");
  SF.push_back("KF");


  for(int i=0; i<kine_simple.size();i++)
    for(int k=0; k<SF.size();k++)
      SinglePlot(config, kine_simple[i], wzchannel, SF[k]);
  
  //for(int i=0; i<kine_simple.size();i++) CombinedPlot(config, kine_simple[i], wzchannel);

  return;
}


void SFPlots::SinglePlot(Config config, std::string kine, std::string boson, std::string SF){

  bool isoFlag=false, idFlag=false, trigFlag=false, ttvaFlag=false, puFlag=false, kFlag=false, jvtFlag=false;

  if(SF=="Iso")    isoFlag=true;
  if(SF=="ID" )    idFlag =true;
  if(SF=="Tri")    trigFlag=true;
  if(SF=="ttva")   ttvaFlag=true;
  if(SF=="PU" )    puFlag=true;
  if(SF=="jvt_PU") jvtFlag=true;
  if(SF=="KF" )    kFlag=true;

  if(!isoFlag && !idFlag && !trigFlag && !ttvaFlag && !puFlag && !jvtFlag && !kFlag) return;

  TString name = "";
  if(isoFlag)  name = "Isolation scale factor";
  if(idFlag)   name = "Identification scale factor";
  if(trigFlag) name = "Trigger scale factor";
  if(ttvaFlag) name = "TTVA scale factor";
  if(puFlag)   name = "Pileup weight";
  if(jvtFlag)  name = "JVT weigth";
  if(kFlag)    name = "kFactor weight";
  
  TH1D *hNom, *hSysUp, *hSysDown, *hStatUp, *hStatDown;
  TProfile *hNomAux, *hSysUpAux, *hSysDownAux;
  TProfile *hStatUpAux, *hStatDownAux;

  double GeV=1./1.;
  if(kine=="mwt"|| kine=="pt" || kine=="m"){GeV=1./1000.; std::cout<<"Perrito"<<std::endl;}

  hNomAux     = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF).c_str());
  hNom     = scaleXaxis(hNomAux,     GeV, kine);

  if(isoFlag || idFlag || trigFlag || ttvaFlag){
    hSysUpAux   = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_sys_UP").c_str());
    hSysDownAux = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_sys_DOWN").c_str());
    hStatUpAux  = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_stat_UP").c_str());
    hStatDownAux= (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_stat_DOWN").c_str());
    hSysUp   = scaleXaxis(hSysUpAux,   GeV, kine);
    hSysDown = scaleXaxis(hSysDownAux, GeV, kine);
    hStatUp  = scaleXaxis(hStatUpAux,  GeV, kine);
    hStatDown= scaleXaxis(hStatDownAux,GeV, kine);
  }else if(puFlag || jvtFlag){
    hSysUpAux   = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_UP").c_str());
    hSysDownAux = (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_"+SF+"_DOWN").c_str());
    hSysUp   = scaleXaxis(hSysUpAux,   GeV, kine);
    hSysDown = scaleXaxis(hSysDownAux, GeV, kine);
    hStatUp  = NULL;
    hStatDown= NULL;
  }else if(kFlag){
    hSysUp   = NULL;
    hSysDown = NULL;
    hStatUp  = NULL;
    hStatDown= NULL;  
  }
  
  if(hNom==NULL) exit(10);

  std::string xlabel,ylabel; double xlow=0., xhigh=0., ylow=0., yhigh=0.;
  TString leglabel;

  if(boson=="zmumu") {xlabel="m_{Z} [GeV]"; xlow=60000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="wminus"|| boson=="wplus"){xlow=40000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="zmumu")  leglabel="Z#rightarrow#mu#mu";
  if(boson=="wplus")  leglabel="W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus") leglabel="W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="combined") leglabel="MC W#rightarrow#mu#nu";
  ylabel="Scale Factor";

  if(kine=="mwt")xlabel="m_{T} [GeV]";
  if(kine=="m")  xlabel="m_{inv} [GeV]";
  if(kine=="pt") xlabel="p_{T} [GeV]";
  if(kine=="eta")xlabel="#eta^{#mu}";

  if(trigFlag && kine=="eta"){ylow=0.8; yhigh=1.1;}
  if(trigFlag && (kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.1;}

  if(idFlag && kine=="eta"){ylow=0.8; yhigh=1.1;}
  if(idFlag && (kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.1;}

  if(isoFlag && kine=="eta"){ylow=0.8; yhigh=1.1;}
  if(isoFlag && (kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.1;}


  if(kFlag && boson=="wminus" && (kine=="pt" || kine=="mwt" || kine=="m")) {ylow=0.5; yhigh=1.5;}
  else{ylow=0.8; yhigh=1.2;}

  if(puFlag){ylow=0.9; yhigh=1.1; hNom->Scale(0.9); hSysUp->Scale(0.9); hSysDown->Scale(0.9);}
  if(jvtFlag){ylow=0.975; yhigh=1.05;}
  if(trigFlag){ylow=0.75; yhigh=1.1;}
  if(ttvaFlag){ylow=0.975; yhigh=1.025;}
  if(isoFlag){ylow=0.975; yhigh=1.025;}


  plotAxisLine(hNom,kBlack,kBlack,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  if(hSysUp!=NULL)
    plotAxisLine(hSysUp,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  if(hSysDown!=NULL)
    plotAxisLine(hSysDown,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  if(hStatUp!=NULL)
    plotAxisLine(hStatUp,kRed,kRed,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  if(hStatDown!=NULL)
    plotAxisLine(hStatDown,kRed,kRed,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);

  TString lumi="";
  TCanvas *c = new TCanvas("canvas","canvas",550,500);
  TLegend* leg = new TLegend(0.175,0.85,0.375,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry(hNom,      "#bf{#scale[1.0]{"+name+"}}","P");  
  if(hSysUp!=NULL)
    leg->AddEntry(hSysUp,   "#bf{#scale[1.0]{Sys.  uncertainty}}","L");  
  if(hStatUp!=NULL)
    leg->AddEntry(hStatUp,  "#bf{#scale[1.0]{Stat. uncertainty}}","L");  

  TPaveText *box;
  box = new TPaveText(0.2,0.865,0.375,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi+" Simulation");
  c->cd();

  TPaveText *box2;
  box2 = new TPaveText(0.7,0.845,0.8,0.965,"NDC");
  box2->SetBorderSize(0); box2->SetTextSize(0.05); box2->SetFillColor(0);
  box2->AddText(""+ leglabel+"");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  if(kine=="mwt" || kine=="pt" || kine=="m") pad->SetLogx();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.1125);

  hNom->SetTitle("");
  hNom->GetXaxis()->SetTitleOffset(1.45);  
  hNom->GetYaxis()->SetTitleOffset(1.65);  
  hNom->DrawCopy("pe");  
  if(hSysUp!=NULL)
    hSysUp->DrawCopy("same hist");
  if(hSysDown!=NULL)
    hSysDown->DrawCopy("same hist");  
  if(hStatUp!=NULL)
    hStatUp->DrawCopy("same hist");
  if(hStatDown!=NULL)
    hStatDown->DrawCopy("same hist");  

  pad->Update();
  pad->Modified();
  box->Draw();
  box2->Draw();
  leg->Draw();
  c->cd();

  if(config.WZSelection=="combined") boson = "combined";
  c->Print((config.OutputFileDir + "Plots/"+ dirYear + "Profiles/"+"scale_factors_" +boson+"_"+kine+"_"+SF+".pdf").c_str());  
  delete hNom;   delete hSysUp;   delete hSysDown; delete hStatUp;   delete hStatDown; delete c;

  return;
}


void SFPlots::CombinedPlot(Config config, std::string kine, std::string boson){

  fnom    = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + nameOfSample + "_" + wzchannel + "_nominal" + puname +total+".root").c_str());  
  
  TH1D *hIso, *hID, *hTri, *hKF, *hPU;
  TProfile *hIsoAux, *hIDAux, *hTriAux, *hKFAux, *hPUAux;
  double GeV=1./1.;
  if(kine=="mwt"|| kine=="pt") GeV=1./1000.;

  hIsoAux  =  (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_Iso").c_str());
  hIDAux   =  (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_ID").c_str());
  hTriAux  =  (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_Tri").c_str());
  hKFAux   =  (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_KF").c_str());
  hPUAux   =  (TProfile*)fnom ->GetObjectUnchecked(("p_"+kine+"_PU").c_str());
  
  hIso  = scaleXaxis(hIsoAux, GeV, kine);
  hID   = scaleXaxis(hIDAux, GeV, kine);
  hTri  = scaleXaxis(hTriAux, GeV, kine);
  hKF   = scaleXaxis(hKFAux, GeV, kine);
  hPU   = scaleXaxis(hPUAux, GeV, kine);

  if(hIso==NULL) exit(10);
  
  std::string xlabel,ylabel; double xlow=0., xhigh=0., ylow=0., yhigh=0.;
  TString leglabel;

  if(boson=="zmumu") {xlabel="m_{Z} [GeV]"; xlow=60000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="wminus"|| boson=="wplus"){xlow=40000.*GeV; xhigh=7000000.*GeV;}
  if(boson=="zmumu")  leglabel="Z#rightarrow#mu#mu";
  if(boson=="wplus")  leglabel="W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus") leglabel="W^{-}#rightarrow#mu^{-}#nu";
  if(config.WZSelection=="combined") leglabel="MC W#rightarrow#mu#nu";
  ylabel="Scale Factor";

  if(kine=="mwt")xlabel="m_{T} [GeV]";
  if(kine=="pt")xlabel="p_{T} [GeV]";
  if(kine=="eta")xlabel="#eta^{#mu}";

  if( kine=="eta"){ylow=0.8; yhigh=1.1;}
  if((kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.175;}

  if( kine=="eta"){ylow=0.8; yhigh=1.175;}
  if( (kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.175;}

  if( kine=="eta"){ylow=0.8; yhigh=1.175;}
  if( (kine=="mwt" || kine=="pt")){ylow=0.8; yhigh=1.175;}

  //if(puFlag){ylow=0.2; yhigh=1.8;}

  plotAxisLine(hIso,kCyan,kCyan,20,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(hID,kCyan-6,kCyan-6,21,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(hTri,kBlue-3,kBlue-3,22,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(hKF,kViolet-6,kViolet-6,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(hPU,kViolet,kViolet,29,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,true,ylow,yhigh);

  TString lumi="";
  TCanvas *c = new TCanvas("canvas","canvas",550,450);
  TLegend* leg = new TLegend(0.6,0.825,0.8,0.67,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry(hIso,  "#bf{#scale[0.85]{Isolation}}","P");  
  leg->AddEntry(hID,  "#bf{#scale[0.85]{ID}}","P");  
  leg->AddEntry(hTri,  "#bf{#scale[0.85]{Trigger}}","P");  
  leg->AddEntry(hKF,  "#bf{#scale[0.85]{KFactor}}","P");  
  leg->AddEntry(hPU,  "#bf{#scale[0.85]{Pileup}}","P");  

  TPaveText *box;
  box = new TPaveText(0.625,0.825,0.8,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi+" Simulation");
  box->AddText(""+ leglabel+"");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  if(kine=="mwt" || kine=="pt" || kine=="m") pad->SetLogx();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.1125);

  hIso->SetTitle("");
  hIso->GetXaxis()->SetTitleOffset(1.45);  
  hIso->DrawCopy("pe");  
  hID ->DrawCopy("samese p");
  hTri->DrawCopy("samese p");
  hKF ->DrawCopy("samese p");
  hPU ->DrawCopy("samese p");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  if(config.WZSelection=="combined") boson = "combined";
  c->Print((config.OutputFileDir + "Plots/"+ dirYear + "Profiles/"+"all_scale_factors_"+kine+"_" +boson+".pdf").c_str());  
  delete hIso;  delete hID;   delete hTri; delete hKF; delete hPU; delete c;

  return;
}


TH1D* SFPlots::scaleXaxis(TProfile* hist, double scale, std::string kine)
{
  if(scale == 1) return hist;
  double binMin = hist->GetXaxis()->GetXmin();
  double binMax = hist->GetXaxis()->GetXmax();
  double nBins  = hist->GetXaxis()->GetNbins();
  double binMinScaled = binMin*scale;
  double binMaxScaled = binMax*scale;
  TH1D* scaled = NULL;

  if(kine=="mwt" || kine=="pt" || kine=="m"){ 
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



void SFPlots::plotAxisLine(TH1D* hist, int lineColor, int markerColor,
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
  if(xRange==true)
    hist->GetXaxis()->SetRangeUser(xlow,xhigh);
  if(yRange==true)
    hist->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}



void SFPlots::setstyle(){

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
