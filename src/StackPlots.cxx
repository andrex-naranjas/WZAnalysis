//author: Andres Ramirez <andres.ramirez.morales@cern.ch> June 2017
//@brief based on plotting functions from Alex Law <atlaw@ucsc.edu>
//update April 2018
//WZplots includes
#ifndef STACKPLOTS_CXX
#define STACKPLOTS_CXX

#include "StackPlots.h"
#include "WZplots.h"

#include "AtlasStyle.h"
#include "AtlasLabels.h"
#include "AtlasUtils.h"

#include <iostream>
#include <vector>
#include <string>


StackPlots::StackPlots()
{
}

StackPlots::~StackPlots(){}


void StackPlots::initialize(Config config){

  // kine.push_back("mu");
  // kine.push_back("mu");//ugly plotting bug
  //kine.push_back("npvtx");
  kine.push_back("rmass");
  kine.push_back("met");
  // kine.push_back("d0sig");
  kine.push_back("d0");
  kine.push_back("pt");
  kine.push_back("eta");
  //kine.push_back("iso30");
  kine.push_back("eta_d");
  kine.push_back("phi");
  kine.push_back("phi_met");
  kine.push_back("delta_phi");
  //kine.push_back("mu");


  // kine.push_back("pt");
  // kine.push_back("eta_d");
  // kine.push_back("phi");

  // kine.push_back("phi_met");
  // kine.push_back("met");
  // kine.push_back("rmass");

  // kine.push_back("SumET");
  // kine.push_back("met_recoil");
  // kine.push_back("rmass_recoil");

  // kine.push_back("pt20Calo");
  // kine.push_back("pt20Track");
  // kine.push_back("delta_phi");
  // kine.push_back("z0");
  // kine.push_back("d0");

  // kine.push_back("pt_boson_reco");
  // kine.push_back("pt_boson_reco_p");

  // kine.push_back("pt_boson_truth");
  // kine.push_back("mass_boson_truth");
  // kine.push_back("eta_boson_truth");

  if(config.Dod0Cut=="True") d0Name="";
  if(config.Dod0Cut!="True") d0Name="_nod0cut";  

  GetHistos(config);

  return;
}


void StackPlots::execute(Config config){

  std::vector<std::string> range, PN; 
  range.clear(); PN.clear();

  range.push_back("");//inclusive region
  //  if(config.OnlyInclusive!="True")
  range.push_back("_high");//high mass region

  if(config.WZSelection=="wplus")
    PN.push_back("P");//Positive Channel
  if(config.WZSelection=="wminus")
    PN.push_back("N");//Negative Channel
  if(config.WZSelection=="combined")
    PN.push_back("B");//Combined Channels
  if(config.WZSelection=="zmumu")//Z Channel
    PN.push_back("Z");

  FILE* fp = fopen("filename.txt","w");
  FILE* sysLog = fopen("syslog.txt","w");

  //ouput directory according year
  std::string dirYear="";
  if(config.DataYears=="2015+2016") dirYear="2015p2016/";
  if(config.DataYears=="2017") dirYear="2017/";
  if(config.DataYears=="2018") dirYear="2018/";
  if(config.DataYears=="Full") dirYear="Full/";

  for(int iKine=0;iKine<(int)kine.size();++iKine){
    for(int iRange=0;iRange<(int)range.size();++iRange){
      for(int iPN=0;iPN<(int)PN.size();++iPN){

	  std::vector<TH1F*> histList;
	  std::string nameOfFile = (config.OutputFileDir+"Plots/"+ dirYear + kine[iKine]+range[iRange]+ d0Name+"_"+PN[iPN]+".pdf").c_str();

	  if(config.PlotData) std::cout<<"   testing... "<<std::endl;

	  bool fit=false;

	  double scale = 1.;
	  if(kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil" || kine[iKine]=="SumET" || kine[iKine]=="met_recoil") scale = 1./1000.;

	  double minx=0., maxx=0.;
	  if((kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="met_recoil") && range[iRange]==""){minx=20; maxx=100;}
	  if((kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil")  && range[iRange]==""){minx=40; maxx=150;}
	  if((kine[iKine]=="pt" || kine[iKine]=="met") && range[iRange]=="_high"){minx=55; maxx=3000;}
	  if((kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil") && range[iRange]=="_high"){minx=110; maxx=6000;}
	  if(kine[iKine]=="eta_d"){minx=-2.6; maxx=2.6;}
	  if(kine[iKine]=="eta"){minx=0; maxx=2.6;}
	  if(kine[iKine]=="phi" || kine[iKine]=="phi_met"){minx=-3.141592; maxx=3.141592;}
	  if(kine[iKine]=="delta_phi"){minx=0.; maxx=3.141592;}
	  if(kine[iKine]=="mu" || kine[iKine]=="npvtx"){minx=0; maxx=100;}
	  if(kine[iKine]=="d0"){minx=-7.; maxx=7.;}
	  if(kine[iKine]=="z0"){minx=-0.5; maxx=0.5;}
	  if(kine[iKine]=="pt20Calo"){minx=0.; maxx=0.1;}
	  if(kine[iKine]=="pt20Track"){minx=0.; maxx=0.5;}
	  if(kine[iKine]=="SumET"){minx=0.; maxx=2000.;}

	  if(kine[iKine]=="rmass" && PN[iPN]=="Z"){minx=66.; maxx=116.;}

	  bool Ylog=false, Xlog=false;
	  bool userDefinedY=false; double miny=0., maxy=3.8e4;
	  if(range[iRange]=="" && (kine[iKine]=="eta_d" || kine[iKine]=="phi" || kine[iKine]=="phi_met"))
	    //{userDefinedY=true; miny=0.1; maxy=6.2e5;}
	    {userDefinedY=true; miny=0.1; maxy=0.6e7;/*maxy=8.1e6;*/}

	  if(range[iRange]=="" && kine[iKine]=="eta")
	    //{userDefinedY=true; miny=0.1; maxy=6.2e5;}
	    {userDefinedY=true; miny=0.1; maxy=1.15e8;/*maxy=8.1e6;*/}

	  if(range[iRange]=="_high"  && (kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="met_recoil" || kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil"))
	    // {Ylog=true; Xlog=true; userDefinedY=true; miny=0.1; maxy=1e4;}
	    {Ylog=true; Xlog=true; userDefinedY=true; miny=0.1; maxy=1e8;}


	  if(range[iRange]==""  && ( kine[iKine]=="rmass" ) && PN[iPN]=="Z")
	    // {Ylog=true; Xlog=true; userDefinedY=true; miny=0.1; maxy=1e4;}
	    {Ylog=true; Xlog=false; userDefinedY=true; miny=2500; maxy=1e9;}


	  if(range[iRange]=="_high" && (kine[iKine]=="eta_d" || kine[iKine]=="eta" || kine[iKine]=="phi" || kine[iKine]=="phi_met"))
	    //{Ylog=true; Xlog=false; userDefinedY=true; miny=0.1; maxy=5e6;}
	  {Ylog=true; Xlog=false; userDefinedY=true; miny=0.1; maxy=5e10;}

	  double ratioYmin = 0.775, ratioYmax = 1.225;
	  if(range[iRange]=="_high"){ratioYmin = 0.6; ratioYmax = 1.4;}
	  //if(range[iRange]=="_high"){ratioYmin = 0.09; ratioYmax = 2.01;}

	  double x1_legend=0.56, y1_legend=0.275;
	  if(kine[iKine]=="eta_d" || kine[iKine]=="eta" || kine[iKine]=="phi" || kine[iKine]=="phi_met" || ((range[iRange]=="_high" || PN[iPN]=="Z") && (kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="met_recoil" || kine[iKine]=="rmass")) )
	    {x1_legend=0.55; y1_legend=0.45;}
	  if(kine[iKine]=="d0"){x1_legend=0.72; y1_legend=0.45;}

	  double Xatlas=0.56, Yatlas=0.8;
	  if(kine[iKine]=="eta_d" || kine[iKine]=="eta" || kine[iKine]=="phi" || kine[iKine]=="phi_met" || ((range[iRange]=="_high" || PN[iPN]=="Z") && (kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="met_recoil" || kine[iKine]=="rmass")) )
	    {Xatlas=0.2; Yatlas=0.77;}	  
	  if(kine[iKine]=="d0"){Xatlas=0.16; Yatlas=0.8;}

	  if(kine[iKine]=="delta_phi"){Xatlas=0.16; Yatlas=0.8; x1_legend=0.2; y1_legend=0.275;}

	  TString units="";
	  if(kine[iKine]=="pt" || kine[iKine]=="met" || kine[iKine]=="met_recoil" || kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil") units="GeV";

	  std::string titleXaxis;
	  if(kine[iKine]=="pt") titleXaxis  = "p_{T}^{#mu} [GeV]";    if(kine[iKine]=="met" || kine[iKine]=="met_recoil" ) titleXaxis = "E_{T}^{miss} [GeV]";
	  if(kine[iKine]=="rmass" || kine[iKine]=="rmass_recoil") titleXaxis = "m_{T} [GeV]";          if(kine[iKine]=="eta_d") titleXaxis = "#eta^{#mu}";
          if(kine[iKine]=="eta") titleXaxis = "#||{#eta^{#mu}}";
	  if(kine[iKine]=="mu") titleXaxis = "<#mu>";                 if(kine[iKine]=="phi") titleXaxis = "#phi^{#mu} [rad]";
	  if(kine[iKine]=="phi_met") titleXaxis = "#phi^{E_{T}} [rad]";     if(kine[iKine]=="delta_phi") titleXaxis = "#||{#Delta#phi} [rad]";
	  if(kine[iKine]=="z0") titleXaxis = "z0sin#theta [mm]";      if(kine[iKine]=="d0") titleXaxis = "d0";
	  if(kine[iKine]=="nptvx")titleXaxis = "Npvtx";               if(kine[iKine]=="pt20Calo") titleXaxis = "ConeE_{T}/p_{T}^{#mu}";
	  if(kine[iKine]=="pt20Track") titleXaxis = "TrackE_{T}20/p_{T}^{#mu}"; 
	  if(kine[iKine]=="SumET") titleXaxis = "#sum E_{T}^{PFO} [GeV]";

	  if(kine[iKine]=="rmass" && PN[iPN]=="Z") titleXaxis = "m_{#mu#mu} [GeV]";
	 
	  std::string YlabelPrecision = "%.0f";// printf format " %.{N]f"
	  if(kine[iKine]=="eta_d" || kine[iKine]=="eta" || kine[iKine]=="phi" || kine[iKine]=="phi_met" || kine[iKine]=="d0"  || kine[iKine]=="delta_phi")
	    YlabelPrecision = "%.1f";
	  if(kine[iKine]=="z0" || kine[iKine]=="pt20Calo" || kine[iKine]=="pt20Track")
	    YlabelPrecision = "%.2f";

	  //Get histograms
	  TFile* inFile = new TFile((config.OutputFileDir+"kine/" + dirYear + "muon_"+kine[iKine]+range[iRange]+".root").c_str());
	  //multijet
	  TFile* inMultijet;
	  if(config.PlotMulti)
	    inMultijet = new TFile((config.OutputFileDir+"Files/" + dirYear + "/Multijet/multijet_final.root").c_str());
	  //TFile* inMultijet = new TFile("/data/morales/atlas/r21ControlPlots/Files/2017/Multijet/multijet_final.root");

	  if(inMultijet==NULL && config.PlotMulti) {std::cout<<"No input file (kine), bye!"<<std::endl; exit(10);}
	  if(inFile==NULL) {std::cout<<"No input file (kine), bye!"<<std::endl; exit(10);}
  
	  TH1F *h_wplus, *h_wminus;
	  if((PN[iPN]=="P" || PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWplusMunu)  h_wplus   = (TH1F*)(inFile->GetObjectUnchecked("wplusmunu"));
	  if((PN[iPN]=="N" || PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWminMunu )  h_wminus  = (TH1F*)(inFile->GetObjectUnchecked("wminmunu"));
	  if((PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWminMunu && config.PlotWplusMunu) h_wplus->Add(h_wminus); //combine W+ and W-
	  
	  TH1F *h_wptv, *h_wmtv;
	  if((PN[iPN]=="P" || PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWplusTaunu)  h_wptv  = (TH1F*)(inFile->GetObjectUnchecked("wplustaunu"));
	  if((PN[iPN]=="N" || PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWminTaunu )  h_wmtv  = (TH1F*)(inFile->GetObjectUnchecked("wmintaunu"));
	  if((PN[iPN]=="B" || PN[iPN]=="Z") && config.PlotWplusTaunu && config.PlotWminTaunu) h_wptv->Add(h_wmtv); //combine W+ and W-

	  TH1F *h_zmumu, *h_zmumuN;
	  if(PN[iPN]=="Z"  && config.PlotZmumu) h_zmumu = (TH1F*)(inFile->GetObjectUnchecked("zmumu"));
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZmumu) h_zmumu   = (TH1F*)(inFile->GetObjectUnchecked("zmumuP"));
	  if((PN[iPN]=="N" || PN[iPN]=="B") && config.PlotZmumu) h_zmumuN  = (TH1F*)(inFile->GetObjectUnchecked("zmumuN"));
	  if(PN[iPN]=="B"  && config.PlotZmumu) h_zmumu->Add(h_zmumuN);//combine W+ and W-
	  
	  TH1F *h_ztautau, *h_ztautauN;
	  if(PN[iPN]=="Z" && config.PlotZtautau)  h_ztautau = (TH1F*)(inFile->GetObjectUnchecked("ztautau"));
	  if(PN[iPN]=="P" || PN[iPN]=="B" && config.PlotZtautau)  h_ztautau = (TH1F*)(inFile->GetObjectUnchecked("ztautauP"));
	  if(PN[iPN]=="N" || PN[iPN]=="B" && config.PlotZtautau)  h_ztautauN = (TH1F*)(inFile->GetObjectUnchecked("ztautauN"));
	  if(PN[iPN]=="B" && config.PlotZtautau) h_ztautau->Add(h_ztautauN);//combine W+ and W-

	  //Top backgrounds
	  TH1F *h_top_1, *h_top_2, *h_top_3, *h_top_4, *h_top_5;
	  if(PN[iPN]=="Z"){
	    if(config.PlotTtbar)           h_top_1   = (TH1F*)(inFile->GetObjectUnchecked("ttbarZ"));
	    if(config.PlotWt_top)         {h_top_2   = (TH1F*)(inFile->GetObjectUnchecked("wt_topZ"));         h_top_1->Add(h_top_2);}
	    if(config.PlotWt_antitop)     {h_top_3   = (TH1F*)(inFile->GetObjectUnchecked("wt_antitopZ"));     h_top_1->Add(h_top_3);}
	    if(config.PlotSingle_top)     {h_top_4   = (TH1F*)(inFile->GetObjectUnchecked("single_topZ"));     h_top_1->Add(h_top_4);}
	    // if(config.PlotSingle_antitop) {h_top_5   = (TH1F*)(inFile->GetObjectUnchecked("single_antitopP")); h_top_1->Add(h_top_5);}	         
	  }

	  if(PN[iPN]=="P" || PN[iPN]=="B"){
	    if(config.PlotTtbar)           h_top_1   = (TH1F*)(inFile->GetObjectUnchecked("ttbarP"));
	    if(config.PlotWt_top)         {h_top_2   = (TH1F*)(inFile->GetObjectUnchecked("wt_topP"));         h_top_1->Add(h_top_2);}
	    if(config.PlotWt_antitop)     {h_top_3   = (TH1F*)(inFile->GetObjectUnchecked("wt_antitopP"));     h_top_1->Add(h_top_3);}
	    if(config.PlotSingle_top)     {h_top_4   = (TH1F*)(inFile->GetObjectUnchecked("single_topP"));     h_top_1->Add(h_top_4);}
	    // if(config.PlotSingle_antitop) {h_top_5   = (TH1F*)(inFile->GetObjectUnchecked("single_antitopP")); h_top_1->Add(h_top_5);}	         
	  }

	  TH1F *h_N_top_1, *h_N_top_2, *h_N_top_3, *h_N_top_4, *h_N_top_5;
	  if(PN[iPN]=="N" || PN[iPN]=="B"){
	    if(config.PlotTtbar)           h_N_top_1   = (TH1F*)(inFile->GetObjectUnchecked("ttbarN"));
	    if(config.PlotWt_top)         {h_N_top_2   = (TH1F*)(inFile->GetObjectUnchecked("wt_topN"));         h_N_top_1->Add(h_N_top_2);}
	    if(config.PlotWt_antitop)     {h_N_top_3   = (TH1F*)(inFile->GetObjectUnchecked("wt_antitopN"));     h_N_top_1->Add(h_N_top_3);}
	    if(config.PlotSingle_top)     {h_N_top_4   = (TH1F*)(inFile->GetObjectUnchecked("single_topN"));     h_N_top_1->Add(h_N_top_4);}
	    // if(config.PlotSingle_antitop) {h_N_top_5   = (TH1F*)(inFile->GetObjectUnchecked("single_antitopN")); h_N_top_1->Add(h_N_top_5);}
	  }
	  if(PN[iPN]=="B" && config.PlotTtbar) h_top_1->Add(h_N_top_1);//combine W+ and W-

	  //Diboson backgrounds
	  TH1F *h_dibo_1, *h_dibo_2, *h_dibo_3, *h_dibo_4, *h_dibo_5, *h_dibo_6, *h_dibo_7, *h_dibo_8, *h_dibo_9;
	  if(PN[iPN]=="Z"){
	    if(config.PlotZzqqll)     h_dibo_1   = (TH1F*)(inFile->GetObjectUnchecked("zzqqllZ"));          			
	    if(config.PlotWwqqll)    {h_dibo_2   = (TH1F*)(inFile->GetObjectUnchecked("wwqqllZ"));        h_dibo_1->Add(h_dibo_2);} 
	    if(config.PlotWwpqqmlnu) {h_dibo_3   = (TH1F*)(inFile->GetObjectUnchecked("wwpqqmlnuZ"));     h_dibo_1->Add(h_dibo_3);} 
	    if(config.PlotWwplnumqq) {h_dibo_4   = (TH1F*)(inFile->GetObjectUnchecked("wwplnumqqZ"));     h_dibo_1->Add(h_dibo_4);} 
	    if(config.PlotWzlnuqq)   {h_dibo_5   = (TH1F*)(inFile->GetObjectUnchecked("wzlnuqqZ"));       h_dibo_1->Add(h_dibo_5);} 
	    if(config.PlotZzllll)    {h_dibo_6   = (TH1F*)(inFile->GetObjectUnchecked("zzllllZ"));        h_dibo_1->Add(h_dibo_6);} 
	    if(config.PlotWzlnull)   {h_dibo_7   = (TH1F*)(inFile->GetObjectUnchecked("wzlnullZ"));       h_dibo_1->Add(h_dibo_7);} 
	    if(config.PlotWzlnununu) {h_dibo_8   = (TH1F*)(inFile->GetObjectUnchecked("wzlnununuZ"));     h_dibo_1->Add(h_dibo_8);} 
	    if(config.PlotZzllnunu)  {h_dibo_9   = (TH1F*)(inFile->GetObjectUnchecked("zzllnunuZ"));      h_dibo_1->Add(h_dibo_9);} 
	    }

	  if(PN[iPN]=="P" || PN[iPN]=="B"){
	    if(config.PlotZzqqll)     h_dibo_1   = (TH1F*)(inFile->GetObjectUnchecked("zzqqllP"));          			
	    if(config.PlotWwqqll)    {h_dibo_2   = (TH1F*)(inFile->GetObjectUnchecked("wwqqllP"));        h_dibo_1->Add(h_dibo_2);} 
	    if(config.PlotWwpqqmlnu) {h_dibo_3   = (TH1F*)(inFile->GetObjectUnchecked("wwpqqmlnuP"));     h_dibo_1->Add(h_dibo_3);} 
	    if(config.PlotWwplnumqq) {h_dibo_4   = (TH1F*)(inFile->GetObjectUnchecked("wwplnumqqP"));     h_dibo_1->Add(h_dibo_4);} 
	    if(config.PlotWzlnuqq)   {h_dibo_5   = (TH1F*)(inFile->GetObjectUnchecked("wzlnuqqP"));       h_dibo_1->Add(h_dibo_5);} 
	    if(config.PlotZzllll)    {h_dibo_6   = (TH1F*)(inFile->GetObjectUnchecked("zzllllP"));        h_dibo_1->Add(h_dibo_6);} 
	    if(config.PlotWzlnull)   {h_dibo_7   = (TH1F*)(inFile->GetObjectUnchecked("wzlnullP"));       h_dibo_1->Add(h_dibo_7);} 
	    if(config.PlotWzlnununu) {h_dibo_8   = (TH1F*)(inFile->GetObjectUnchecked("wzlnununuP"));     h_dibo_1->Add(h_dibo_8);} 
	    if(config.PlotZzllnunu)  {h_dibo_9   = (TH1F*)(inFile->GetObjectUnchecked("zzllnunuP"));      h_dibo_1->Add(h_dibo_9);} 
	    }

	  if(h_dibo_1==NULL) exit(10);

	  TH1F *h_N_dibo_1, *h_N_dibo_2, *h_N_dibo_3, *h_N_dibo_4, *h_N_dibo_5, *h_N_dibo_6, *h_N_dibo_7, *h_N_dibo_8, *h_N_dibo_9;
	  if(PN[iPN]=="N" || PN[iPN]=="B"){
	   if(config.PlotZzqqll)       h_N_dibo_1   = (TH1F*)(inFile->GetObjectUnchecked("zzqqllN"));
	   if(config.PlotWwqqll)      {h_N_dibo_2   = (TH1F*)(inFile->GetObjectUnchecked("wwqqllN"));        h_N_dibo_1->Add(h_N_dibo_2);}  
	   if(config.PlotWwpqqmlnu)   {h_N_dibo_3   = (TH1F*)(inFile->GetObjectUnchecked("wwpqqmlnuN"));     h_N_dibo_1->Add(h_N_dibo_3);} 
	   if(config.PlotWwplnumqq)   {h_N_dibo_4   = (TH1F*)(inFile->GetObjectUnchecked("wwplnumqqN"));     h_N_dibo_1->Add(h_N_dibo_4);} 
	   if(config.PlotWzlnuqq)     {h_N_dibo_5   = (TH1F*)(inFile->GetObjectUnchecked("wzlnuqqN"));       h_N_dibo_1->Add(h_N_dibo_5);} 
	   if(config.PlotZzllll)      {h_N_dibo_6   = (TH1F*)(inFile->GetObjectUnchecked("zzllllN"));        h_N_dibo_1->Add(h_N_dibo_6);}  
	   if(config.PlotWzlnull)     {h_N_dibo_7   = (TH1F*)(inFile->GetObjectUnchecked("wzlnullN"));       h_N_dibo_1->Add(h_N_dibo_7);}  
	   if(config.PlotWzlnununu)   {h_N_dibo_8   = (TH1F*)(inFile->GetObjectUnchecked("wzlnununuN"));     h_N_dibo_1->Add(h_N_dibo_8);}  
	   if(config.PlotZzllnunu)    {h_N_dibo_9   = (TH1F*)(inFile->GetObjectUnchecked("zzllnunuN"));      h_N_dibo_1->Add(h_N_dibo_9);}  
	  }
	  if(PN[iPN]=="B" && config.PlotZzqqll) h_dibo_1->Add(h_N_dibo_1);//combine W+ and W-

	  //multijet background
	  TH1F *h_N_multi,*h_multi;
	  if(config.PlotMulti){
	    if(kine[iKine]=="d0"){
	      h_multi = (TH1F*)inMultijet->GetObjectUnchecked("d0sig_wplus");	     
	      h_N_multi = (TH1F*)inMultijet->GetObjectUnchecked("d0sig_wminus");
	    }else{
	      h_multi = (TH1F*)inMultijet->GetObjectUnchecked((kine[iKine]+"_wplus").c_str());	     
	      h_N_multi = (TH1F*)inMultijet->GetObjectUnchecked((kine[iKine]+"_wminus").c_str());
	    }
	  }//plotmulti

	  if(h_multi==NULL && config.PlotMulti){ std::cout<<"Multijet enabled, but not found in file"<<std::endl; exit(1);}
	  
	  //Data
	  TH1F *data=NULL;
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotData) data = (TH1F*)(inFile->GetObjectUnchecked("dataP"));
	  if(PN[iPN]=="N" && config.PlotData)                   data = (TH1F*)(inFile->GetObjectUnchecked("dataN"));
	  if(PN[iPN]=="B" && config.PlotData)                   data->Add((TH1F*)(inFile->GetObjectUnchecked("dataN")));//combine W+ and W-
	  if(PN[iPN]=="Z" && config.PlotData)                   data = (TH1F*)(inFile->GetObjectUnchecked("dataZ"));

	  //scale
	  // if(PN[iPN]=="P"){ h_wplus->Scale(1.065); h_zmumu->Scale(1.065); h_top_1->Scale(1.065); h_ztautau->Scale(1.065); h_dibo_1->Scale(1.065);}
	  // if(PN[iPN]=="N"){ h_wminus->Scale(1.045); h_zmumuN->Scale(1.045); h_N_top_1->Scale(1.045); h_ztautauN->Scale(1.045); h_N_dibo_1->Scale(1.045);}
	  // if(config.PlotMulti){
	  //   h_multi ->Scale(1.065);
	  //   h_N_multi->Scale(1.045);
	  // }

	  if(kine[iKine]=="eta"){

	    for(int k=1;k<data->GetXaxis()->GetNbins()+1;k++)
	      data->SetBinContent(k,(data   ->GetBinContent(k)/data   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZzqqll && kine[iKine]=="eta")
	      for(int k=1;k<h_dibo_1->GetXaxis()->GetNbins()+1;k++)
		h_dibo_1->SetBinContent(k,(h_dibo_1   ->GetBinContent(k)/h_dibo_1   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N"  && config.PlotZzqqll && kine[iKine]=="eta")
	      for(int k=1;k<h_N_dibo_1->GetXaxis()->GetNbins()+1;k++)
		h_N_dibo_1->SetBinContent(k,(h_N_dibo_1   ->GetBinContent(k)/h_N_dibo_1   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZtautau && kine[iKine]=="eta")
	      for(int k=1;k<h_ztautau->GetXaxis()->GetNbins()+1;k++)
		h_ztautau->SetBinContent(k,(h_ztautau   ->GetBinContent(k)/h_ztautau   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N" && config.PlotZtautau && kine[iKine]=="eta")
	      for(int k=1;k<h_ztautauN->GetXaxis()->GetNbins()+1;k++)
		h_ztautauN->SetBinContent(k,(h_ztautauN   ->GetBinContent(k)/h_ztautauN   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotWplusTaunu && kine[iKine]=="eta")
	      for(int k=1;k<h_wptv->GetXaxis()->GetNbins()+1;k++)
		h_wptv->SetBinContent(k,(h_wptv   ->GetBinContent(k)/h_wptv   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N" && config.PlotWminTaunu && kine[iKine]=="eta")
	      for(int k=1;k<h_wmtv->GetXaxis()->GetNbins()+1;k++)
		h_wmtv->SetBinContent(k,(h_wmtv   ->GetBinContent(k)/h_wmtv   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZmumu && kine[iKine]=="eta")
	      for(int k=1;k<h_zmumu->GetXaxis()->GetNbins()+1;k++)
		h_zmumu->SetBinContent(k,(h_zmumu   ->GetBinContent(k)/h_zmumu   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N" && config.PlotZmumu && kine[iKine]=="eta")
	      for(int k=1;k<h_zmumuN->GetXaxis()->GetNbins()+1;k++)
		h_zmumuN->SetBinContent(k,(h_zmumuN   ->GetBinContent(k)/h_zmumuN   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotTtbar && kine[iKine]=="eta")
	      for(int k=1;k<h_top_1->GetXaxis()->GetNbins()+1;k++)
		h_top_1->SetBinContent(k,(h_top_1   ->GetBinContent(k)/h_top_1   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N" && config.PlotTtbar && kine[iKine]=="eta")
	      for(int k=1;k<h_N_top_1->GetXaxis()->GetNbins()+1;k++)
		h_N_top_1->SetBinContent(k,(h_N_top_1   ->GetBinContent(k)/h_N_top_1   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotMulti && kine[iKine]=="eta")
	      for(int k=1;k<h_multi->GetXaxis()->GetNbins()+1;k++)
		h_multi->SetBinContent(k,(h_multi   ->GetBinContent(k)/h_multi   ->GetBinWidth(k)));
	    
	    if(PN[iPN]=="N" && config.PlotMulti && kine[iKine]=="eta")
	      for(int k=1;k<h_N_multi->GetXaxis()->GetNbins()+1;k++)
		h_N_multi->SetBinContent(k,(h_N_multi   ->GetBinContent(k)/h_N_multi   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotWplusMunu && kine[iKine]=="eta")
	      for(int k=1;k<h_wplus->GetXaxis()->GetNbins()+1;k++)
		h_wplus->SetBinContent(k,(h_wplus   ->GetBinContent(k)/h_wplus   ->GetBinWidth(k)));
	    
	    if((PN[iPN]=="N") && config.PlotWminMunu && kine[iKine]=="eta")
	      for(int k=1;k<h_wminus->GetXaxis()->GetNbins()+1;k++)
		h_wminus->SetBinContent(k,(h_wminus   ->GetBinContent(k)/h_wminus   ->GetBinWidth(k)));
	  }
	 
	  //std::cout<<h_sum_error->GetBinError(k)<<"   "<<h_sys->GetBinContent(k)<<std::endl;//check consistency of stats errors!! it matches :)
	  //h_sum_error->SetBinError(k,h_sys->GetBinContent(k));
	  //Systematics
	  TH1F *h_sum_error=NULL, *h_sys=NULL;
	  if(config.SysInplots){
	    double sumSysP=0.;
	    if(PN[iPN]=="P" || PN[iPN]=="B"){ 
	      h_sum_error  = (TH1F*)(inFile->GetObjectUnchecked("sumP"));
	      h_sys = (TH1F*)(inFile->GetObjectUnchecked("sysupP"));

	      if(kine[iKine]=="eta"){
		for(int k=1;k<h_sum_error->GetXaxis()->GetNbins()+1; k++){
		  h_sum_error->SetBinContent(k,h_sum_error->GetBinContent(k)/h_sum_error->GetBinWidth(k));
		  h_sys->SetBinContent(k,h_sys->GetBinContent(k)/h_sys->GetBinWidth(k));		  
		}
	      }

	      //h_sum_error->Scale(1.065);
	      h_sum_error->Add(h_dibo_1);
	      h_sum_error->Add(h_zmumu);
	      h_sum_error->Add(h_ztautau);
	      h_sum_error->Add(h_wptv);
	      h_sum_error->Add(h_top_1);
	      if(config.PlotMulti) h_sum_error->Add(h_multi);
	      for(int k=1; k< h_sum_error->GetNbinsX()+1; k++){
		double error = std::sqrt(std::pow(h_sys->GetBinContent(k),2)+std::pow(h_sum_error->GetBinError(k),2));
		//sumSysP+=h_sys->GetBinContent(k);
		//error = h_sys->GetBinContent(k);
		h_sum_error->SetBinError(k,error);		
	      }
	      sumSysP=h_sys->Integral(0,h_sys->GetNbinsX()+1);
	    }//positive channel

	    double sumSysN=0.;	    	    
	    if(PN[iPN]=="N" || PN[iPN]=="B"){
	      h_sum_error  = (TH1F*)(inFile->GetObjectUnchecked("sumN"));
	      h_sys = (TH1F*)(inFile->GetObjectUnchecked("sysupN"));

	      if(kine[iKine]=="eta"){
		for(int k=1;k<h_sum_error->GetXaxis()->GetNbins()+1; k++){
		  h_sum_error->SetBinContent(k,h_sum_error->GetBinContent(k)/h_sum_error->GetBinWidth(k));
		  h_sys->SetBinContent(k,h_sys->GetBinContent(k)/h_sys->GetBinWidth(k));		  
		}
	      }

	      //h_sum_error->Scale(1.045);
	      h_sum_error->Add(h_N_dibo_1);
	      h_sum_error->Add(h_zmumuN);
	      h_sum_error->Add(h_ztautauN);
	      h_sum_error->Add(h_wmtv);
	      h_sum_error->Add(h_N_top_1);
	      if(config.PlotMulti) h_sum_error->Add(h_N_multi);	     
	      for(int k=1; k< h_sum_error->GetNbinsX()+1; k++){
		double error = std::sqrt(std::pow(h_sys->GetBinContent(k),2)+std::pow(h_sum_error->GetBinError(k),2));
		if(h_sum_error->GetBinContent(k)!=0)
		  sumSysN+=std::pow(h_sys->GetBinContent(k) ,2);
		//error = h_sys->GetBinContent(k);
		h_sum_error->SetBinError(k,error);
	      }
	      sumSysN=std::sqrt(sumSysN);

	      //sumSysN=std::sqrt(sumSysN);	
	    }//negative channel	

	    std::string nameva=kine[iKine];
	    if(range[iRange]==""){    
	      //fprintf(sysLog,"*********************************************\n");
	      fprintf(sysLog,"variable: %s  total sys minus: %g  total plus: %g \n",nameva.c_str(),sumSysN,sumSysP);
	    }
	    sumSysP=0; sumSysN=0;
	  }//if sys plots
         
	  //Instantiation of the WZplots class.
	  WZplots *m_wzplots = new WZplots();

	  TString latexTitle;
	  if(PN[iPN]=="P") latexTitle =  m_wzplots->getProcessLegendEntry("wplus_munu");
	  if(PN[iPN]=="N") latexTitle =  m_wzplots->getProcessLegendEntry("wminus_munu");
	  if(PN[iPN]=="B") latexTitle =  m_wzplots->getProcessLegendEntry("w_munu");
	  if(PN[iPN]=="Z") latexTitle =  m_wzplots->getProcessLegendEntry("z_mumu");

	  //The most important bit, set correctly the parameters for each plot
	  m_wzplots->SetPlotParameters(config,nameOfFile,scale,minx,maxx,Ylog,Xlog,ratioYmin,
				       ratioYmax, userDefinedY, miny,maxy, x1_legend, 
				       y1_legend,Xatlas, Yatlas, latexTitle, units,
				       titleXaxis, YlabelPrecision, h_sum_error, data, fit);
	  //We use addHistToList() and addSignalToStack() to build the stacked plot

	  if(h_top_1==NULL){
	    std::cout<<kine[iKine]<<"     range:  "<<range[iRange]<<std::endl;
	    exit(10);
	  }

	  std::string partLeg="";
	  if(PN[iPN]=="B")partLeg="";
	  if(PN[iPN]=="P")partLeg="plus";
	  if(PN[iPN]=="N")partLeg="minus";

	  //WBoson plots
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZzqqll) m_wzplots->addHistToList(&histList, h_dibo_1,   "diboson",scale);
	  if(PN[iPN]=="N"  && config.PlotZzqqll) m_wzplots->addHistToList(&histList, h_N_dibo_1,   "diboson",scale);
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZtautau) m_wzplots->addHistToList(&histList, h_ztautau,"z_tautau",scale);
	  if(PN[iPN]=="N" && config.PlotZtautau) m_wzplots->addHistToList(&histList, h_ztautauN,"z_tautau",scale);
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotWplusTaunu)m_wzplots->addHistToList(&histList, h_wptv,  "w"+partLeg+"_taunu",scale);
	  if(PN[iPN]=="N" && config.PlotWminTaunu )m_wzplots->addHistToList(&histList, h_wmtv,  "w"+partLeg+"_taunu",scale);
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotZmumu) m_wzplots->addHistToList(&histList,h_zmumu,"z_mumu",scale);
	  if(PN[iPN]=="N" && config.PlotZmumu) m_wzplots->addHistToList(&histList,h_zmumuN,"z_mumu",scale);
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotTtbar) m_wzplots->addHistToList(&histList, h_top_1, "top",scale);
	  if(PN[iPN]=="N" && config.PlotTtbar) m_wzplots->addHistToList(&histList, h_N_top_1, "top",scale);

	  if(config.PlotMulti){
	    if((PN[iPN]=="P" || PN[iPN]=="B") ) m_wzplots->addHistToList(&histList, h_multi, "multijet",scale);
	    if(PN[iPN]=="N" ) m_wzplots->addHistToList(&histList, h_N_multi, "multijet",scale);
	  }
	  
	  if((PN[iPN]=="P" || PN[iPN]=="B") && config.PlotWplusMunu) m_wzplots->addSignalToStack(&histList, h_wplus, "w"+partLeg+"_munu",scale);//top of stack
	  if((PN[iPN]=="N") && config.PlotWminMunu) m_wzplots->addSignalToStack(&histList, h_wminus, "w"+partLeg+"_munu",scale);// top of stack


	  //ZBoson plots
	  if(PN[iPN]=="Z" && config.PlotWplusMunu) m_wzplots->addHistToList(&histList, h_wplus, "w"+partLeg+"_munu",scale);
	  if(PN[iPN]=="Z" && config.PlotWplusTaunu)m_wzplots->addHistToList(&histList, h_wptv,  "w"+partLeg+"_taunu",scale);
	  if(PN[iPN]=="Z" && config.PlotTtbar) m_wzplots->addHistToList(&histList, h_top_1, "top",scale);
	  if(PN[iPN]=="Z" && config.PlotZtautau) m_wzplots->addHistToList(&histList, h_ztautau,"z_tautau",scale);
	  if(PN[iPN]=="Z" && config.PlotZzqqll) m_wzplots->addHistToList(&histList, h_dibo_1,   "diboson",scale);
	  if(PN[iPN]=="Z" && config.PlotZmumu) m_wzplots->addSignalToStack(&histList, h_zmumu, "z_mumu",scale);// top of stack

	  SetAtlasStyle();//in AtlasStyle.cxx
	  m_wzplots->SetControlPlotStyle();
	  m_wzplots->printPlot(fp,kine[iKine]+range[iRange]+PN[iPN]);
	  std::cout<<"NEXT PLOT!"<<std::endl;
	  delete m_wzplots;	  
	}//channel loop
      }//range loop
    }//kine loop
  fclose(fp);
  fclose(sysLog);

  return;
}


void StackPlots::GetHistos(Config config){

  //names for output files
  std::string wzchannel;
  std::string calib;
  if(config.SETCalibration=="True" && config.InsituCorrection=="True"){
    calib="_set_insitu";
  }else if(config.SETCalibration=="True" && config.InsituCorrection=="False"){
    calib="_set";
    }else if(config.SETCalibration=="False" && config.InsituCorrection=="True"){
    calib="_insitu";
  }else if(config.SETCalibration=="False" && config.InsituCorrection=="False"){
    calib="";
  }
  std::string puname="";
  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
  std::string systematic;
  systematic="_nominal";
  
  std::string dirInclusive, total;
  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}

  //ouput directory according year
  std::string dirYear="", year="";
  if(config.DataYears=="2015+2016"){ dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017"){ dirYear="2017/"; year="17";}
  if(config.DataYears=="2018"){ dirYear="2018/"; year="18";}
  if(config.DataYears=="Full"){ dirYear="Full/"; year="Full";}


  if(config.WZSelection=="wplus" || config.WZSelection=="combined"){
    wzchannel="wplus";
    if(config.PlotZmumu)      grzmumuP     = new TFile((config.OutputFileDir+"Files/" + dirYear + dirInclusive+"zmumu_"      + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotZtautau)    grztautauP   = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"ztautau_"    + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWplusTaunu) grwplustaunu = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wplustaunu_" + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWplusMunu)  grwplusmunu  = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wplusmunu_"  + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());

    std::cout<<(config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wplusmunu_"  + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str()<<std::endl;
    std::cin.get();

    if(config.PlotTtbar)          grttbarP          = new TFile((config.OutputFileDir+"Files/" + dirYear + "ttbar_"          + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_top)         grwt_topP         = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_top_"         + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_antitop)     grwt_antitopP     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_antitop_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_top)     grsingle_topP     = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_top_s_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_antitop) grsingle_antitopP = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_antitop_" + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

    if(config.PlotZzqqll)          grzzqqllP     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwqqll)          grwzqqllP     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwpqqmlnu)       grwwpqqmlnuP  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwpqqmlnu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwplnumqq)       grwwplnumqqP  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwplnumqq_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnuqq)         grwzlnuqqP    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnuqq_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllll)          grzzllllP     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnull)         grwzlnullP    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnull_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnununu)       grwzlnununuP  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnununu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllnunu)	   grzzllnunuP   = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllnunu_"   + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

    grDataP = new TFile((config.OutputFileDir+"Files/" + dirYear + "data"+year+"_wplus"+d0Name+".root").c_str());
    
    if(config.SysInplots)
      grSysP = new TFile((config.OutputFileDir+"Sys/" + dirYear + dirInclusive + "sys_wminus" +total+".root").c_str());
  }
  
  if(config.WZSelection=="wminus" || config.WZSelection=="combined"){
    wzchannel="wminus";
    if(config.PlotZmumu)     grzmumuN    = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"zmumu_"     + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotZtautau)   grztautauN  = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"ztautau_"   + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWminTaunu) grwmintaunu = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wmintaunu_" + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWminMunu)  grwminmunu  = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wminmunu_"  + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());

    if(config.PlotZzqqll)         grzzqqllN     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwqqll)         grwzqqllN     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwpqqmlnu)      grwwpqqmlnuN  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwpqqmlnu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwplnumqq)      grwwplnumqqN  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwplnumqq_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnuqq)        grwzlnuqqN    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnuqq_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllll)         grzzllllN     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnull)        grwzlnullN    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnull_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnununu)      grwzlnununuN  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnununu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllnunu)       grzzllnunuN   = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllnunu_"   + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

    if(config.PlotTtbar)          grttbarN          = new TFile((config.OutputFileDir+"Files/" + dirYear + "ttbar_"          + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_top)         grwt_topN         = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_top_"         + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_antitop)     grwt_antitopN     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_antitop_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_top)     grsingle_topN     = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_top_s_"   + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_antitop) grsingle_antitopN = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_antitop_" + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

  
    grDataN = new TFile((config.OutputFileDir+"Files/" + dirYear + "data"+year+"_wminus"+d0Name+".root").c_str());

    if(config.SysInplots)
      grSysN = new TFile((config.OutputFileDir+"Sys/" + dirYear +dirInclusive + "sys_wminus"+total+".root").c_str());
  }
 
  if(config.WZSelection=="zmumu"){
    wzchannel="z";
    if(config.PlotZmumu)      grzmumuZ     = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"zmumu_"      + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotZtautau)    grztautauZ   = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"ztautau_"    + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWplusTaunu) grwplustaunu = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wplustaunu_" + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWplusMunu)  grwplusmunu  = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wplusmunu_"  + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWminTaunu)  grwmintaunu  = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wmintaunu_"  + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());
    if(config.PlotWminMunu)   grwminmunu   = new TFile((config.OutputFileDir+"Files/" + dirYear +dirInclusive+"wminmunu_"   + wzchannel + systematic + calib + puname + d0Name +total+".root").c_str());

    if(config.PlotTtbar)          grttbarZ          = new TFile((config.OutputFileDir+"Files/" + dirYear + "ttbar_"          + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_top)         grwt_topZ         = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_top_"         + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWt_antitop)     grwt_antitopZ     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wt_antitop_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_top)     grsingle_topZ     = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_top_s_"   + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotSingle_antitop) grsingle_antitopZ = new TFile((config.OutputFileDir+"Files/" + dirYear + "single_antitop_" + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

    if(config.PlotZzqqll)          grzzqqllZ     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwqqll)          grwzqqllZ     = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzqqll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwpqqmlnu)       grwwpqqmlnuZ  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwpqqmlnu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWwplnumqq)       grwwplnumqqZ  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wwplnumqq_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnuqq)         grwzlnuqqZ    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnuqq_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllll)          grzzllllZ     = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllll_"     + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnull)         grwzlnullZ    = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnull_"    + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotWzlnununu)       grwzlnununuZ  = new TFile((config.OutputFileDir+"Files/" + dirYear + "wzlnununu_"  + wzchannel + systematic + calib + puname + d0Name +".root").c_str());
    if(config.PlotZzllnunu)	   grzzllnunuZ   = new TFile((config.OutputFileDir+"Files/" + dirYear + "zzllnunu_"   + wzchannel + systematic + calib + puname + d0Name +".root").c_str());

    grDataZ = new TFile((config.OutputFileDir+"Files/" + dirYear + "data"+year+"_z.root").c_str());
    
    if(config.SysInplots)
      grSysZ = new TFile((config.OutputFileDir+"Sys/" + dirYear + dirInclusive + "sys_z.root").c_str());
  }


  std::vector<std::string> range = {"","_high"};

  for(auto ran : range){
    for(Int_t iKine = 0; iKine < (int)kine.size(); ++iKine){
      
      TFile newfile((config.OutputFileDir + "kine/" + dirYear + "muon_" + kine[iKine] + ran +".root").c_str(),"RECREATE");
      
      TH1D *h_zmumuP, *h_ztautauP, *h_ttbarP, *h_wplustaunu, *h_wplusmunu, *sumP, *h_DataP;
      TH1D *h_wt_topP, *h_wt_antitopP, *h_single_topP, *h_single_antitopP;
      TH1D *h_zzqqllP, *h_wzqqllP, *h_wwpqqmlnuP, *h_wwplnumqqP, *h_wzlnuqqP, *h_zzllllP, *h_wzlnullP, *h_wzlnununuP, *h_zzllnunuP;
      TH1D *h_sysDownP, *h_sysDownN;

      TH1D *h_zmumuN, *h_ztautauN, *h_ttbarN, *h_wmintaunu, *h_wminmunu, *sumN, *h_DataN;
      TH1D *h_wt_topN, *h_wt_antitopN, *h_single_topN, *h_single_antitopN;
      TH1D *h_zzqqllN, *h_wzqqllN, *h_wwpqqmlnuN, *h_wwplnumqqN, *h_wzlnuqqN, *h_zzllllN, *h_wzlnullN, *h_wzlnununuN, *h_zzllnunuN;
      TH1D *h_sysUpP, *h_sysUpN;

      TH1D *h_zmumuZ, *h_ztautauZ, *h_ttbarZ, *sumZ, *h_DataZ;
      TH1D *h_wt_topZ, *h_wt_antitopZ, *h_single_topZ, *h_single_antitopZ;
      TH1D *h_zzqqllZ, *h_wzqqllZ, *h_wwpqqmlnuZ, *h_wwplnumqqZ, *h_wzlnuqqZ, *h_zzllllZ, *h_wzlnullZ, *h_wzlnununuZ, *h_zzllnunuZ;
      TH1D *h_sysDownZ, *h_sysUpZ;

      if(config.WZSelection=="wplus" || config.WZSelection=="combined"){
	
	if(config.PlotZmumu)	  h_zmumuP      = (TH1D*)grzmumuP    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZtautau)    h_ztautauP    = (TH1D*)grztautauP  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWplusTaunu) h_wplustaunu  = (TH1D*)grwplustaunu->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWplusMunu)  h_wplusmunu   = (TH1D*)grwplusmunu ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWplusMunu)  sumP = (TH1D*)h_wplusmunu->Clone("sumP");
 
	h_DataP = (TH1D*)grDataP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotTtbar)          h_ttbarP           = (TH1D*)grttbarP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());	
	if(config.PlotWt_top)         h_wt_topP          = (TH1D*)grwt_topP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWt_antitop)     h_wt_antitopP      = (TH1D*)grwt_antitopP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_top)     h_single_topP      = (TH1D*)grsingle_topP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_antitop) h_single_antitopP  = (TH1D*)grsingle_antitopP->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotZzqqll)         h_zzqqllP     = (TH1D*)grzzqqllP    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwqqll)         h_wzqqllP     = (TH1D*)grwzqqllP    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwpqqmlnu)      h_wwpqqmlnuP  = (TH1D*)grwwpqqmlnuP ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwplnumqq)      h_wwplnumqqP  = (TH1D*)grwwplnumqqP ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnuqq)        h_wzlnuqqP    = (TH1D*)grwzlnuqqP   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllll)         h_zzllllP     = (TH1D*)grzzllllP    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnull)        h_wzlnullP    = (TH1D*)grwzlnullP   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnununu)      h_wzlnununuP  = (TH1D*)grwzlnununuP ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllnunu)	      h_zzllnunuP   = (TH1D*)grzzllnunuP  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.SysInplots) h_sysDownP    = (TH1D*)grSysP       ->GetObjectUnchecked(("h_"+kine[iKine]+"_down").c_str());
	if(config.SysInplots) h_sysUpP      = (TH1D*)grSysP       ->GetObjectUnchecked(("h_"+kine[iKine]+"_up").c_str());

      }
      
      if(config.WZSelection=="wminus" || config.WZSelection=="combined"){
	
	if(config.PlotZmumu)	  h_zmumuN   = (TH1D*)grzmumuN    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZtautau)    h_ztautauN = (TH1D*)grztautauN  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWminTaunu)  h_wmintaunu = (TH1D*)grwmintaunu->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWminMunu )  h_wminmunu  = (TH1D*)grwminmunu ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWminMunu )  sumN        = (TH1D*)grwminmunu ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	h_DataN = (TH1D*)grDataN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotTtbar)          h_ttbarN = (TH1D*)grttbarN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());	
	if(config.PlotWt_top)         h_wt_topN = (TH1D*)grwt_topN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWt_antitop)     h_wt_antitopN = (TH1D*)grwt_antitopN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_top)     h_single_topN = (TH1D*)grsingle_topN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_antitop) h_single_antitopN = (TH1D*)grsingle_antitopN->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotZzqqll)         h_zzqqllN     = (TH1D*)grzzqqllN    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwqqll)         h_wzqqllN     = (TH1D*)grwzqqllN    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwpqqmlnu)      h_wwpqqmlnuN  = (TH1D*)grwwpqqmlnuN ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwplnumqq)      h_wwplnumqqN  = (TH1D*)grwwplnumqqN ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnuqq)        h_wzlnuqqN    = (TH1D*)grwzlnuqqN   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllll)         h_zzllllN     = (TH1D*)grzzllllN    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnull)        h_wzlnullN    = (TH1D*)grwzlnullN   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnununu)      h_wzlnununuN  = (TH1D*)grwzlnununuN ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllnunu)	      h_zzllnunuN   = (TH1D*)grzzllnunuN  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.SysInplots) h_sysDownN    = (TH1D*)grSysN       ->GetObjectUnchecked(("h_"+kine[iKine]+"_down").c_str());
	if(config.SysInplots) h_sysUpN      = (TH1D*)grSysN       ->GetObjectUnchecked(("h_"+kine[iKine]+"_up").c_str());
      }
 
      if(config.WZSelection=="zmumu"){
	
	if(config.PlotZmumu)	  h_zmumuZ      = (TH1D*)grzmumuZ    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZtautau)    h_ztautauZ    = (TH1D*)grztautauZ  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWplusTaunu) h_wplustaunu  = (TH1D*)grwplustaunu->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWplusMunu)  h_wplusmunu   = (TH1D*)grwplusmunu ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWminTaunu)  h_wmintaunu   = (TH1D*)grwmintaunu ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWminMunu )  h_wminmunu    = (TH1D*)grwminmunu  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZmumu)      sumZ = (TH1D*)h_zmumuZ->Clone("sumZ");
 
	h_DataZ = (TH1D*)grDataZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotTtbar)          h_ttbarZ           = (TH1D*)grttbarZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());	
	if(config.PlotWt_top)         h_wt_topZ          = (TH1D*)grwt_topZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWt_antitop)     h_wt_antitopZ      = (TH1D*)grwt_antitopZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_top)     h_single_topZ      = (TH1D*)grsingle_topZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotSingle_antitop) h_single_antitopZ  = (TH1D*)grsingle_antitopZ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());

	if(config.PlotZzqqll)         h_zzqqllZ     = (TH1D*)grzzqqllZ    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwqqll)         h_wzqqllZ     = (TH1D*)grwzqqllZ    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwpqqmlnu)      h_wwpqqmlnuZ  = (TH1D*)grwwpqqmlnuZ ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWwplnumqq)      h_wwplnumqqZ  = (TH1D*)grwwplnumqqZ ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnuqq)        h_wzlnuqqZ    = (TH1D*)grwzlnuqqZ   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllll)         h_zzllllZ     = (TH1D*)grzzllllZ    ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnull)        h_wzlnullZ    = (TH1D*)grwzlnullZ   ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotWzlnununu)      h_wzlnununuZ  = (TH1D*)grwzlnununuZ ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.PlotZzllnunu)	      h_zzllnunuZ   = (TH1D*)grzzllnunuZ  ->GetObjectUnchecked(("h_"+kine[iKine]+ran).c_str());
	if(config.SysInplots) h_sysDownZ    = (TH1D*)grSysZ       ->GetObjectUnchecked(("h_"+kine[iKine]+"_down").c_str());
	if(config.SysInplots) h_sysUpZ      = (TH1D*)grSysZ       ->GetObjectUnchecked(("h_"+kine[iKine]+"_up").c_str());

      }


  if(config.WZSelection=="wplus" || config.WZSelection=="combined"){

    newfile.cd();
    if(config.PlotZmumu)      h_zmumuP      ->Write("zmumuP");
    if(config.PlotZtautau)    h_ztautauP    ->Write("ztautauP");
    if(config.PlotWplusTaunu) h_wplustaunu  ->Write("wplustaunu");
    if(config.PlotWplusMunu)  h_wplusmunu   ->Write("wplusmunu");
    if(config.PlotWplusMunu)  sumP          ->Write("sumP");
    
    h_DataP ->Write("dataP");
    
    if(config.PlotTtbar)          h_ttbarP          ->Write("ttbarP");
    if(config.PlotWt_top)         h_wt_topP         ->Write("wt_topP");
    if(config.PlotWt_antitop)     h_wt_antitopP     ->Write("wt_antitopP");
    if(config.PlotSingle_top)     h_single_topP     ->Write("single_topP");
    if(config.PlotSingle_antitop) h_single_antitopP ->Write("single_antitopP");

    if(config.PlotZzqqll)         h_zzqqllP     ->Write("zzqqllP");
    if(config.PlotWwqqll)         h_wzqqllP     ->Write("wzqqllP");
    if(config.PlotWwpqqmlnu)      h_wwpqqmlnuP  ->Write("wwpqqmlnuP");
    if(config.PlotWwplnumqq)      h_wwplnumqqP  ->Write("wwplnumqqP");
    if(config.PlotWzlnuqq)        h_wzlnuqqP    ->Write("wzlnuqqP");
    if(config.PlotZzllll)         h_zzllllP     ->Write("zzllllP");
    if(config.PlotWzlnull)        h_wzlnullP    ->Write("wzlnullP");
    if(config.PlotWzlnununu)      h_wzlnununuP  ->Write("wzlnununuP");
    if(config.PlotZzllnunu)       h_zzllllP     ->Write("zzllllP");
    if(config.SysInplots) h_sysDownP    ->Write("sysdownP");
    if(config.SysInplots) h_sysUpP      ->Write("sysupP");
    
    newfile.Write();
  }

  if(config.WZSelection=="wminus" || config.WZSelection=="combined"){

    newfile.cd();
    if(h_zmumuN==NULL){
      std::cout<<kine[iKine]<<"     StackPlots::Message"<<std::endl;
      exit(10);
    }

    if(config.PlotZmumu)      h_zmumuN     ->Write("zmumuN");	
    if(config.PlotZtautau)    h_ztautauN   ->Write("ztautauN");	
    if(config.PlotWminTaunu)  h_wmintaunu  ->Write("wmintaunu");
    if(config.PlotWminMunu)   h_wminmunu   ->Write("wminmunu");	
    if(config.PlotWminMunu)   sumN         ->Write("sumN");      
    
    h_DataN ->Write("dataN");
   
    if(config.PlotTtbar)          h_ttbarN          ->Write("ttbarN");	       
    if(config.PlotWt_top)         h_wt_topN         ->Write("wt_topN");	       
    if(config.PlotWt_antitop)     h_wt_antitopN     ->Write("wt_antitopN");    
    if(config.PlotSingle_top)     h_single_topN     ->Write("single_topN");    
    if(config.PlotSingle_antitop) h_single_antitopN ->Write("single_antitopN");
   
    if(config.PlotZzqqll)         h_zzqqllN     ->Write("zzqqllN");   
    if(config.PlotWwqqll)         h_wzqqllN     ->Write("wzqqllN");   
    if(config.PlotWwpqqmlnu)      h_wwpqqmlnuN  ->Write("wwpqqmlnuN");
    if(config.PlotWwplnumqq)      h_wwplnumqqN  ->Write("wwplnumqqN");
    if(config.PlotWzlnuqq)        h_wzlnuqqN    ->Write("wzlnuqqN");  
    if(config.PlotZzllll)         h_zzllllN     ->Write("zzllllN");   
    if(config.PlotWzlnull)        h_wzlnullN    ->Write("wzlnullN");  
    if(config.PlotWzlnununu)      h_wzlnununuN  ->Write("wzlnununuN");
    if(config.PlotZzllnunu)	  h_zzllllN     ->Write("zzllllN");   
    if(config.SysInplots) h_sysDownN    ->Write("sysdownN");  
    if(config.SysInplots) h_sysUpN      ->Write("sysupN");      

    newfile.Write();      
  }


  if(config.WZSelection=="zmumu"){
    newfile.cd();

    if(config.PlotZmumu)      h_zmumuZ     ->Write("zmumu");	
    if(config.PlotZtautau)    h_ztautauZ   ->Write("ztautau");	 
    if(config.PlotWplusTaunu) h_wplustaunu ->Write("wplustaunu"); 
    if(config.PlotWplusMunu)  h_wplusmunu  ->Write("wplusmunu");	 
    if(config.PlotWminTaunu)  h_wmintaunu  ->Write("wmintaunu"); 
    if(config.PlotWminMunu )  h_wminmunu   ->Write("wminmunu"); 
    if(config.PlotZmumu)      sumZ         ->Write("sumZ");
    
    h_DataZ->Write("dataZ");
    
    if(config.PlotTtbar)          h_ttbarZ           ->Write("ttbarZ");	       
    if(config.PlotWt_top)         h_wt_topZ          ->Write("wt_topZ");	       
    if(config.PlotWt_antitop)     h_wt_antitopZ      ->Write("wt_antitopZ");    
    if(config.PlotSingle_top)     h_single_topZ      ->Write("single_topZ");    
    if(config.PlotSingle_antitop) h_single_antitopZ  ->Write("single_antitopZ");
    
    if(config.PlotZzqqll)         h_zzqqllZ     ->Write("zzqqllZ");   
    if(config.PlotWwqqll)         h_wzqqllZ     ->Write("wzqqllZ");   
    if(config.PlotWwpqqmlnu)      h_wwpqqmlnuZ  ->Write("wwpqqmlnuZ");
    if(config.PlotWwplnumqq)      h_wwplnumqqZ  ->Write("wwplnumqqZ");
    if(config.PlotWzlnuqq)        h_wzlnuqqZ    ->Write("wzlnuqqZ");  
    if(config.PlotZzllll)         h_zzllllZ     ->Write("zzllllZ");   
    if(config.PlotWzlnull)        h_wzlnullZ    ->Write("wzlnullZ");  
    if(config.PlotWzlnununu)      h_wzlnununuZ  ->Write("wzlnununuZ");
    if(config.PlotZzllnunu)	  h_zzllnunuZ   ->Write("zzllllZ");   
    if(config.SysInplots) h_sysDownZ    ->Write("sysdownZ");  
    if(config.SysInplots) h_sysUpZ      ->Write("sysupZ");    
   
    newfile.Write();      
  }



    // delete h_zmumuP;
    // delete h_ztautauP;
    // delete h_ttbarP;
    // delete h_wplustaunuP;
    // delete h_wplusmunuP;
    // delete sumP;
    // delete h_DataP;
    
    // delete h_wt_topP;
    // delete h_wt_antitopP;
    // delete h_single_topP;
    // delete h_single_antitopP;
    // delete h_zzqqllP;
    // delete h_wzqqllP;
    // delete h_wzlnuqqP;
    
    // delete h_wwlnuqqP;
    // delete h_wwpqqmlvP;
    
    // delete h_zzllnunuP;
    // delete h_zwlllnu_SFMinusP;
    // delete h_zwlllnu_OFMinusP;
    // delete h_zwlllnu_SFPlusP;
    // delete h_zwlllnu_OFPlusP;
    // delete h_zzllllP;

    // delete h_zmumuN;
    // delete h_ztautauN;
    // delete h_ttbarN;
    // delete h_wmintaunuN;
    // delete h_wminmunuN;
    // delete sumN;
    // delete h_DataN;
    
    // delete h_wt_topN;
    // delete h_wt_antitopN;
    // delete h_single_topN;
    // delete h_single_antitopN;
    // delete h_zzqqllN;
    // delete h_wzqqllN;
    // delete h_wzlnuqqN;
    
    // delete h_wwlnuqqN;
    // delete h_wwpqqmlvN;
    
    // delete h_zzllnunuN;
    // delete h_zwlllnu_SFMinusN;
    // delete h_zwlllnu_OFMinusN;
    // delete h_zwlllnu_SFPlusN;
    // delete h_zwlllnu_OFPlusN;
    // delete h_zzllllN;

    newfile.Close(); 
    }//kine loop
  }//range loop
 
  return; 
}

#endif


// double IntW=0., IntData=0., IntZMUMU=0., IntZTAUTAU=0., IntWTAU =0., IntTop =0., IntDiboson =0., Total=0.;

//   if(config.WZSelection=="wplus" || config.WZSelection=="combined"){
//     IntW=h_wplusmunuP->Integral(0,10000000);
//     IntData=h_DataP->Integral(0,1000000000);
//     IntZMUMU = h_zmumuP->Integral(0,10000000);
//     IntZTAUTAU = h_ztautauP->Integral(0,10000000);
//     IntWTAU = h_wplustaunuP->Integral(0,10000000);
//     IntTop = h_ttbarP->Integral(0,10000000)+ h_wt_antitopP->Integral(0,10000000)+h_single_topP->Integral(0,10000000)+h_single_antitopP->Integral(0,10000000)+h_wt_topP->Integral(0,10000000);
//     IntDiboson =  h_wzqqllP->Integral(0,10000000)+h_wzlnuqqP->Integral(0,10000000)+h_zzqqllP->Integral(0,10000000)+/*h_wwplvmqqP->Integral(0,10000000)*/+h_wwpqqmlvP->Integral(0,10000000)+h_zzllnunuP->Integral(0,10000000)+h_zwlllnu_SFMinusP->Integral(0,10000000)+h_zwlllnu_OFMinusP->Integral(0,10000000)+h_zwlllnu_SFPlusP->Integral(0,10000000)+h_zwlllnu_OFPlusP->Integral(0,10000000)+h_zzllllP->Integral(0,10000000);

//     Total = IntW + IntZMUMU + IntZTAUTAU + IntWTAU + IntTop + IntDiboson;

//     std::cout<<"************ Positive Channel ************"<<std::endl;
    
//     std::cout<<"Signal:     "<<"Events:  "<<IntW<<"     "<<IntW/Total<<std::endl;
//     std::cout<<"ZMUMU:      "<<"Events:  "<<IntZMUMU<<"     "<<IntZMUMU/Total<<std::endl;
//     std::cout<<"ZTAUTAU:    "<<"Events:  "<<IntZTAUTAU<<"     "<<IntZTAUTAU/Total<<std::endl;
//     std::cout<<"WTAU:       "<<"Events:  "<<IntWTAU<<"     "<<IntWTAU/Total<<std::endl;
//     std::cout<<"Top:        "<<"Events:  "<<IntTop<<"     "<<IntTop/Total<<std::endl;
//     std::cout<<"Diboson:    "<<"Events:  "<<IntDiboson<<"     "<<IntDiboson/Total<<std::endl;
 
//     std::cout<<"************ Event Yields ************"<<std::endl;
//     std::cout<<"Total MC    "<<Total<<"     Total Data    "<<IntData<<std::endl;
//   }

  
//   IntW=0.; IntData=0.; IntZMUMU=0.; IntZTAUTAU=0.; IntWTAU=0., IntTop=0.; IntDiboson=0.; Total=0.; 
 
//   if(config.WZSelection=="wminus" || config.WZSelection=="combined"){
//     IntW=h_wminmunuN->Integral(0,10000000);
//     IntData=h_DataN->Integral(0,1000000000);
//     IntZMUMU = h_zmumuN->Integral(0,10000000);
//     IntZTAUTAU = h_ztautauN->Integral(0,10000000);
//     IntWTAU = h_wmintaunuN->Integral(0,10000000);
//     IntTop = h_ttbarN->Integral(0,10000000)+ h_wt_antitopN->Integral(0,10000000)+h_single_topN->Integral(0,10000000)+h_single_antitopN->Integral(0,10000000)+h_wt_topN->Integral(0,10000000);
//     IntDiboson =  h_wzqqllN->Integral(0,10000000)+h_wzlnuqqN->Integral(0,10000000)+h_zzqqllN->Integral(0,10000000)/*+wwplvmqqN->Integral(0,10000000)+wwpqqmlvN->Integral(0,10000000)*/+h_zzllnunuN->Integral(0,10000000)+h_zwlllnu_SFMinusN->Integral(0,10000000)+h_zwlllnu_OFMinusN->Integral(0,10000000)+h_zwlllnu_SFPlusN->Integral(0,10000000)+h_zwlllnu_OFPlusN->Integral(0,10000000)+h_zzllllN->Integral(0,10000000);

//     Total = IntW + IntZMUMU + IntZTAUTAU + IntWTAU + IntTop + IntDiboson;
    
//     std::cout<<"************ Negative Channel ************"<<std::endl;
   
//     std::cout<<"Signal:     "<<"Events:  "<<IntW<<"     "<<IntW/Total<<std::endl;
//     std::cout<<"ZMUMU:      "<<"Events:  "<<IntZMUMU<<"     "<<IntZMUMU/Total<<std::endl;
//     std::cout<<"ZTAUTAU:    "<<"Events:  "<<IntZTAUTAU<<"     "<<IntZTAUTAU/Total<<std::endl;
//     std::cout<<"WTAU:       "<<"Events:  "<<IntWTAU<<"     "<<IntWTAU/Total<<std::endl;
//     std::cout<<"Top:        "<<"Events:  "<<IntTop<<"     "<<IntTop/Total<<std::endl;
//     std::cout<<"Diboson:    "<<"Events:  "<<IntDiboson<<"     "<<IntDiboson/Total<<std::endl;
   
//     std::cout<<"************ Event Yields ************"<<std::endl;
//     std::cout<<"Total MC    "<<Total<<"     Total Data    "<<IntData<<std::endl;
//   }
