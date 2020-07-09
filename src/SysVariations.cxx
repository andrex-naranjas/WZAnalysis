//Systematic Uncertainties class
#ifndef SYSVARIATIONS_CXX
#define SYSVARIATIONS_CXX

#include "SysVariations.h"

#include <cmath>
#include <string>
#include <iostream>

SysVariations::SysVariations()
{
}

SysVariations::~SysVariations(){}

void SysVariations::initialize(Config config){

  //names for input Files
  if(config.WZSelection=="zmumu"){
    wzchannel="z";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus";}

  if(config.SETCalibration=="True" && config.InsituCorrection=="True"){
    calib="_set_insitu";
  }else if(config.SETCalibration=="True" && config.InsituCorrection=="False"){
    calib="_set";
    }else if(config.SETCalibration=="False" && config.InsituCorrection=="True"){
    calib="_insitu";
  }else if(config.SETCalibration=="False" && config.InsituCorrection=="False"){
    calib="";
  }

  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";

  //ouput directory according year
  if(config.DataYears=="2015+2016") dirYear="2015p2016/";
  if(config.DataYears=="2017") dirYear="2017/";
  if(config.DataYears=="2018") dirYear="2018/";
  if(config.DataYears=="Full") dirYear="Full/";

  if(config.OnlyInclusive!="True"){ dirInclusive="Add/"; tot="_Total";}

  newfile = new TFile((config.OutputFileDir +"Sys/" +dirYear + dirInclusive + "sys_" + wzchannel + tot + ".root").c_str(),"RECREATE");

  kine.clear();
  kine.push_back("mu");   
  kine.push_back("pt");   
  kine.push_back("eta");  
  kine.push_back("eta_d");  
  kine.push_back("phi");  
  kine.push_back("phi_met");  
  kine.push_back("met");      
  kine.push_back("rmass");
  kine.push_back("delta_phi");
  kine.push_back("d0");
  
  // kine.push_back("mu_high"); 
  // kine.push_back("pt_high"); 
  // kine.push_back("eta_high");
  // kine.push_back("phi_high");
 
  // kine.push_back("phi_met_high"); 
  // kine.push_back("met_high");    
  // kine.push_back("rmass_high");  

  // kine.push_back("SumET");         kine.push_back("SumET_high");
  // kine.push_back("met_recoil");    kine.push_back("met_recoil_high");
  // kine.push_back("rmass_recoil");  kine.push_back("rmass_recoil_high");

  // kine.push_back("pt20Calo");   kine.push_back("pt20Calo_high");
  // kine.push_back("pt20Track");  kine.push_back("pt20Track");
  // kine.push_back("delta_phi");  kine.push_back("delta_phi");
  // kine.push_back("z0");         kine.push_back("z0_high");
  // kine.push_back("d0");         kine.push_back("d0_high");



  samples.clear();
  if(wzchannel=="wplus")  samples.push_back("wplusmunu");
  if(wzchannel=="wminus") samples.push_back("wminmunu");
  // samples.push_back("wplustaunu");
  // samples.push_back("wmintaunu");
  // samples.push_back("zmumu");
  // samples.push_back("ztautau");
  // samples.push_back("ttbar");
  // samples.push_back("wt_top");
  // samples.push_back("wt_antitop");
  // samples.push_back("single_top");
  // samples.push_back("single_antitop");
  // samples.push_back("zzqqll");
  // samples.push_back("wzqqll");
  // samples.push_back("wzlnuqq");
  // samples.push_back("wwplvmqq");
  // samples.push_back("wwpqqmlv");
  // samples.push_back("wwlnuqq");
  // samples.push_back("zzllnunu");
  // samples.push_back("zwlllnu_SFMinus");
  // samples.push_back("zwlllnu_OFMinus");
  // samples.push_back("zwlllnu_SFPlus");
  // samples.push_back("zwlllnu_OFPlus");
  // samples.push_back("zzllll");

  systematics_down.clear();//Systematics down;
  //systematics_down.push_back("_nominal");

  systematics_down.push_back("MUON_ID__1down");//0
  systematics_down.push_back("MUON_MS__1down");
  systematics_down.push_back("MUON_SAGITTA_RESBIAS__1down");
  systematics_down.push_back("MUON_SAGITTA_RHO__1down");
  systematics_down.push_back("MUON_SCALE__1down");

  systematics_down.push_back("MET_SoftTrk_ResoPara");//5
  systematics_down.push_back("MET_SoftTrk_ResoPerp");
  systematics_down.push_back("MET_SoftTrk_ScaleDown");
  systematics_down.push_back("PileUpDown");

  systematics_down.push_back("IDStatD");//9
  systematics_down.push_back("IDSysD");
  systematics_down.push_back("IsoStatD");
  systematics_down.push_back("IsoSysD");	      

  systematics_down.push_back("TrigStatD");//13
  systematics_down.push_back("TrigSysD");
  systematics_down.push_back("TTVAStatD");
  systematics_down.push_back("TTVASysD");

  systematics_down.push_back("kfactor_ALPHAS__1down");//17
  systematics_down.push_back("kfactor_BEAM_ENERGY__1down");
  systematics_down.push_back("kfactor_PDF_EW__1down");
  systematics_down.push_back("kfactor_PDF__1down");
  systematics_down.push_back("kfactor_PI__1down");

  systematics_down.push_back("kfactor_SCALE_W__1down");//22
  systematics_down.push_back("kfactor_SCALE_Z__1down");
  systematics_down.push_back("kfactor_CHOICE_HERAPDF20");
  systematics_down.push_back("kfactor_CHOICE_NNPDF30");

  systematics_down.push_back("kfactor_PDF_EV1");//26
  systematics_down.push_back("kfactor_PDF_EV2");
  systematics_down.push_back("kfactor_PDF_EV3");
  systematics_down.push_back("kfactor_PDF_EV4");
  systematics_down.push_back("kfactor_PDF_EV5");
  systematics_down.push_back("kfactor_PDF_EV6");
  systematics_down.push_back("kfactor_PDF_EV7");

  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_EtaIntercalibration_NonClosure__1down");//33
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_DataVsMC__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_1__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_2__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_3__1down");

  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_1__1down");//38
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_2__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_3__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_4__1down");

  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_5__1down");//42
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_6__1down");
  systematics_down.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_7restTerm__1down");

  systematics_up.clear(); //Systematics UP
  // systematics_up.push_back("_nominal");

  systematics_up.push_back("MUON_ID__1up");
  systematics_up.push_back("MUON_MS__1up");
  systematics_up.push_back("MUON_SAGITTA_RESBIAS__1up");
  systematics_up.push_back("MUON_SAGITTA_RHO__1up");
  systematics_up.push_back("MUON_SCALE__1up");

  systematics_up.push_back("MET_SoftTrk_ResoPara");//5
  systematics_up.push_back("MET_SoftTrk_ResoPerp");
  systematics_up.push_back("MET_SoftTrk_ScaleUp");
  systematics_up.push_back("PileUpUp");

  systematics_up.push_back("IDStatU");//9
  systematics_up.push_back("IDSysU");  
  systematics_up.push_back("IsoStatU");
  systematics_up.push_back("IsoSysU"); 

  systematics_up.push_back("TrigStatU");//13
  systematics_up.push_back("TrigSysU");
  systematics_up.push_back("TTVAStatU");
  systematics_up.push_back("TTVASysU");

  systematics_up.push_back("kfactor_ALPHAS__1up");//17
  systematics_up.push_back("kfactor_BEAM_ENERGY__1up");
  systematics_up.push_back("kfactor_PDF_EW__1up");
  systematics_up.push_back("kfactor_PDF__1up");
  systematics_up.push_back("kfactor_PI__1up");

  systematics_up.push_back("kfactor_SCALE_W__1up");//22
  systematics_up.push_back("kfactor_SCALE_Z__1up");
  systematics_up.push_back("kfactor_CHOICE_HERAPDF20");
  systematics_up.push_back("kfactor_CHOICE_NNPDF30");

  systematics_up.push_back("kfactor_PDF_EV1");//26
  systematics_up.push_back("kfactor_PDF_EV2");
  systematics_up.push_back("kfactor_PDF_EV3");
  systematics_up.push_back("kfactor_PDF_EV4");
  systematics_up.push_back("kfactor_PDF_EV5");
  systematics_up.push_back("kfactor_PDF_EV6");
  systematics_up.push_back("kfactor_PDF_EV7");

  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_EtaIntercalibration_NonClosure__1up");//33
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_DataVsMC__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_1__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_2__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_3__1up");

  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_1__1up");//38
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_2__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_3__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_4__1up");

  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_5__1up");//42
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_6__1up");
  systematics_up.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_7restTerm__1up");

  return;
}


void SysVariations::execute(Config config){
  
  for(int iKine = 0; iKine<(int)kine.size();iKine++){

    std::cout<<"*** Extracting Uncertainties for variable:   "<<kine[iKine]<<std::endl;    
    double total=0.; double total_noSF=0.;

    TH1D *h_nom[100]={};     
    TH1D *h_nom_noSF[100]={};
    TH1D *h_down_a[100][100]={};   TH1D *h_up_a[100][100]={};

    for (int iSamples = 0; iSamples < (int)samples.size(); iSamples++){

      TFile *GR_nom = new TFile((config.OutputFileDir +"Files/" + dirYear + dirInclusive + samples[iSamples] + "_"+wzchannel + "_nominal" + calib + puname + tot +".root").c_str());      
      TFile *GR_nom_noSF = new TFile((config.OutputFileDir +"Files/" + dirYear + dirInclusive + samples[iSamples] + "_"+wzchannel + "_nominal_noSF" + calib + puname + tot +".root").c_str());

      h_nom[iSamples] = (TH1D*)GR_nom->Get(("h_"+kine[iKine]).c_str());
      total=h_nom[iSamples]->Integral(0,h_nom[iSamples]->GetNbinsX()+1);

      h_nom_noSF[iSamples] = (TH1D*)GR_nom_noSF->Get(("h_"+kine[iKine]).c_str());
      total_noSF=h_nom_noSF[iSamples]->Integral(0,h_nom_noSF[iSamples]->GetNbinsX()+1);
    
      //binning for the final histos
      unsigned int nx = h_nom[iSamples]->GetXaxis()->GetNbins();
      Double_t *xbins = new Double_t[h_nom[iSamples]->GetXaxis()->GetNbins()+1];
      for(unsigned int i=0; i< h_nom[iSamples]->GetXaxis()->GetNbins()+1; i++) xbins[i]=h_nom[iSamples]->GetXaxis()->GetBinLowEdge(i+1);
      xbins[h_nom[iSamples]->GetXaxis()->GetNbins()]=h_nom[iSamples]->GetXaxis()->GetBinUpEdge(h_nom[iSamples]->GetXaxis()->GetNbins());
      
      h_down     = new TH1D(("h_"+kine[iKine]+"_down").c_str(),"",nx,xbins);
      h_up       = new TH1D(("h_"+kine[iKine]+"_up").c_str(),"",nx,xbins);
      h_stat     = new TH1D(("h_"+kine[iKine]+"_stat").c_str(),"",nx,xbins);
      
      h_statTEMP = (TH1D*)h_nom[iSamples]->Clone(("h_statTEMP_"+kine[iKine]).c_str()); h_statTEMP->Add(h_statTEMP,-1.);
      h_downTEMP = (TH1D*)h_nom[iSamples]->Clone(("h_downTEMP_"+kine[iKine]).c_str()); h_downTEMP->Add(h_downTEMP,-1.);
      h_upTEMP   = (TH1D*)h_nom[iSamples]->Clone(("h_upTEMP_"  +kine[iKine]).c_str()); h_upTEMP  ->Add(h_upTEMP,-1.);
      h_nomTEMP  = (TH1D*)h_nom[iSamples]->Clone(("h_nomTEMP_" +kine[iKine]).c_str()); h_nomTEMP ->Add(h_nomTEMP,-1.);
  
      //Getting the statistical errors. Alex Law macros dont get these numbers, so we have to get them ourselves.
      //bin = 0 underflow bin, bin = 1 first bin, bin = nbins last bin, bin=nbins+1 overflow
      Int_t nbins = h_nom[iSamples]->GetNbinsX();
      for(Int_t k = 0; k<nbins+2; k++) h_stat->SetBinContent(k,h_nom[iSamples]->GetBinError(k));

      h_statTEMP->Add(h_stat);
      h_nomTEMP->Add(h_nom[iSamples]);

      //histos to store the total systematic
      h_downInd  = new TH1D("h_downInd_Total", "", (int)systematics_down.size(), 0, (int)systematics_down.size());
      h_upInd    = new TH1D("h_upInd_Total",   "", (int)systematics_up.size(),   0, (int)systematics_up.size());

      if(kine[iKine]=="eta"){
	h_downInd_eta  = new TH2D("h_downInd_eta", "", nx, xbins, (int)systematics_down.size(), 0, (int)systematics_down.size());
	h_upInd_eta    = new TH2D("h_upInd_eta",   "", nx, xbins, (int)systematics_up.size(),   0, (int)systematics_up.size());
      }
     
      double sumDown=0., sumUp=0.;
      std::vector<std::vector<double>> value_up_tot, value_down_tot; value_up_tot.clear(); value_down_tot.clear();

      for(int iSysDown = 0; iSysDown < (int)systematics_down.size(); ++iSysDown){
	
	TFile *GR_down = new TFile((config.OutputFileDir +"Files/" + dirYear + dirInclusive + samples[iSamples] +"_"+ wzchannel +"_"+ systematics_down[iSysDown]+calib+puname+tot + ".root").c_str());

	h_down_a[iSamples][iSysDown] = (TH1D*)GR_down->Get(("h_"+kine[iKine]).c_str());

	//Individual systematics
	h_downInd->GetXaxis()->SetBinLabel(iSysDown+1, systematics_down[iSysDown].c_str());
	double value_down = 0;
	if(isScalekFactor(systematics_down[iSysDown]) || isScaleEffFactor(systematics_down[iSysDown]))
	  value_down = (h_down_a[iSamples][iSysDown]->Integral(0,h_down_a[iSamples][iSysDown]->GetNbinsX()+1)-total) / total;
	else
	  value_down = (h_down_a[iSamples][iSysDown]->Integral(0,h_down_a[iSamples][iSysDown]->GetNbinsX()+1)-total_noSF) / total_noSF;

	if(std::isnan(value_down) || std::isinf(value_down) || isScalekFactor(systematics_down[iSysDown])) value_down=0;//no kfactor sys(change)!
	h_downInd->SetBinContent(iSysDown+1,value_down);

	//std::cout<<systematics_down[iSysDown]<<"  test: "<<h_down_a[iSamples][iSysDown]->Integral(0,h_down_a[iSamples][iSysDown]->GetNbinsX()+1) - total<<std::endl;
      
	//Individual eta bin systematics
	if(kine[iKine]=="eta"){
	  //std::cout<<"systematic   "<<systematics_down[iSysDown]<<std::endl;
	  for(int k=1;k<nx+1;k++){

	    double value_eta_down=0.;
	    if(isScalekFactor(systematics_down[iSysDown]) || isScaleEffFactor(systematics_down[iSysDown]))
	      value_eta_down = (h_down_a[iSamples][iSysDown]->GetBinContent(k)-h_nom[iSamples]->GetBinContent(k))/h_nom[iSamples]->GetBinContent(k);
	    else
	      value_eta_down = (h_down_a[iSamples][iSysDown]->GetBinContent(k)-h_nom_noSF[iSamples]->GetBinContent(k))/h_nom_noSF[iSamples]->GetBinContent(k);
	      
	    if(std::isnan(value_eta_down) || std::isinf(value_eta_down) || isScalekFactor(systematics_down[iSysDown])) value_eta_down = 0.;//no kfactor sys(change)!
	    h_downInd_eta->GetYaxis()->SetBinLabel(iSysDown+1, systematics_down[iSysDown].c_str());
	    h_downInd_eta->SetBinContent(k,iSysDown+1,value_eta_down);
	    //std::cout<<h_down_a[iSamples][iSysDown]->GetBinContent(k)<<"    "<<h_nom[iSamples]->GetBinContent(k)<<"    "<<100*value_eta_down<<std::endl;
	  }
	}
	
	sumDown += std::pow((h_down_a[iSamples][iSysDown]->Integral(0,h_down_a[iSamples][iSysDown]->GetNbinsX()+1)-total)/total,2);//total uncertainty

	std::vector<double> value_down_vec; value_down_vec.clear();
	for(int k=1;k<nx+1;k++){
	  double dummyDown=0;
	  if(isScalekFactor(systematics_down[iSysDown]) || isScaleEffFactor(systematics_down[iSysDown]))
	    dummyDown = h_down_a[iSamples][iSysDown]->GetBinContent(k) - h_nom[iSamples]->GetBinContent(k);
	  else
	    dummyDown = h_down_a[iSamples][iSysDown]->GetBinContent(k) - h_nom_noSF[iSamples]->GetBinContent(k);
	
	  if(std::isnan(dummyDown)||std::isinf(dummyDown) || isScalekFactor(systematics_down[iSysDown])) dummyDown = 0;//no kfactor sys(change)!
	  value_down_vec.push_back(std::pow(dummyDown,2));
	}
	value_down_tot.push_back(value_down_vec);

	// if(isScalekFactor(systematics_down[iSysDown]) || isScaleEffFactor(systematics_down[iSysDown]))
	//   h_down_a[iSamples][iSysDown]->Add(h_nom[iSamples],-1.);
	// else
	//   h_down_a[iSamples][iSysDown]->Add(h_nom_noSF[iSamples],-1.);
	// h_down_a[iSamples][iSysDown]->Multiply(h_down_a[iSamples][iSysDown]);

	// if(isScalekFactor(systematics_down[iSysDown]))//remove these lines to include the kfactor sys
	//   h_down_a[iSamples][iSysDown]->Scale(0.);
		  
	//   if(!(std::isnan(h_down_a[iSamples][iSysDown]->Integral()) || std::isinf(h_down_a[iSamples][iSysDown]->Integral())))
	//     h_downTEMP->Add(h_down_a[iSamples][iSysDown]);

	delete GR_down;
      }//loop down systematics

      for(int iSysUp = 0; iSysUp < (int)systematics_up.size(); ++iSysUp){

	TFile *GR_up = new TFile((config.OutputFileDir +"Files/" + dirYear + dirInclusive + samples[iSamples] +"_"+ wzchannel +"_"+ systematics_up[iSysUp] + calib + puname + tot +".root").c_str());
	
	h_up_a[iSamples][iSysUp] = (TH1D*)GR_up->Get(("h_"+kine[iKine]).c_str());

	//Individual systematics
	h_upInd->GetXaxis()->SetBinLabel(iSysUp+1, systematics_up[iSysUp].c_str());

	double value_up=0.;
	if(isScalekFactor(systematics_up[iSysUp]) || isScaleEffFactor(systematics_up[iSysUp]))
	  value_up = (h_up_a[iSamples][iSysUp]->Integral(0,h_up_a[iSamples][iSysUp]->GetNbinsX()+1)-total)/total;
	else
	  value_up = (h_up_a[iSamples][iSysUp]->Integral(0,h_up_a[iSamples][iSysUp]->GetNbinsX()+1)-total_noSF)/total_noSF;
	  
	if(std::isnan(value_up) || std::isinf(value_up) || isScalekFactor(systematics_up[iSysUp])) value_up =0;//no kfactor sys(change)!
	h_upInd->SetBinContent(iSysUp+1,value_up);

	//Individual eta bin systematics
	if(kine[iKine]=="eta"){//underflow and overflow bins are empty!
	  for(int k=1;k<nx+1;k++){
	    h_upInd_eta->GetYaxis()->SetBinLabel(iSysUp+1, systematics_up[iSysUp].c_str());

	    double value_eta_up=0.;
	    if(isScalekFactor(systematics_up[iSysUp]) || isScaleEffFactor(systematics_up[iSysUp]))
	      value_eta_up = (h_up_a[iSamples][iSysUp]->GetBinContent(k)-h_nom[iSamples]->GetBinContent(k))/h_nom[iSamples]->GetBinContent(k);
	    else
	      value_eta_up = (h_up_a[iSamples][iSysUp]->GetBinContent(k)-h_nom_noSF[iSamples]->GetBinContent(k))/h_nom_noSF[iSamples]->GetBinContent(k);

	    if(std::isnan(value_eta_up)||std::isinf(value_eta_up) || isScalekFactor(systematics_up[iSysUp])) value_eta_up = 0;//no kfactor sys(change)!
	    h_upInd_eta->SetBinContent(k,iSysUp+1,value_eta_up);
	  }
	}

	sumUp += std::pow((h_up_a[iSamples][iSysUp]->Integral(0,h_up_a[iSamples][iSysUp]->GetNbinsX()+1)-total)/total,2);//total uncertainty

	std::vector<double> value_up_vec; value_up_vec.clear();
	for(int k=1;k<nx+1;k++){
	  double dummyUp=0;
	  if(isScalekFactor(systematics_up[iSysUp]) || isScaleEffFactor(systematics_up[iSysUp]))
	    dummyUp = h_up_a[iSamples][iSysUp]->GetBinContent(k) - h_nom[iSamples]->GetBinContent(k);
	  else
	    dummyUp = h_up_a[iSamples][iSysUp]->GetBinContent(k) - h_nom_noSF[iSamples]->GetBinContent(k);
	
	  if(std::isnan(dummyUp)||std::isinf(dummyUp) || isScalekFactor(systematics_up[iSysUp])) dummyUp = 0;//no kfactor sys(change)!
	  value_up_vec.push_back(std::pow(dummyUp,2));
	}
	value_up_tot.push_back(value_up_vec);	      	

	// if(isScalekFactor(systematics_up[iSysUp]) || isScaleEffFactor(systematics_up[iSysUp]))
	//   h_up_a[iSamples][iSysUp]->Add(h_nom[iSamples],-1.);
	// else
	//   h_up_a[iSamples][iSysUp]->Add(h_nom_noSF[iSamples],-1.);	  

	// h_up_a[iSamples][iSysUp]->Multiply(h_up_a[iSamples][iSysUp]);

	// if(isScalekFactor(systematics_up[iSysUp]))//remove these lines to include the kfactor sys
	//   h_up_a[iSamples][iSysUp]->Scale(0);
	
	// if(!(std::isnan(h_up_a[iSamples][iSysUp]->Integral()) || std::isinf(h_up_a[iSamples][iSysUp]->Integral())))
	//   h_upTEMP->Add(h_up_a[iSamples][iSysUp]);
	
	delete GR_up;
      }//loop up systematics

      std::cout<<"miami horror"<<std::endl;
      std::cout<<value_up_tot.size()<<std::endl; //number of systematics
      std::cout<<value_up_tot.at(1).size()<<std::endl; //number of bins

      for(int i=0; i<value_down_tot.at(0).size(); i++){
	double binvalueDown=0., binvalueUp=0.;
	for(int j=0; j<value_down_tot.size(); j++){
	  binvalueDown += (value_down_tot.at(j)).at(i);	
	  binvalueUp += (value_up_tot.at(j)).at(i);	
	}
      	h_upTEMP->SetBinContent(i+1,binvalueDown);
      	h_downTEMP->SetBinContent(i+1,binvalueUp);
      }
	
      //delete GR_nom;
      h_nom[iSamples]=NULL;
      
    }//loop samples


    unsigned int gbin = h_downTEMP->GetNbinsX();
    for(int i=1;i<gbin+1;i++){

      double downError = h_downTEMP->GetBinContent(i);//+h_statTEMP->GetBinContent(i)*h_statTEMP->GetBinContent(i);
      //double downError = h_statTEMP->GetBinContent(i)*h_statTEMP->GetBinContent(i);
      double downDelta =std::sqrt(downError);
      
      double upError = h_upTEMP->GetBinContent(i);//+h_statTEMP->GetBinContent(i)*h_statTEMP->GetBinContent(i);
      //double upError = h_statTEMP->GetBinContent(i)*h_statTEMP->GetBinContent(i);
      double upDelta=std::sqrt(upError);
      
      double statError = h_statTEMP->GetBinContent(i);

      double nomValue  = h_nomTEMP->GetBinContent(i);

      //std::cout<<(upDelta-statError)/statError<<"     "<<(downDelta-statError)/statError<<"    "<<(statError-statError)/statError<<std::endl;
      
      // Symmetric Sys + Stat
      Double_t deltaAv=(downDelta+upDelta)*(0.5);
      //std::cout<<deltaAv<<"     "<<nomValue<<"      "<<deltaAv/nomValue<<std::endl;

      h_down->SetBinContent(i,(deltaAv)*(-1.));
      h_up->SetBinContent(i,deltaAv);
    }

    h_downTEMP=NULL; h_downTEMP=NULL; h_upTEMP=NULL; h_statTEMP=NULL;
   
    newfile->cd();
    h_down->Write(); h_up->Write();
    
  }//kine loop

  //Plots
  setstyle();
  SysPlots(config,h_upInd_eta,h_downInd_eta);

  newfile->cd();
  h_upInd->Write(); h_downInd->Write();
  h_upInd_eta->Write(); h_downInd_eta->Write();
  newfile->Close();

  return;
}

void SysVariations::SysPlots(Config config,TH2D* h_up_eta,TH2D* h_dn_eta){

  TString boson;
  if(config.WZSelection=="wplus")  boson="W^{+}#rightarrow#mu^{+}#nu";
  if(config.WZSelection=="wminus") boson="W^{-}#rightarrow#mu^{-}#nu";

  std::vector<TH1D*> sysHistosUp, sysHistosDn; sysHistosUp.clear(); sysHistosDn.clear();
  std::vector<std::string> nameSys, nameSysUp, nameSysDn; nameSys.clear(); nameSysUp.clear(); nameSysDn.clear();

  unsigned int nx = h_up_eta->GetXaxis()->GetNbins();
  Double_t *xbins = new Double_t[h_up_eta->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< h_up_eta->GetXaxis()->GetNbins()+1; i++) xbins[i]=h_up_eta->GetXaxis()->GetBinLowEdge(i+1);
  xbins[h_up_eta->GetXaxis()->GetNbins()]=h_up_eta->GetXaxis()->GetBinUpEdge(h_up_eta->GetXaxis()->GetNbins());
      
  for(int j=1; j<h_up_eta->GetYaxis()->GetNbins()+1; j++){

    std::string dummy =std::to_string(j);
    TH1D *hAuxUp  = new TH1D(("hAuxUp"+dummy).c_str(), "", nx, xbins);
    TH1D *hAuxDn  = new TH1D(("hAuxDn"+dummy).c_str(), "", nx, xbins);

    for(int i=1; i<h_up_eta->GetXaxis()->GetNbins()+1; i++){
      //do symmetric errors
      double up_sys = h_up_eta->GetBinContent(i,j);
      double dn_sys = h_dn_eta->GetBinContent(i,j);
      double sy_sys =(std::fabs(up_sys)+std::fabs(dn_sys))*0.5;

      if(std::isinf(sy_sys) || std::isnan(sy_sys))//change this
	sy_sys=0;

      float sign= std::fabs(up_sys)/up_sys;
      sign=1;//check this

      hAuxUp->SetBinContent(i,sign*sy_sys*100);//percentage values (%)
      hAuxUp->SetBinError(i,0);
      hAuxDn->SetBinContent(i,-1.*sign*sy_sys*100);//percentage values (%)
      hAuxDn->SetBinError(i,0);      
    }
    nameSysUp.push_back(h_up_eta->GetYaxis()->GetBinLabel(j));
    nameSysDn.push_back(h_dn_eta->GetYaxis()->GetBinLabel(j));
    sysHistosUp.push_back(hAuxUp); sysHistosDn.push_back(hAuxDn); 
    //delete hAuxUp; delete hAuxDn;
  }

  std::vector<int> markers, colors;//10 groups of periodic color in plots
  for(int k=1;k<10+1;k++){
    colors.push_back(kBlue);     markers.push_back(20);     //colors.push_back(kCyan);   
    colors.push_back(kRed);   markers.push_back(21);        //colors.push_back(kCyan-6);
    colors.push_back(kGreen);   markers.push_back(22);      //colors.push_back(kBlue-3);
    if(k==1 || k==5 || k==8){                               //colors.push_back(kViolet-6);
      colors.push_back(kViolet-6); markers.push_back(23);   //colors.push_back(kViolet);  
      colors.push_back(kOrange);   markers.push_back(29); 
      continue;
    }else if(k==2 || k==3 || k==4 || k==6 || k==9){
      colors.push_back(kViolet-6); markers.push_back(23);
      continue;
    }else if(k==7){
	colors.push_back(kViolet-6); markers.push_back(23); 
	colors.push_back(kOrange);   markers.push_back(29); 
	colors.push_back(kBlack);      markers.push_back(33);
	colors.push_back(kGray);    markers.push_back(38);
	continue;
    }
  }

  std::vector<std::string> up_down; up_down.clear();
  up_down.push_back("sysUp"); up_down.push_back("sysDn");

  for(int UpDown=0;UpDown<2;UpDown++){

    std::vector<TH1D*> histVec; histVec.clear();
    if(up_down[UpDown]=="sysUp"){histVec=sysHistosUp; nameSys=nameSysUp;}
    if(up_down[UpDown]=="sysDn"){histVec=sysHistosDn; nameSys=nameSysDn;}

    std::vector<std::string> sysGroup; sysGroup.clear();
    sysGroup.push_back("MUON");
    sysGroup.push_back("MET");
    sysGroup.push_back("SFac_1");
    sysGroup.push_back("SFac_2");
    sysGroup.push_back("kFac_1");
    sysGroup.push_back("kFac_2");
    sysGroup.push_back("kFac_3");
    sysGroup.push_back("JET_1");
    sysGroup.push_back("JET_2");
    sysGroup.push_back("JET_3");

    
    for(int k=0;k<(int)sysGroup.size();k++){
      
      TString lumi="Simulation";
      TCanvas *c = new TCanvas("canvas","canvas",550,450);
      TLegend* leg = new TLegend(0.175,0.85,0.375,0.6,"");
      leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
      TString naming="";//config.DataYears;
      
      TPaveText *box;
      box = new TPaveText(0.2,0.865,0.375,0.965,"NDC");
      box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
      //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
      box->AddText("CCDY Analysis");
      box->AddText("13 TeV, "+lumi);

      TPaveText *box2;
      box2 = new TPaveText(0.7,0.845,0.8,0.965,"NDC");
      box2->SetBorderSize(0); box2->SetTextSize(0.05); box2->SetFillColor(0);
      box2->AddText(boson);
      c->cd();
      
      //main panel
      TPad *pad = new TPad("pad","pad",0,0,1,1);
      pad->SetTicks(1,1);
      pad->Draw();
      pad->cd();
      pad->SetLeftMargin(0.1125);
      pad->SetRightMargin(0.03);
      pad->SetTopMargin(0.015);
      pad->SetBottomMargin(0.1125);
      
      int nSys=0, InitnSys=0; double ymin=0., ymax=0.;
      std::cout<<"   HistVec size  "<<histVec.size()<<"   color size"<<colors.size()<<std::endl;
      
      if(sysGroup[k]=="MUON"    && up_down[UpDown]=="sysUp"){nSys=5;  InitnSys=0; ymin=-0.1; ymax=0.5;}
      if(sysGroup[k]=="MUON"    && up_down[UpDown]=="sysDn"){nSys=5;  InitnSys=0; ymin=-0.5; ymax=0.5;}
      if(sysGroup[k]=="MET"     && up_down[UpDown]=="sysUp"){nSys=9;  InitnSys=5; ymin=-0.1; ymax=1.5;}
      if(sysGroup[k]=="MET"     && up_down[UpDown]=="sysDn"){nSys=9;  InitnSys=5; ymin=-1; ymax=1;}
      if(sysGroup[k]=="SFac_1"  && up_down[UpDown]=="sysUp"){nSys=13; InitnSys=9; ymin=-0.25; ymax=1.75;}
      if(sysGroup[k]=="SFac_1"  && up_down[UpDown]=="sysDn"){nSys=13; InitnSys=9; ymin=-1.25; ymax=1.0;}
      if(sysGroup[k]=="SFac_2"  && up_down[UpDown]=="sysUp"){nSys=17; InitnSys=13; ymin=-0.25; ymax=1.5;}
      if(sysGroup[k]=="SFac_2"  && up_down[UpDown]=="sysDn"){nSys=17; InitnSys=13; ymin=-0.75; ymax=1.0;}
      if(sysGroup[k]=="kFac_1"  && up_down[UpDown]=="sysUp"){nSys=22; InitnSys=17; ymin=-1; ymax=10;}
      if(sysGroup[k]=="kFac_1"  && up_down[UpDown]=="sysDn"){nSys=22; InitnSys=17; ymin=-5; ymax=6;}
      if(sysGroup[k]=="kFac_2"  && up_down[UpDown]=="sysUp"){nSys=26; InitnSys=22; ymin=-1; ymax=10;}
      if(sysGroup[k]=="kFac_2"  && up_down[UpDown]=="sysDn"){nSys=26; InitnSys=22; ymin=-5; ymax=6;}
      if(sysGroup[k]=="kFac_3"  && up_down[UpDown]=="sysUp"){nSys=33; InitnSys=26; ymin=-1; ymax=10;}
      if(sysGroup[k]=="kFac_3"  && up_down[UpDown]=="sysDn"){nSys=33; InitnSys=26; ymin=-5; ymax=6;}
      if(sysGroup[k]=="JET_1"   && up_down[UpDown]=="sysUp"){nSys=38; InitnSys=33; ymin=-0.1; ymax=2.75;}
      if(sysGroup[k]=="JET_1"   && up_down[UpDown]=="sysDn"){nSys=38; InitnSys=33; ymin=-1.75; ymax=1.75;}
      if(sysGroup[k]=="JET_2"   && up_down[UpDown]=="sysUp"){nSys=42; InitnSys=38; ymin=-0.1; ymax=1.25;}
      if(sysGroup[k]=="JET_2"   && up_down[UpDown]=="sysDn"){nSys=42; InitnSys=38; ymin=-0.6; ymax=0.75;}
      if(sysGroup[k]=="JET_3"   && up_down[UpDown]=="sysUp"){nSys=45; InitnSys=42; ymin=-0.1; ymax=0.75;}
      if(sysGroup[k]=="JET_3"   && up_down[UpDown]=="sysDn"){nSys=45; InitnSys=42; ymin=-0.4; ymax=0.4;}
      
      for(int iPlot=InitnSys;iPlot<nSys;iPlot++){
	histVec[iPlot]->SetLineColor(colors[iPlot]);
	histVec[iPlot]->SetMarkerStyle(markers[iPlot]);
	histVec[iPlot]->SetMarkerColor(colors[iPlot]);
	histVec[iPlot]->SetMarkerSize(0.75);
	histVec[iPlot]->SetTitle("");
	histVec[iPlot]->SetMinimum(ymin);
	histVec[iPlot]->SetMaximum(ymax);
	histVec[iPlot]->GetXaxis()->SetTitle("#||{#eta^{#mu}}");
	histVec[iPlot]->GetYaxis()->SetTitle("#delta^{sys} (%)");
	histVec[iPlot]->SetLineWidth(2.5);
	
	TString name_sys = nameSys[iPlot].c_str();
	leg->AddEntry(histVec[iPlot],"#scale[0.85]{#bf{"+name_sys+"}}","L");
	
	if(iPlot==InitnSys){
	  histVec[iPlot]->DrawCopy("hist");
	  //histVec[iPlot]->DrawCopy("hist p sames");
	}else{
	  histVec[iPlot]->DrawCopy("hist sames");
	  //histVec[iPlot]->DrawCopy("hist p sames");
	}	
      }
      
      pad->Update();
      pad->Modified();
      
      // TLine *line;
      // line = new TLine(-2.5,1.,2.5,1.);
      // line->SetLineColor(kBlack);
      // line->SetLineStyle(2);
      //line->Draw();
      leg->Draw();
      box->Draw();
      box2->Draw();
      c->cd();
      //c->Draw();
      std::cout<<("./systematics_"+sysGroup[k]+"_"+up_down[UpDown]+"_"+config.WZSelection+".pdf").c_str()<<std::endl;
      c->Print(("./systematics_"+sysGroup[k]+"_"+up_down[UpDown]+"_"+config.WZSelection+".pdf").c_str());
      //histVec.clear();

      delete c;      
      delete leg;
      // delete line;
      delete box;
    }
  }

  return;
}

bool SysVariations::isScalekFactor(std::string sys){

  bool isKfactor=false;

  if (sys=="kfactor_ALPHAS__1up" || sys=="kfactor_BEAM_ENERGY__1up" || sys=="kfactor_PDF_EW__1up" || sys=="kfactor_PDF__1up" || sys=="kfactor_PI__1up"  || sys=="kfactor_SCALE_W__1up" || sys=="kfactor_SCALE_Z__1up"  || sys=="kfactor_CHOICE_HERAPDF20"  || sys=="kfactor_CHOICE_NNPDF30"  || sys=="kfactor_PDF_EV1"  || sys=="kfactor_PDF_EV2"  || sys=="kfactor_PDF_EV3"  || sys=="kfactor_PDF_EV4"  || sys=="kfactor_PDF_EV5"  || sys=="kfactor_PDF_EV6" || sys=="kfactor_PDF_EV7")
    isKfactor = true;

  if(sys=="kfactor_ALPHAS__1down" || sys=="kfactor_BEAM_ENERGY__1down" || sys=="kfactor_PDF_EW__1down" || sys=="kfactor_PDF__1down"	   || sys=="kfactor_PI__1down"  || sys=="kfactor_SCALE_W__1down" || sys=="kfactor_SCALE_Z__1down"  || sys=="kfactor_CHOICE_HERAPDF20"  || sys=="kfactor_CHOICE_NNPDF30"  || sys=="kfactor_PDF_EV1"  || sys=="kfactor_PDF_EV2"  || sys=="kfactor_PDF_EV3"  || sys=="kfactor_PDF_EV4"  || sys=="kfactor_PDF_EV5"  || sys=="kfactor_PDF_EV6" || sys=="kfactor_PDF_EV7")
    isKfactor=true;

  return isKfactor;
}

bool SysVariations::isScaleEffFactor(std::string sys){

  bool isfactor=false;

  if(sys=="PileUpDown" ||  sys=="IDStatD" ||  sys=="IDSysD" ||  sys=="IsoStatD" ||  sys=="IsoSysD" ||  sys=="TrigStatD" || sys=="TrigSysD" ||  sys=="TTVAStatD" ||  sys=="TTVASysD" ||  sys=="PileUpUp" ||  sys=="IDStatU" ||  sys=="IDSysU" ||    sys=="IsoStatU" ||  sys=="IsoSysU" ||   sys=="TrigStatU" ||  sys=="TrigSysU" ||  sys=="TTVAStatU" || sys=="TTVASysU") 
    isfactor=true;

  return isfactor;
}

void SysVariations::setstyle(){

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

  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Detector1__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical1__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical4__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed2__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical3__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling1__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical6__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling3__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling4__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed3__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling2__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_BJES_Response__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical2__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed1__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical5__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_negEta__1down");

  // systematics_down.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_7__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_5__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_4__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_negEta__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_2__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_1__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_6__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_3__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_EffectiveNP_8restTerm__1down");
  // systematics_down.push_back("_JET_GlobalReduction_JET_BJES_Response__1down");

  // systematics_down.push_back("_JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1down");
  // systematics_down.push_back("_JET_SR_Scenario1_JET_GroupedNP_1__1down");
  // systematics_down.push_back("_JET_SR_Scenario1_JET_GroupedNP_2__1down");
  // systematics_down.push_back("_JET_SR_Scenario1_JET_GroupedNP_3__1down");

  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical1__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed2__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_BJES_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling2__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical5__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling1__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical4__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Detector1__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed1__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling3__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical6__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical3__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical2__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_negEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed3__1down");
  // systematics_down.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling4__1down");

  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_1__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_negEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_BJES_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_3__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_2__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_5__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_7__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_8restTerm__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_6__1down");
  // systematics_down.push_back("_JET_SR_Scenario3_JET_EffectiveNP_4__1down");

  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat11__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_MC__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Jvt__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat10__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat9__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Veto__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat8__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat7__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat4__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat1__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat12__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat11__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat3__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat10__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Generator__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat6__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_GamESZee__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Jvt__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_dPhi__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat8__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat13__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat6__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Veto__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat12__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat5__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat9__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat1__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_OOC__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat3__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat2__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_dPhi__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat7__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat13__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_BJES_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Purity__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_ElecESZee__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_Stat4__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat14__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat2__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_negEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_ElecEsmear__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Zjet_KTerm__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_GamEsmear__1down");
  // systematics_down.push_back("_JET_SR_Scenario4_JET_Gjet_Stat5__1down");

  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat3__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat1__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_KTerm__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Jvt__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_ElecESZee__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_SingleParticle_HighPt__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_PunchThrough_MC16__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_highE__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat7__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat7__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_dPhi__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_GamEsmear__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Pileup_OffsetNPV__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_negEta__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Generator__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Flavor_Composition__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Purity__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Pileup_PtTerm__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat13__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat4__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Pileup_RhoTopology__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat2__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat12__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Pileup_OffsetMu__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat11__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat5__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat9__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_ElecEsmear__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat1__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_MC__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat10__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Veto__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat6__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_Modelling__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat13__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_OOC__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat11__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat8__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_BJES_Response__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat2__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat9__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_posEta__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Veto__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat10__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat3__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat6__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat12__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_TotalStat__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat5__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Flavor_Response__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat14__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat4__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_GamESZee__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat8__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_dPhi__1down");
  // systematics_down.push_back("_JET_AllNuisanceParameters_JET_Gjet_Jvt__1down");


  // systematics_up.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling1__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical5__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical4__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical3__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_negEta__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical2__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed2__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_SingleParticle_HighPt__1up");
  // //systematics_up.push_back("JET_JER_SINGLE_NP__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_PunchThrough_MC16__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling4__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical1__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Statistical6__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Detector1__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed1__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling3__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Modelling2__1up");
  // systematics_up.push_back("_JET_CategoryReduction_JET_EffectiveNP_Mixed3__1up");

  // systematics_up.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_6__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_2__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_negEta__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_8restTerm__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_SingleParticle_HighPt__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_3__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_4__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_PunchThrough_MC16__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_7__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_5__1up");
  // systematics_up.push_back("_JET_GlobalReduction_JET_EffectiveNP_1__1up");

  // systematics_up.push_back("_JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1up");
  // systematics_up.push_back("_JET_SR_Scenario1_JET_GroupedNP_1__1up");
  // systematics_up.push_back("_JET_SR_Scenario1_JET_GroupedNP_2__1up");
  // systematics_up.push_back("_JET_SR_Scenario1_JET_GroupedNP_3__1up");


  // systematics_up.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed1__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_SingleParticle_HighPt__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling3__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical3__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_PunchThrough_MC16__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Detector1__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical6__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical4__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling1__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical2__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling2__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed2__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical1__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Mixed3__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Modelling4__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EffectiveNP_Statistical5__1up");
  // //JET_JER_SINGLE_NP__1up");
  // systematics_up.push_back("_JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_negEta__1up");

  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_3__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_SingleParticle_HighPt__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_4__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_5__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_8restTerm__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_1__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_2__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_7__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EffectiveNP_6__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_negEta__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_SR_Scenario3_JET_PunchThrough_MC16__1up");

  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_MC__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Jvt__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_SingleParticle_HighPt__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_dPhi__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Veto__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_GamEsmear__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat7__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_KTerm__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat3__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat1__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat13__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat3__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat10__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat14__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat5__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_dPhi__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat8__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat11__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat11__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Veto__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat9__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat4__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat2__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat6__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat10__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat4__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat13__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat2__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_ElecEsmear__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat12__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_negEta__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_PunchThrough_MC16__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_GamESZee__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat5__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat6__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_ElecESZee__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat12__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_OOC__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat7__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat8__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Stat9__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Zjet_Stat1__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Generator__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Jvt__1up");
  // systematics_up.push_back("_JET_SR_Scenario4_JET_Gjet_Purity__1up");

  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat1__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_TotalStat__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat11__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Flavor_Response__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat7__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_KTerm__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat10__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_ElecEsmear__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Flavor_Composition__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_ElecESZee__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat11__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Jvt__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Pileup_RhoTopology__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Pileup_PtTerm__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat13__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat14__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Pileup_OffsetNPV__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Pileup_OffsetMu__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Generator__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat12__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat1__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_PunchThrough_MC16__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat3__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat2__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat7__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat6__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat13__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_Modelling__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Purity__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat3__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat2__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat5__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat4__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat8__1up");
  //systematics_up.push_back("_JET_JER_SINGLE_NP__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_dPhi__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_posEta__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat4__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat6__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat10__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_dPhi__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_SingleParticle_HighPt__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_BJES_Response__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat8__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_negEta__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_GamEsmear__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Jvt__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Veto__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_highE__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_OOC__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_GamESZee__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat12__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat5__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_MC__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Zjet_Stat9__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Stat9__1up");
  // systematics_up.push_back("_JET_AllNuisanceParameters_JET_Gjet_Veto__1up");
