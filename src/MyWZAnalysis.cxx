#define MyWZAnalysis_cxx
#include "MyWZAnalysis.h"
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include <vector>
#include <cmath>

#include <unordered_map>
#include <map>
#include <memory>

#include <iostream>
#include <string>
#include <fstream>
#include <TTree.h>
#include <TFile.h>
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include <TRandom3.h>

#include "TSystemDirectory.h"
#include "TSystemFile.h" 
#include "TString.h"
#include "TCollection.h"

#include "TLorentzVector.h"
#include "BosonTruth.h"

//RecoilCalibration Inc
#include <TProfile.h>
#include "RecoilCalibHistos.h"
#include "RecoilCalibration.h"


void MyWZAnalysis::Loop()
{
  if (fChain == 0) return;

  //set configuration
  SetConfiguration(config);
 
    bool isData=false;
    if(nameOfSample=="data16" || nameOfSample=="data15" || nameOfSample=="data17" || nameOfSample=="data18") isData=true;

    //names for output files
    std::string wzchannel;
    if(BosonZSelection){
      wzchannel="z";
    }else if(BosonWplusSelection){
      wzchannel="wplus";
    }else if(BosonWminusSelection){
      wzchannel="wminus";}

    std::string calib;
    if(SETCalibrationAnalysis && InsituCorrectionAnalysis){
      calib="_set_insitu";
    }else if(SETCalibrationAnalysis && !InsituCorrectionAnalysis){
      calib="_set";
    }else if(config.SETCalibration=="False" && InsituCorrectionAnalysis){
      calib="_insitu";
    }else if(config.SETCalibration=="False" && !InsituCorrectionAnalysis){
      calib="";
    }
    
    if(SETCalibrationAnalysis && InsituCorrectionAnalysis && ResolResponseAnalysis)
      calib="_set_insitu_resolresponse";

    std::string puname="";
    if(!OnTheFlyPileUpAnalysis) puname="_nonewPU";
    if(OnTheFlyPileUpAnalysis)  puname=config.PRWName;

    std::string multiname="", multiDir="";
    if(config.DoMultijet!="True") multiname="";
    if(DoMultijetAnalysis){
      multiDir="Multijet/";
      if(MJFitRegionSignal) multiname="_multijet_signal";
      if(MJFitRegionFitR1)  multiname="_multijet_fr1";
      if(MJFitRegionFitR2)  multiname="_multijet_fr2";
    }

    std::string d0Name="";
    if(Dod0CutAnalysis &&  !DoMultijetAnalysis) d0Name="";
    if(!Dod0CutAnalysis && !DoMultijetAnalysis) d0Name="_nod0cut";

    //if(DoMultijetAnalysis && !isData/*&& Applyd0WeightAnalysis*/) 
    Getd0sigFit(config);
    
    //Hadronic Recoil Calibration Histograms
    RecoilCalibHistos *m_recoilHistos;
    std::vector<TLorentzVector> recoLeptons;
    std::vector<int> truthCharge, recoCharge;
    TLorentzVector recoil(0.,0.,0.,0.), recoLeptons_temp_1(0.,0.,0.,0.), recoLeptons_temp_2(0.,0.,0.,0.);

    //Hadronic Recoil calibration
    RecoilCalibration *m_recoilCalibration;
    float calibWeight=1.;

    if(HasRecoilInfoAnalysis){
      //Hadronic Recoil Calibration Histograms
      m_recoilHistos = new RecoilCalibHistos();
      m_recoilHistos -> initialize(config.OutputFileDir+"RecoilHistos/",(nameOfSample+"_"+wzchannel+calib+puname).c_str());

      //Hadronic Recoil calibration
      m_recoilCalibration = new RecoilCalibration();
      m_recoilCalibration->initialize((config.RecoilCalibFileDir).c_str(),wzchannel.c_str());
    }

    //on-the-fly pileup weight calculation
    float NewPileupWeight;

    if(OnTheFlyPileUpAnalysis && !isData){      
      TChain* prwChain = new TChain("prwTree"); 
      prwChain->Add(config.PRWFile.c_str());
      prwChain->BuildIndex("weight_pileup_hash");
      fChain->AddFriend(prwChain);
      fChain->SetBranchAddress("NewPileupWeight",&NewPileupWeight);
    }

    //ouput directory according year
    std::string dirYear="";
    if(DataYears2015p2016) dirYear="2015p2016/";
    if(DataYears2017) dirYear="2017/";
    if(DataYears2018) dirYear="2018/";
    

    //skim normalisation
    float tt=0.; double tt2=0.;   
    if(!isData){
      fChain_sumWeights->SetBranchAddress("totalEventsWeighted",&tt);
      int ntree = fChain_sumWeights->GetEntries();
      for(Int_t i=0;i<ntree;i++){
	fChain_sumWeights->GetEntry(i);
	tt2=tt2+tt;
      }
    }
    float tt2p=tt2/1.;
        
    double skimNorm=1;
    
    std::cout<<" Skimming Normalisation  "<<skimNorm<<"    Sample:   "<<nameOfSample<<std::endl;
    double weightsSumWZ=0., weightsSumWZHigh=0.;
    double weightsSumWZRecoil=0., weightsSumWZRecoilHigh=0.;
    double weightsSum_1=0., weightsSum_2=0., weightsSum_3=0., weightsSum_4=0., weightsSum_5=0.;
    double weightsSum_6=0., weightsSum_7=0., weightsSum_8=0., weightsSum_9=0., weightsSum_10=0.;
    Long64_t nTruth = 0, mapCounter=0, nReco=0;
    double cutFlowIncl1 = 0, cutFlowIncl2 = 0, cutFlowIncl3 = 0, cutFlowIncl4 = 0, cutFlowIncl5 = 0, cutFlowIncl6=0, cutFlowIncl7=0, cutFlowIncl8=0;
    double cutFlowHigh1 = 0, cutFlowHigh2 = 0, cutFlowHigh3 = 0, cutFlowHigh4 = 0, cutFlowHigh5 = 0, cutFlowHigh6=0, cutFlowHigh7=0, cutFlowHigh8=0;
    
    TRandom3 *randomSmearing = new TRandom3();
    double quadraticdiff =  0.319131634;
    
    //eta binning, for multijet estimation
    std::vector<double> eta_bin_low;  eta_bin_low.clear();
    std::vector<double> eta_bin_high; eta_bin_high.clear();  
    std::vector<std::string> eta_bin_name; eta_bin_name.clear();
    eta_bin_low.push_back(10.);  eta_bin_high.push_back(40.0); eta_bin_name.push_back(""); 
    eta_bin_low.push_back(0.0);  eta_bin_high.push_back(0.21); eta_bin_name.push_back("_eta1"); 
    eta_bin_low.push_back(0.21);	eta_bin_high.push_back(0.42); eta_bin_name.push_back("_eta2");
    eta_bin_low.push_back(0.42);	eta_bin_high.push_back(0.63); eta_bin_name.push_back("_eta3");
    eta_bin_low.push_back(0.63);	eta_bin_high.push_back(0.84); eta_bin_name.push_back("_eta4");
    eta_bin_low.push_back(0.84);	eta_bin_high.push_back(1.05); eta_bin_name.push_back("_eta5");
    eta_bin_low.push_back(1.05);	eta_bin_high.push_back(1.37); eta_bin_name.push_back("_eta6");
    eta_bin_low.push_back(1.37);	eta_bin_high.push_back(1.52); eta_bin_name.push_back("_eta7");
    eta_bin_low.push_back(1.52);	eta_bin_high.push_back(1.74); eta_bin_name.push_back("_eta8");
    eta_bin_low.push_back(1.74);	eta_bin_high.push_back(1.95); eta_bin_name.push_back("_eta9");
    eta_bin_low.push_back(1.95);	eta_bin_high.push_back(2.18); eta_bin_name.push_back("_eta10");
    eta_bin_low.push_back(2.18);	eta_bin_high.push_back(2.4);  eta_bin_name.push_back("_eta11");
    
   bool sampleFlag = false;//This flag is to know for which samples the truth mass cut is needed
   if( nameOfSample == "wplusmunu" || nameOfSample == "wplustaunu" || nameOfSample == "ztautau" || nameOfSample == "zmumu") sampleFlag = true;
   if( nameOfSample == "wminmunu" || nameOfSample == "wmintaunu") sampleFlag = true;

   //truth mapping
   //creating a map for truth info and its iterator
   std::unordered_map<int,std::shared_ptr<BosonTruth>> m_myTruthMap;
   typedef std::unordered_map<int,std::shared_ptr<BosonTruth>>::const_iterator Itr;
   
   if(!isData && TruthMatchingAnalysis && !OnlyInclusiveAnalysis && !DoMultijetAnalysis && sampleFlag){

     ULong64_t T_eventNumber;
     UInt_t    T_runNumber, T_mcChannelNumber;
     
     Double_t  T_KFactor_weight;
     Float_t   T_weight_mc, T_weight_pileup;
     
     Float_t   T_WZ_m,T_WZ_pt,T_WZ_eta,T_WZ_phi;
     Float_t   T_lepton_m,T_lepton_pt,T_lepton_eta,T_lepton_phi;
     Float_t   T_antilepton_m,T_antilepton_pt,T_antilepton_eta,T_antilepton_phi;

     nTruth   = fChain_truth->GetEntries();
     fChain_truth->SetBranchAddress("eventNumber", &T_eventNumber);
     fChain_truth->SetBranchAddress("runNumber",   &T_runNumber);
     fChain_truth->SetBranchAddress("mcChannelNumber", &T_mcChannelNumber);
     
     fChain_truth->SetBranchAddress("weight_mc", &T_weight_mc);
     fChain_truth->SetBranchAddress("weight_pileup", &T_weight_pileup);
     fChain_truth->SetBranchAddress("KFactor_weight_truth", &T_KFactor_weight);
     
     fChain_truth->SetBranchAddress("MC_WZ_m",   &T_WZ_m);
     fChain_truth->SetBranchAddress("MC_WZ_pt",  &T_WZ_pt);
     fChain_truth->SetBranchAddress("MC_WZ_eta", &T_WZ_eta);
     fChain_truth->SetBranchAddress("MC_WZ_phi", &T_WZ_phi);
     
     fChain_truth->SetBranchAddress("MC_WZneutrino_m_born",   &T_lepton_m);
     fChain_truth->SetBranchAddress("MC_WZneutrino_pt_born",  &T_lepton_pt);
     fChain_truth->SetBranchAddress("MC_WZneutrino_eta_born", &T_lepton_eta);
     fChain_truth->SetBranchAddress("MC_WZneutrino_phi_born", &T_lepton_phi);
     
     fChain_truth->SetBranchAddress("MC_WZmu_el_m_born",   &T_antilepton_m);
     fChain_truth->SetBranchAddress("MC_WZmu_el_pt_born",  &T_antilepton_pt);
     fChain_truth->SetBranchAddress("MC_WZmu_el_eta_born", &T_antilepton_eta);
     fChain_truth->SetBranchAddress("MC_WZmu_el_phi_born", &T_antilepton_phi);
     
     if(config.NumberOfEvents!=-1) nTruth=config.NumberOfEvents;//custom # events
     
     for(Long64_t truth_entry=0; truth_entry <nTruth; truth_entry++){//loop over truth tree
       
       fChain_truth->GetEntry(truth_entry);

       if(truth_entry % 100000 == 0) std::cout<<"Processing truth event  "<<truth_entry<<"/"<<nTruth<<std::endl;
       
       std::shared_ptr<BosonTruth> truthEvent = std::shared_ptr<BosonTruth>(new BosonTruth());
       
       truthEvent->eventNumberTruth = T_eventNumber;
       truthEvent->runNumberTruth = T_runNumber;
       truthEvent->mcChannelNumber = T_mcChannelNumber;
       
       truthEvent->weight_mc = T_weight_mc;
       truthEvent->weight_pileup = T_weight_pileup;	    
       truthEvent->KFactor_weight_truth = T_KFactor_weight;
       	
       truthEvent->MC_WZ_m = T_WZ_m;
       truthEvent->MC_WZ_pt = T_WZ_pt;
       truthEvent->MC_WZ_eta = T_WZ_eta;
       truthEvent->MC_WZ_phi = T_WZ_phi;
       
       truthEvent->MC_lepton_m = T_lepton_m;
       truthEvent->MC_lepton_pt = T_lepton_pt;
       truthEvent->MC_lepton_eta = T_lepton_eta;
       truthEvent->MC_lepton_phi = T_lepton_phi;
       
       truthEvent->MC_antilepton_m = T_antilepton_m;
       truthEvent->MC_antilepton_pt = T_antilepton_pt;
       truthEvent->MC_antilepton_eta = T_antilepton_eta;
       truthEvent->MC_antilepton_phi = T_antilepton_phi;

       m_myTruthMap.insert(std::make_pair(T_eventNumber,truthEvent));//fill the map

      }//truth loop
    }//if not Data and sampleFlag

   //create a loop for every systematic given by scale factors
   Int_t SFloopSize = 0;
   std::vector<std::string> SF_variations;//Scale Factors
   SF_variations.push_back("nominal");
   if(SFVariationsAnalysis && nameOfSystematic == "nominal"){
     SF_variations.push_back("IDStatD");
     SF_variations.push_back("IDSysD");
     SF_variations.push_back("IDStatU");
     SF_variations.push_back("IDSysU");
     SF_variations.push_back("IsoStatD");
     SF_variations.push_back("IsoSysD");
     SF_variations.push_back("IsoStatU");
     SF_variations.push_back("IsoSysU");
     SF_variations.push_back("TrigStatD");
     SF_variations.push_back("TrigSysD");
     SF_variations.push_back("TrigStatU");
     SF_variations.push_back("TrigSysU");
     SF_variations.push_back("TTVAStatD");
     SF_variations.push_back("TTVASysD");
     SF_variations.push_back("TTVAStatU");
     SF_variations.push_back("TTVASysU");
     SF_variations.push_back("PileUpUp");
     SF_variations.push_back("PileUpDown");

     SF_variations.push_back("JVTUp");
     SF_variations.push_back("JVTDown");

     SF_variations.push_back("kfactor_ALPHAS__1down");
     SF_variations.push_back("kfactor_ALPHAS__1up");
     SF_variations.push_back("kfactor_BEAM_ENERGY__1down");
     SF_variations.push_back("kfactor_BEAM_ENERGY__1up");
     SF_variations.push_back("kfactor_CHOICE_HERAPDF20");
     SF_variations.push_back("kfactor_CHOICE_NNPDF30");
     SF_variations.push_back("kfactor_PDF_EV1");
     SF_variations.push_back("kfactor_PDF_EV2");
     SF_variations.push_back("kfactor_PDF_EV3");
     SF_variations.push_back("kfactor_PDF_EV4");
     SF_variations.push_back("kfactor_PDF_EV5");
     SF_variations.push_back("kfactor_PDF_EV6");
     SF_variations.push_back("kfactor_PDF_EV7");
     SF_variations.push_back("kfactor_PDF_EW__1down");
     SF_variations.push_back("kfactor_PDF_EW__1up");
     SF_variations.push_back("kfactor_PDF__1down");
     SF_variations.push_back("kfactor_PDF__1up");
     SF_variations.push_back("kfactor_PI__1down");
     SF_variations.push_back("kfactor_PI__1up");
     SF_variations.push_back("kfactor_REDCHOICE_NNPDF30");
     SF_variations.push_back("kfactor_SCALE_W__1down");
     SF_variations.push_back("kfactor_SCALE_W__1up");
     SF_variations.push_back("kfactor_SCALE_Z__1down");
     SF_variations.push_back("kfactor_SCALE_Z__1up");
   }
   
    if(nameOfSystematic == "nominal" &&  !isData){
      SFloopSize = (int)SF_variations.size();
    }else  SFloopSize = 1;


    for(int iSF_variations = 0; iSF_variations < SFloopSize; ++iSF_variations){//Scale Factor loop
     
      if(nameOfSystematic == "nominal" && !isData){
	newfile.push_back( new TFile((config.OutputFileDir + "Files/" + dirYear + multiDir + nameOfSample + "_"+ wzchannel +  "_" + SF_variations[iSF_variations] + calib + puname +multiname+d0Name+".root").c_str(),"RECREATE") );
      }else if(nameOfSystematic != "nominal" && !isData){
	newfile.push_back( new TFile((config.OutputFileDir +"Files/" + dirYear + multiDir + nameOfSample +"_"+ wzchannel + "_" + nameOfSystematic + calib + puname + multiname+d0Name+".root").c_str(),"RECREATE") );
      }else if (nameOfSystematic == "nominal" && isData){
	newfile.push_back( new TFile((config.OutputFileDir + "Files/" + dirYear + multiDir + nameOfSample +"_"+ wzchannel + multiname + d0Name+ ".root").c_str(),"RECREATE") );
      }

    }//scale factor loop

    Long64_t nentries = fChain->GetEntries();
    
    gStyle->SetOptStat("nemruo");
    TH1::SetDefaultSumw2(kTRUE);
    gROOT->ForceStyle();
    
    //Load Histograms
    if(!DoMultijetAnalysis){
      for(int iSF_variations = 0; iSF_variations < SFloopSize; ++iSF_variations)//Scale Factor loop
	setHist(config,SF_variations[iSF_variations]);
    }else if(DoMultijetAnalysis){ 
      for(int iBin=0; iBin<(int)eta_bin_name.size();iBin++)
	setHistMultijet(config,eta_bin_name[iBin]);
    }
    
    //reco variables
    pt=0.; lpt=0.; calo20=0.; track20=0.; eta=0.; phi=0.; dPhi=0.;  z0=0.;
    pt2=0.;  calo202=0.; track202=0.; eta2=0.; phi2=0.; dPhi2=0.; z02=0.;
    iso30=0.; d0sig=0.; d0sig2=0.;
    met=0.; met_recoil=0.; phimet=0.; lrmass=0.; recoMass=0.; recoMass_recoil=0.;
    pileup = 1.; skimming = 1.; kFactor=1.; mc_weight=1.;
    SF=1.; SF_isol=1.; SF_trig=1.; SF_ID=1.; SF_ttva=1.;
    d0sigWeight=1.;
    mu_scaled=0.;
    px=0.; py=0.;
    px_met=0.; py_met=0.;
    px_sum=0.; py_sum=0.;
    pT_tot=0.;
    
    //Cw factor error calculation
    int iRecBin = 0, iGenBin=0;
    int iRecBinEta = 0, iGenBinEta=0;
    //define lepton charge
    float chargeSelection = 0.;
    if(config.WZSelection == "wplus"){
      chargeSelection = 1.;
    }else if(config.WZSelection == "wminus"){
      chargeSelection = -1.;
    }
        
    Long64_t nbytes = 0, nb = 0;
    
    if(config.NumberOfEvents!=-1) nentries=config.NumberOfEvents;//custom # events
    nReco= nentries;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      
      if(jentry % 100000 == 0) std::cout<<"Processing event  "<<jentry<<"/"<<nentries<<std::endl;
      
      bool truthmassMatch  = false;
      bool truthSelection  = false;
      TLorentzVector truthWZboson(0.,0.,0.,0.), diLepton(0.,0.,0.,0.);	
      
      //truth-reco event matching
      if(TruthMatchingAnalysis && !OnlyInclusiveAnalysis && !DoMultijetAnalysis && sampleFlag){
	Itr i = m_myTruthMap.find(eventNumber);
	if(i != m_myTruthMap.end()){ 
	  std::shared_ptr<BosonTruth> myTruthEvent = (*i).second;
	  mapCounter=mapCounter+1;
	  
	  truthWZboson.SetPtEtaPhiM(myTruthEvent->MC_WZ_pt,myTruthEvent->MC_WZ_eta,myTruthEvent->MC_WZ_phi,myTruthEvent->MC_WZ_m);
	  
	  TLorentzVector temp1(0.,0.,0.,0.), temp2(0.,0.,0.,0.);
	  temp1.SetPtEtaPhiM(myTruthEvent->MC_lepton_pt,myTruthEvent->MC_lepton_eta,myTruthEvent->MC_lepton_phi,myTruthEvent->MC_lepton_m);
	  temp2.SetPtEtaPhiM(myTruthEvent->MC_antilepton_pt,myTruthEvent->MC_antilepton_eta,myTruthEvent->MC_antilepton_phi,myTruthEvent->MC_antilepton_m);
	  diLepton=temp1+temp2;
	  
	  if( myTruthEvent->MC_WZ_m < 120000.) truthmassMatch = true;//avoid double counting
	  
	  //Cw factor calculation
	  float ptvTruth  = myTruthEvent->MC_lepton_pt;//neutrino
	  float etavTruth = myTruthEvent->MC_lepton_eta;
	  float phivTruth = myTruthEvent->MC_lepton_phi;
	  
	  float ptlTruth  = myTruthEvent->MC_antilepton_pt;
	  float etalTruth = myTruthEvent->MC_antilepton_eta;
	  float philTruth = myTruthEvent->MC_antilepton_phi;
	  
	  float dPhiTruth=std::fabs(phivTruth - philTruth);
	  if (dPhiTruth>3.14159265) dPhiTruth=2.*3.14159265-dPhiTruth;
	  double wtmTruth = std::sqrt( 2.*(ptlTruth*ptvTruth)*(1.-std::cos(dPhiTruth)));
	  
	  iGenBin = hMassCgen->FindBin(wtmTruth);
	  iGenBinEta = hEtaCgen->FindBin(std::fabs(etalTruth));
	  
	  if(ptlTruth > cut_pt1 && ptvTruth > cut_met1){
	    if(std::fabs(etalTruth) < cut_eta /*&& std::fabs(etavTruth) < cut_eta*/){
	      if(wtmTruth > cut_mwt1){
		truthSelection = true;
	      }
	    }
	  }
	  
	} //check reco event has a truth event associated with it 
      }//sample flag, truth matching
      
      if( !sampleFlag || isData )  truthmassMatch = true;
      if( OnlyInclusiveAnalysis ) truthmassMatch = true;
      
      SF_isol=weight_indiv_SF_MU_Isol;
      //SF_trig=weight_indiv_SF_MU_Trigger;
      SF_trig=weight_globalLeptonTriggerSF;
      SF_ID=weight_indiv_SF_MU_ID;
      SF_ttva=weight_indiv_SF_MU_TTVA;
      SF=SF_isol*SF_trig*SF_ID*SF_ttva;//final eff scale factor
      
      if(/*iSF_variations == 0 && nameOfSystematic != "nominal"*/ config.Systematics=="True")//In the systematic trees there is only one SF
	if(!isData) SF=1.;//weight_oldTriggerSF;//SF=weight_leptonSF; //change this, bug coming from AT!
      
      if(isData){ SF=1.; skimming = 1.; pileup = 1.; kFactor = 1.; mc_weight=1.; weight_leptonSF =1.; mu_scaled=mu/1.03; lumi_sample=1.;} 
      if(!isData){ skimming = skimNorm;  mu_scaled=mu; mc_weight=weight_mc; pileup = weight_pileup;}
            
      if(nameOfSample == "ttbar" || nameOfSample == "wt_top"  || nameOfSample == "wt_antitop"
	 || nameOfSample == "single_top"  || nameOfSample == "single_antitop"  || nameOfSample == "zzqqll"
	 || nameOfSample == "wzqqll"  || nameOfSample == "wzlnuqq"  || nameOfSample == "zzllnunu"
	 || nameOfSample == "zwlllnu_SFMinus"  || nameOfSample == "zwlllnu_OFMinus"  || nameOfSample == "zwlllnu_SFPlus" 
	 || nameOfSample == "zwlllnu_OFPlus"  || nameOfSample == "zzllll" ){
	  kFactor = 1.;
	}else if(!isData){
	  kFactor =  weight_KFactor;
	}
	
	//wzselection depends on year/dataset
	bool WZSelection=false;
	if(BosonWminusSelection || BosonWplusSelection){
	  if(DataYears2015p2016) WZSelection = munu_2015==1 || munu_2016==1;
	  else if(DataYears2017)      WZSelection = munu_2017==1;
	  else if(DataYears2018)      WZSelection = munu_2018==1;
	}	
	else if(BosonZSelection){
	  if(DataYears2015p2016) WZSelection = mumu_2015==1 || mumu_2016==1;
	  if(DataYears2017) WZSelection = mumu_2017==1;
	  if(DataYears2018) WZSelection = mumu_2018==1;
	}	  


	if(WZSelection && truthmassMatch){
	  
	basicZEvent=false; basicWEvent=false;
	goodZEvent=false; goodWEvent=false;
	goodZEventHigh=false; goodWEventHigh=false;
	goodWEventRecoil=false; goodWEventRecoilHigh=false;
	Zboson.SetPtEtaPhiM(0.,0.,0.,0.); Zboson1.SetPtEtaPhiM(0.,0.,0.,0.);
	Zboson2.SetPtEtaPhiM(0.,0.,0.,0.);
	Wboson.SetPtEtaPhiM(0.,0.,0.,0.);  wtm=0.;
	tboson.SetPtEtaPhiM(0.,0.,0.,0.); upfo.SetPtEtaPhiM(0.,0.,0.,0.);
	upfo_high.SetPtEtaPhiM(0.,0.,0.,0.);
	d0sigCut=true;
	
	tboson = truthWZboson;
	if(!isData && OnTheFlyPileUpAnalysis) pileup = NewPileupWeight;
	
	//ZSelection
	if(BosonZSelection){
 
	    if(mu_pt->size()>=2 && el_pt->size()==0){   
	      Zboson1.SetPtEtaPhiE(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_e->at(0));
	      Zboson2.SetPtEtaPhiE(mu_pt->at(1),mu_eta->at(1),mu_phi->at(1),mu_e->at(1));
	      Zboson=Zboson1+Zboson2;
	      if(mu_charge->at(0)!=mu_charge->at(1)){
		if( fabs(mu_eta->at(0)) < cut_eta && fabs(mu_eta->at(1)) < cut_eta){
		  basicZEvent=true;		  
		}
	      }
	    }
	    
	    if(basicZEvent){
	      if(mu_pt->at(0) > cut_pt1 || mu_pt->at(1) > cut_pt1){
		if(Zboson.M() > 66000. && Zboson.M()< 116000.){
		  goodZEvent=true;
		  //weightsSumWZ+=SF*pileup*kFactor*mc_weight;
		  weightsSumWZ+=pileup*mc_weight;	  
		}
	      }
	      
	      if(mu_pt->at(0) > cut_pt2 || mu_pt->at(1) > cut_pt2){		      
		if(Zboson.M() > 66000. && Zboson.M()< 116000.){
		  goodZEventHigh=true;
		  weightsSumWZHigh+=pileup*mc_weight;
		}
	      }		
	    }//basic Z event
			   	  
	  }//Zselection

	  //WSelection
	  if(BosonWminusSelection || BosonWplusSelection){

	    if(mu_pt->size()==1 && el_pt->size()==0){
	      //cutFlowIncl1=cutFlowIncl1+1; cutFlowHigh1=cutFlowHigh1+1;
	      cutFlowIncl1+=pileup*mc_weight; cutFlowHigh1+=pileup*mc_weight;
	      float dPhi=std::fabs(met_phi-mu_phi->at(0));
	      if (dPhi>3.14159265) dPhi=2.*3.14159265-dPhi;
	      Wboson.SetPtEtaPhiE(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_e->at(0));	   
	      wtm = std::sqrt(2.*((mu_pt->at(0))*met_met)*(1.-std::cos(dPhi)));
	      //cw factor calc
	      if(!DoMultijetAnalysis) iRecBin    = hMassCrec->FindBin(wtm);
	      if(!DoMultijetAnalysis) iRecBinEta = hEtaCrec->FindBin(std::abs(mu_eta->at(0)));

	      if((mu_charge->at(0)) == chargeSelection){
		//cutFlowIncl2=cutFlowIncl2+1; cutFlowHigh2=cutFlowHigh2+1;
		cutFlowIncl2+=pileup*mc_weight; cutFlowHigh2+=pileup*mc_weight;
		if(std::fabs(mu_eta->at(0)) < cut_eta && std::fabs(mu_delta_z0_sintheta->at(0))<0.5){
		  //cutFlowIncl3=cutFlowIncl3+1; cutFlowHigh3=cutFlowHigh3+1;
		  cutFlowIncl3+=pileup*mc_weight; cutFlowHigh3+=pileup*mc_weight;
		  basicWEvent=true;
		}
	      }
	    }
	
	    if(basicWEvent){
	      d0sig=mu_d0sig->at(0); d0sigCut=true;
	      if(isData && config.Nod0Shift=="True") shift=0.;
	      d0sig = d0sig + shift;// + randomSmearing->Gaus(0,quadraticdiff);

	      //d0sig reweighting
	      if(Applyd0WeightAnalysis && !isData){

		d0sigWeight = fit_d0sig_gaussratio->Eval(d0sig);

	        double sigmaVar=0.;
		if(config.MultijetVariation=="20" && config.DoMultijet=="True") sigmaVar = -0.5;
		if(config.MultijetVariation=="30" && config.DoMultijet=="True") sigmaVar =  0.5;
		//if(d0sigWeight>2.)   d0sigWeight = 2.;
		if(d0sig >= 5.+sigmaVar)
		  d0sigWeight = fit_d0sig_gaussratio->Eval(5.+sigmaVar);
		else if(d0sig <= (-1.)*(5.+sigmaVar))
		  d0sigWeight = fit_d0sig_gaussratio->Eval((-1.)*(5.+sigmaVar));

		if(d0sigWeight<0.001 || std::isnan(d0sigWeight)|| std::isinf(d0sigWeight)) d0sigWeight = 1.;
	      }else
		d0sigWeight = 1;
		     
	      if(!DoMultijetAnalysis && Dod0CutAnalysis) d0sigCut = std::fabs(d0sig) < 3. && std::fabs(mu_delta_z0_sintheta->at(0))<0.5;

	      if(d0sigCut){
		if(mu_pt->at(0) > cut_pt1){
		  //cutFlowIncl4=cutFlowIncl4+1;
		  cutFlowIncl4+=pileup*mc_weight;
		  if(met_met > cut_met1){
		    //cutFlowIncl5=cutFlowIncl5+1;
		    cutFlowIncl5+=pileup*mc_weight;
		    if(wtm > cut_mwt1){
		      //cutFlowIncl6=cutFlowIncl6+1;
		      cutFlowIncl6+=pileup*mc_weight;
		      goodWEvent=true;
		      weightsSumWZ+=pileup*mc_weight;
		      //weightsSum_1+=pileup;
		      //weightsSum_2+=mc_weight;
		      weightsSum_3+=SF; weightsSum_4+=kFactor; weightsSum_5+=SF_ID;
		      weightsSum_6+=SF_isol;// weightsSum_7+=SF_trig;
		    }
		  }
		}
	      }
		
	      if(mu_pt->at(0) > cut_pt2 && d0sigCut){
		//cutFlowHigh4=cutFlowHigh4+1;
		cutFlowHigh4+=pileup*mc_weight;
		if(met_met > cut_met2){
		  //cutFlowHigh5=cutFlowHigh5+1;
		  cutFlowHigh5+=pileup*mc_weight;
		  if(wtm > cut_mwt2){
		    //cutFlowHigh6=cutFlowHigh6+1;
		    cutFlowHigh6+=pileup*mc_weight;
		    goodWEventHigh=true;
		    weightsSumWZHigh+=pileup*mc_weight;
		  }
		}
	      }

	      //peak region recoil
	      if(mu_pt->at(0)>cut_pt1){
		TLorentzVector tempRec(0.,0.,0.,0.), lepTemp(0.,0.,0.,0.);
		tempRec.SetPtEtaPhiE(hadronic_pt,hadronic_eta,hadronic_phi,hadronic_e);
		lepTemp.SetPtEtaPhiE(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_e->at(0));
		tempRec+=lepTemp;
		tempRec=-1*tempRec;
		double dPhiRec=std::fabs(tempRec.Phi()-mu_phi->at(0));
		if(dPhiRec>3.14159265) dPhiRec=2.*3.14159265-dPhiRec;
		double RecoilMass = std::sqrt(2.*((mu_pt->at(0))*tempRec.Pt())*(1.-std::cos(dPhiRec)));
		if(tempRec.Pt()>cut_met1){
		  if(RecoilMass>cut_mwt1){
		    goodWEventRecoil=true;
		    weightsSumWZRecoil+=pileup*mc_weight;
		  }
		}
	      }

	      //high mass recoil
	      if(mu_pt->at(0)>cut_pt2){
		TLorentzVector tempRec(0.,0.,0.,0.), lepTemp(0.,0.,0.,0.);
		tempRec.SetPtEtaPhiE(hadronic_pt,hadronic_eta,hadronic_phi,hadronic_e);
		lepTemp.SetPtEtaPhiE(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_e->at(0));
		tempRec+=lepTemp;
		tempRec=-1*tempRec;
		double dPhiRec=std::fabs(tempRec.Phi()-mu_phi->at(0));
		if(dPhiRec>3.14159265) dPhiRec=2.*3.14159265-dPhiRec;
		double RecoilMass = std::sqrt(2.*((mu_pt->at(0))*tempRec.Pt())*(1.-std::cos(dPhiRec)));
		if(tempRec.Pt()>cut_met2){
		  if(RecoilMass>cut_mwt2){
		    goodWEventRecoilHigh=true;
		    weightsSumWZRecoilHigh+=pileup*mc_weight;
		  }
		}
	      }

	    }//basic W event
	  }//WSelection

	  if(HasRecoilInfoAnalysis /*&& SF_variations[iSF_variations]=="nominal"*/ && nameOfSystematic== "nominal"){
	    if(goodZEvent || goodWEventRecoil){
	      
	      recoil.SetPtEtaPhiE(hadronic_pt,hadronic_eta,hadronic_phi,hadronic_e);
	      recoLeptons_temp_1.SetPtEtaPhiE(mu_pt->at(0),mu_eta->at(0),mu_phi->at(0),mu_e->at(0));
	      recoLeptons.push_back(recoLeptons_temp_1);
	      if(basicZEvent){
		recoLeptons_temp_2.SetPtEtaPhiE(mu_pt->at(1),mu_eta->at(1),mu_phi->at(1),mu_e->at(1));
		recoLeptons.push_back(recoLeptons_temp_2);}
	      
	      upfo=recoil;	      
	      double sumet = sumET_PFO,  recoWeight=SF*pileup*kFactor*mc_weight*skimming;
	      
	      TLorentzVector calibBoson(0.,0.,0.,0.);
	      if(goodZEvent)       calibBoson=Zboson;
	      if(goodWEventRecoil) calibBoson=Wboson;

	      //set calibration
	      if(SETCalibrationAnalysis)
		m_recoilCalibration->applySumetReweighting(calibBoson.Pt()/1000.,sumet/1000.,calibWeight);
	      
	      //insitu correction
	      if(!isData && InsituCorrectionAnalysis)
		m_recoilCalibration->applyUX(sumet/1000., upfo);
	      if(!isData && InsituCorrectionAnalysis)
		m_recoilCalibration->applyUY(sumet/1000., upfo);

	      //response and resolution correction
	      if(!isData && ResolResponseAnalysis)
		m_recoilCalibration->applyResponseAndResolCorr(tboson/*calibBoson*/,sumet/1000.,upfo);

	      //final recoil histos  
	      if(!SETCalibrationAnalysis || isData ) calibWeight=1.;
	      m_recoilHistos->execute(tboson,recoLeptons,upfo,sumet,mu,recoWeight*calibWeight);	
	      
	    }//good WZ event
	  }//do recoil calibration
	  
	  if(basicZEvent || basicWEvent){
	    
	    if(!HasRecoilInfoAnalysis) calibWeight=1.;
	    
	    //setup our variables
	    met=met_met; phimet=met_phi;	   

	    pt=mu_pt->at(0); eta=mu_eta->at(0); phi=mu_phi->at(0); z0=mu_delta_z0_sintheta->at(0);
	    calo20=mu_topoetcone20->at(0);
	    lpt=std::log10(pt*0.001);
	    //track20->mu_ptvamWTrcone30->at(1);
	    dPhi=std::fabs(phimet-phi);
	    if(dPhi>3.14159265) dPhi=2.*3.14159265-dPhi;

	    px=pt*std::cos(phi); py=pt*std::sin(phi);
	    if(basicWEvent){px_met=met*std::cos(met_phi); py_met=met*std::sin(met_phi);}
	    if(basicZEvent){px_met=(mu_pt->at(1))*std::cos((mu_phi->at(1))); py_met=(mu_pt->at(1))*std::sin((mu_phi->at(1)));}

	    px_sum=px+px_met;  py_sum=py+py_met;
	    pT_tot=std::sqrt(px_sum*px_sum+py_sum*py_sum);

	    iso30=mu_ptvarcone30->at(0)/pt;
	    d0sig=mu_d0sig->at(0)+shift;

	    bool passFitRegion=false;
	    if(DoMultijetAnalysis){
	      if(MJFitRegionSignal) passFitRegion = goodWEvent;
	      if(MJFitRegionFitR1)  passFitRegion = pT_tot<30000. && pt>cut_pt1 && basicWEvent;
	      if(MJFitRegionFitR2)  passFitRegion = pt>cut_pt1 && basicWEvent;
	    }

	    //recoil variales 
	    TLorentzVector v_recoil_met(0.,0.,0.,0.), lep1(0.,0.,0.,0.), lep2(0.,0.,0.,0.);
	    v_recoil_met.SetPtEtaPhiE(hadronic_pt,hadronic_eta,hadronic_phi,hadronic_e);
	    lep1.SetPtEtaPhiE(pt,eta,phi,mu_e->at(0));

	    if(goodZEvent){
	      recoMass=Zboson.M();
	      pt2=mu_pt->at(1); eta2=mu_eta->at(1); phi2=mu_phi->at(1); z02=mu_delta_z0_sintheta->at(1);
	      d0sig2=mu_d0sig->at(1)+shift; calo202=mu_topoetcone20->at(1); //track202->mu_ptvarcone30->at(1);
	      dPhi2=std::fabs(phimet-phi2);
	      if(dPhi2>3.14159265) dPhi2=2.*3.14159265-dPhi2;
	      lep2.SetPtEtaPhiE(pt2,eta2,phi2,mu_e->at(1));
	      v_recoil_met+=lep1+lep2;
	      v_recoil_met= -1* v_recoil_met;
	      met_recoil=v_recoil_met.Pt();	      
	    }

	    if(goodWEvent) recoMass=wtm;//wtm transversmass contructed before any selection

	    lrmass=std::log10(recoMass*0.001);

	    //pileup=1.;
	    //skimming=1.;
	    //SF=weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA;
	    //SF=1.;
	    // if(DoMultijetAnalysis) SF=1.;
	    // SF=1.;
	    //kFactor=1;
	    finalWeight=SF*pileup*kFactor*mc_weight*d0sigWeight*weight_jvt;//*calibWeight*skimming;
	    if(isData) finalWeight=1.;

	    //trigger normalisation
	    weightsSum_1+=(SF*pileup*kFactor*mc_weight*d0sigWeight*weight_jvt)/weight_globalLeptonTriggerSF;
	    weightsSum_2+=SF*pileup*kFactor*mc_weight*d0sigWeight*weight_jvt;

	    //pileup weight normalisation
	    weightsSum_7+=SF*kFactor*mc_weight*d0sigWeight*weight_jvt;
	    weightsSum_8+=SF*pileup*kFactor*mc_weight*d0sigWeight*weight_jvt;

	    //d0sigWeight normalisation	(pending)
	    weightsSum_9 +=SF*pileup*kFactor*mc_weight*weight_jvt;
	    weightsSum_10+=SF*pileup*kFactor*mc_weight*weight_jvt*d0sigWeight;

	    std::vector<double> variation; variation.clear();
	    variation.push_back(1.);
	    if(SFVariationsAnalysis && nameOfSystematic == "nominal"){
	      variation.push_back(weight_indiv_SF_MU_ID_STAT_DOWN/weight_indiv_SF_MU_ID);
	      variation.push_back(weight_indiv_SF_MU_ID_SYST_DOWN/weight_indiv_SF_MU_ID);
	      variation.push_back(weight_indiv_SF_MU_ID_STAT_UP/weight_indiv_SF_MU_ID);
	      variation.push_back(weight_indiv_SF_MU_ID_SYST_UP/weight_indiv_SF_MU_ID);

	      variation.push_back(weight_indiv_SF_MU_Isol_STAT_DOWN/weight_indiv_SF_MU_Isol);
	      variation.push_back(weight_indiv_SF_MU_Isol_SYST_DOWN/weight_indiv_SF_MU_Isol);
	      variation.push_back(weight_indiv_SF_MU_Isol_STAT_UP/weight_indiv_SF_MU_Isol);
	      variation.push_back(weight_indiv_SF_MU_Isol_SYST_UP/weight_indiv_SF_MU_Isol);

	      variation.push_back(weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN/weight_globalLeptonTriggerSF);
	      variation.push_back(weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN/weight_globalLeptonTriggerSF);
	      variation.push_back(weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP/  weight_globalLeptonTriggerSF);
	      variation.push_back(weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP/  weight_globalLeptonTriggerSF);

	      variation.push_back(weight_indiv_SF_MU_TTVA_STAT_DOWN/weight_indiv_SF_MU_TTVA);
	      variation.push_back(weight_indiv_SF_MU_TTVA_SYST_DOWN/weight_indiv_SF_MU_TTVA);
	      variation.push_back(weight_indiv_SF_MU_TTVA_STAT_UP/weight_indiv_SF_MU_TTVA);
	      variation.push_back(weight_indiv_SF_MU_TTVA_SYST_UP/weight_indiv_SF_MU_TTVA);

	      variation.push_back(weight_pileup_UP/weight_pileup);
	      variation.push_back(weight_pileup_DOWN/weight_pileup);

	      variation.push_back(weight_jvt_UP/weight_jvt);
	      variation.push_back(weight_jvt_DOWN/weight_jvt);
	      //sys kfactor
	      variation.push_back(kFactor/weight_kfactor_sys_ALPHAS__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_ALPHAS__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_BEAM_ENERGY__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_BEAM_ENERGY__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_CHOICE_HERAPDF20);
	      variation.push_back(kFactor/weight_kfactor_sys_CHOICE_NNPDF30);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV1);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV2);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV3);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV4);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV5);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV6);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EV7);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EW__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF_EW__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_PDF__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_PI__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_PI__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_REDCHOICE_NNPDF30);
	      variation.push_back(kFactor/weight_kfactor_sys_SCALE_W__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_SCALE_W__1up);
	      variation.push_back(kFactor/weight_kfactor_sys_SCALE_Z__1down);
	      variation.push_back(kFactor/weight_kfactor_sys_SCALE_Z__1up);
	    }
	    
	    if(/*!isData && */!DoMultijetAnalysis)//Scale factors systematics, creates sys files
	      for(int iSF_variations = 0; iSF_variations < SFloopSize; ++iSF_variations)//Scale Factor loop
	    	FillHistos(iSF_variations, variation[iSF_variations]);
	  	    
	    //Cw factor calculation
	    if(!DoMultijetAnalysis){
	      hMassCrec->Fill(recoMass,finalWeight);
	      if( iGenBin == iRecBin && truthSelection  ) hMassCstay->Fill(recoMass,finalWeight);
	      if( iGenBin == iRecBin && !truthSelection ) hMassCcome->Fill(recoMass,finalWeight);
	      if( iGenBin != iRecBin                    ) hMassCcome->Fill(recoMass,finalWeight);
	      
	      hEtaCrec->Fill(std::abs(eta),finalWeight);
	      //std::cout<<iGenBinEta<<"    "<<iRecBinEta<<"    "<<truthSelection<<std::endl;
	      if( iGenBinEta == iRecBinEta && truthSelection  ) hEtaCstay->Fill(std::abs(eta),finalWeight);
	      if( iGenBinEta == iRecBinEta && !truthSelection ) hEtaCcome->Fill(std::abs(eta),finalWeight);
	      if( iGenBinEta != iRecBinEta                    ) hEtaCcome->Fill(std::abs(eta),finalWeight);
	    }

	    if((goodWEvent || goodZEvent) && !DoMultijetAnalysis && nameOfSystematic== "nominal"){
	      p_mwt_Tot->Fill(recoMass,finalWeight);								     
	      p_mwt_Rec->Fill(recoMass,weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);    	     

	      p_mwt_KF ->Fill(recoMass,kFactor);
	      p_mwt_MC ->Fill(recoMass,mc_weight);    

	      p_mwt_Iso->Fill(recoMass,weight_indiv_SF_MU_Isol); 
	      p_mwt_Iso_sys_UP->Fill(recoMass,weight_indiv_SF_MU_Isol_SYST_UP);
	      p_mwt_Iso_sys_DOWN->Fill(recoMass,weight_indiv_SF_MU_Isol_SYST_DOWN);
	      p_mwt_Iso_stat_UP->Fill(recoMass,weight_indiv_SF_MU_Isol_STAT_UP);
	      p_mwt_Iso_stat_DOWN->Fill(recoMass,weight_indiv_SF_MU_Isol_STAT_DOWN);

	      p_mwt_ID ->Fill(recoMass,weight_indiv_SF_MU_ID);
	      p_mwt_ID_sys_UP ->Fill(recoMass,weight_indiv_SF_MU_ID_SYST_UP);
	      p_mwt_ID_sys_DOWN ->Fill(recoMass,weight_indiv_SF_MU_ID_SYST_DOWN);
	      p_mwt_ID_stat_UP ->Fill(recoMass,weight_indiv_SF_MU_ID_STAT_UP);
	      p_mwt_ID_stat_DOWN ->Fill(recoMass,weight_indiv_SF_MU_ID_STAT_DOWN);

	      p_mwt_Tri->Fill(recoMass,weight_globalLeptonTriggerSF);
	      p_mwt_Tri_sys_UP->Fill(recoMass,weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP);
	      p_mwt_Tri_sys_DOWN->Fill(recoMass,weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN);
	      p_mwt_Tri_stat_UP->Fill(recoMass,weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP);
	      p_mwt_Tri_stat_DOWN->Fill(recoMass,weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN);

	      p_mwt_ttv->Fill(recoMass,weight_indiv_SF_MU_TTVA);
	      p_mwt_ttv_sys_UP->Fill(recoMass,weight_indiv_SF_MU_TTVA_SYST_UP);
	      p_mwt_ttv_sys_DOWN->Fill(recoMass,weight_indiv_SF_MU_TTVA_SYST_DOWN);
	      p_mwt_ttv_stat_UP->Fill(recoMass,weight_indiv_SF_MU_TTVA_STAT_UP);
	      p_mwt_ttv_stat_DOWN->Fill(recoMass,weight_indiv_SF_MU_TTVA_STAT_DOWN);

	      p_mwt_PU ->Fill(recoMass,pileup);
	      p_mwt_PU_UP ->Fill(recoMass,weight_pileup_UP);
	      p_mwt_PU_DOWN ->Fill(recoMass,weight_pileup_DOWN);

	      p_mwt_jvt->Fill(recoMass,weight_jvt);    
	      p_mwt_jvt_UP->Fill(recoMass,weight_jvt_UP);    
	      p_mwt_jvt_DOWN->Fill(recoMass,weight_jvt_DOWN);    

	      p_m_Tot->Fill(bornMass_KFactor,finalWeight);											     
	      p_m_Rec->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);    	     

	      p_m_KF ->Fill(bornMass_KFactor,kFactor);
	      p_m_MC ->Fill(bornMass_KFactor,mc_weight);    

	      p_m_Iso->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol); 
	      p_m_Iso_sys_UP->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol_SYST_UP);
	      p_m_Iso_sys_DOWN->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol_SYST_DOWN);
	      p_m_Iso_stat_UP->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol_STAT_UP);
	      p_m_Iso_stat_DOWN->Fill(bornMass_KFactor,weight_indiv_SF_MU_Isol_STAT_DOWN);

	      p_m_ID ->Fill(bornMass_KFactor,weight_indiv_SF_MU_ID);
	      p_m_ID_sys_UP ->Fill(bornMass_KFactor,weight_indiv_SF_MU_ID_SYST_UP);
	      p_m_ID_sys_DOWN ->Fill(bornMass_KFactor,weight_indiv_SF_MU_ID_SYST_DOWN);
	      p_m_ID_stat_UP ->Fill(bornMass_KFactor,weight_indiv_SF_MU_ID_STAT_UP);
	      p_m_ID_stat_DOWN ->Fill(bornMass_KFactor,weight_indiv_SF_MU_ID_STAT_DOWN);

	      p_m_Tri->Fill(bornMass_KFactor,weight_globalLeptonTriggerSF);
	      p_m_Tri_sys_UP->Fill(bornMass_KFactor,weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP);
	      p_m_Tri_sys_DOWN->Fill(bornMass_KFactor,weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN);
	      p_m_Tri_stat_UP->Fill(bornMass_KFactor,weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP);
	      p_m_Tri_stat_DOWN->Fill(bornMass_KFactor,weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN);

	      p_m_ttv->Fill(bornMass_KFactor,weight_indiv_SF_MU_TTVA);
	      p_m_ttv_sys_UP->Fill(bornMass_KFactor,weight_indiv_SF_MU_TTVA_SYST_UP);
	      p_m_ttv_sys_DOWN->Fill(bornMass_KFactor,weight_indiv_SF_MU_TTVA_SYST_DOWN);
	      p_m_ttv_stat_UP->Fill(bornMass_KFactor,weight_indiv_SF_MU_TTVA_STAT_UP);
	      p_m_ttv_stat_DOWN->Fill(bornMass_KFactor,weight_indiv_SF_MU_TTVA_STAT_DOWN);

	      p_m_PU ->Fill(bornMass_KFactor,pileup);
	      p_m_PU_UP ->Fill(bornMass_KFactor,weight_pileup_UP);
	      p_m_PU_DOWN ->Fill(bornMass_KFactor,weight_pileup_DOWN);

	      p_m_jvt->Fill(bornMass_KFactor,weight_jvt);    
	      p_m_jvt_UP->Fill(bornMass_KFactor,weight_jvt_UP);    
	      p_m_jvt_DOWN->Fill(bornMass_KFactor,weight_jvt_DOWN);


	      p_pt_Tot->Fill(pt,finalWeight);											     
	      p_pt_Rec->Fill(pt,weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);    	     

	      p_pt_KF ->Fill(pt,kFactor);
	      p_pt_MC ->Fill(pt,mc_weight);    

	      p_pt_Iso->Fill(pt,weight_indiv_SF_MU_Isol); 
	      p_pt_Iso_sys_UP->Fill(pt,weight_indiv_SF_MU_Isol_SYST_UP);
	      p_pt_Iso_sys_DOWN->Fill(pt,weight_indiv_SF_MU_Isol_SYST_DOWN);
	      p_pt_Iso_stat_UP->Fill(pt,weight_indiv_SF_MU_Isol_STAT_UP);
	      p_pt_Iso_stat_DOWN->Fill(pt,weight_indiv_SF_MU_Isol_STAT_DOWN);

	      p_pt_ID ->Fill(pt,weight_indiv_SF_MU_ID);
	      p_pt_ID_sys_UP ->Fill(pt,weight_indiv_SF_MU_ID_SYST_UP);
	      p_pt_ID_sys_DOWN ->Fill(pt,weight_indiv_SF_MU_ID_SYST_DOWN);
	      p_pt_ID_stat_UP ->Fill(pt,weight_indiv_SF_MU_ID_STAT_UP);
	      p_pt_ID_stat_DOWN ->Fill(pt,weight_indiv_SF_MU_ID_STAT_DOWN);

	      p_pt_Tri->Fill(pt,weight_globalLeptonTriggerSF);
	      p_pt_Tri_sys_UP->Fill(pt,weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP);
	      p_pt_Tri_sys_DOWN->Fill(pt,weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN);
	      p_pt_Tri_stat_UP->Fill(pt,weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP);
	      p_pt_Tri_stat_DOWN->Fill(pt,weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN);

	      p_pt_ttv->Fill(pt,weight_indiv_SF_MU_TTVA);
	      p_pt_ttv_sys_UP->Fill(pt,weight_indiv_SF_MU_TTVA_SYST_UP);
	      p_pt_ttv_sys_DOWN->Fill(pt,weight_indiv_SF_MU_TTVA_SYST_DOWN);
	      p_pt_ttv_stat_UP->Fill(pt,weight_indiv_SF_MU_TTVA_STAT_UP);
	      p_pt_ttv_stat_DOWN->Fill(pt,weight_indiv_SF_MU_TTVA_STAT_DOWN);

	      p_pt_PU ->Fill(pt,pileup);
	      p_pt_PU_UP ->Fill(pt,weight_pileup_UP);
	      p_pt_PU_DOWN ->Fill(pt,weight_pileup_DOWN);

	      p_pt_jvt->Fill(pt,weight_jvt);    
	      p_pt_jvt_UP->Fill(pt,weight_jvt_UP);    
	      p_pt_jvt_DOWN->Fill(pt,weight_jvt_DOWN);    


	      p_eta_Tot->Fill(std::fabs(eta),finalWeight);											     
	      p_eta_Rec->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);    	     

	      p_eta_KF ->Fill(std::fabs(eta),kFactor);
	      p_eta_MC ->Fill(std::fabs(eta),mc_weight);    

	      p_eta_Iso->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol); 
	      p_eta_Iso_sys_UP->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol_SYST_UP);
	      p_eta_Iso_sys_DOWN->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol_SYST_DOWN);
	      p_eta_Iso_stat_UP->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol_STAT_UP);
	      p_eta_Iso_stat_DOWN->Fill(std::fabs(eta),weight_indiv_SF_MU_Isol_STAT_DOWN);

	      p_eta_ID ->Fill(std::fabs(eta),weight_indiv_SF_MU_ID);
	      p_eta_ID_sys_UP ->Fill(std::fabs(eta),weight_indiv_SF_MU_ID_SYST_UP);
	      p_eta_ID_sys_DOWN ->Fill(std::fabs(eta),weight_indiv_SF_MU_ID_SYST_DOWN);
	      p_eta_ID_stat_UP ->Fill(std::fabs(eta),weight_indiv_SF_MU_ID_STAT_UP);
	      p_eta_ID_stat_DOWN ->Fill(std::fabs(eta),weight_indiv_SF_MU_ID_STAT_DOWN);

	      p_eta_Tri->Fill(std::fabs(eta),weight_globalLeptonTriggerSF);
	      p_eta_Tri_sys_UP->Fill(std::fabs(eta),weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP);
	      p_eta_Tri_sys_DOWN->Fill(std::fabs(eta),weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN);
	      p_eta_Tri_stat_UP->Fill(std::fabs(eta),weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP);
	      p_eta_Tri_stat_DOWN->Fill(std::fabs(eta),weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN);

	      p_eta_ttv->Fill(std::fabs(eta),weight_indiv_SF_MU_TTVA);
	      p_eta_ttv_sys_UP->Fill(std::fabs(eta),weight_indiv_SF_MU_TTVA_SYST_UP);
	      p_eta_ttv_sys_DOWN->Fill(std::fabs(eta),weight_indiv_SF_MU_TTVA_SYST_DOWN);
	      p_eta_ttv_stat_UP->Fill(std::fabs(eta),weight_indiv_SF_MU_TTVA_STAT_UP);
	      p_eta_ttv_stat_DOWN->Fill(std::fabs(eta),weight_indiv_SF_MU_TTVA_STAT_DOWN);

	      p_eta_PU ->Fill(std::fabs(eta),pileup);
	      p_eta_PU_UP ->Fill(std::fabs(eta),weight_pileup_UP);
	      p_eta_PU_DOWN ->Fill(std::fabs(eta),weight_pileup_DOWN);

	      p_eta_jvt->Fill(std::fabs(eta),weight_jvt);    
	      p_eta_jvt_UP->Fill(std::fabs(eta),weight_jvt_UP);    
	      p_eta_jvt_DOWN->Fill(std::fabs(eta),weight_jvt_DOWN);    

	    }

	    if(goodZEvent && !DoMultijetAnalysis){	    	    
	      p_pt_Tot->Fill(pt2,finalWeight);
	      p_pt_Rec->Fill(pt2,weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);
	      p_pt_Iso->Fill(pt2,weight_indiv_SF_MU_Isol);	     
	      p_pt_ID ->Fill(pt2,weight_indiv_SF_MU_ID);
	      p_pt_Tri->Fill(pt2,weight_globalLeptonTriggerSF);
	      p_pt_ttv->Fill(pt2,weight_indiv_SF_MU_TTVA);
	      p_pt_KF ->Fill(pt2,kFactor);
	      p_pt_PU ->Fill(pt2,pileup);
	      p_pt_MC ->Fill(pt2,mc_weight);    
	      p_eta_Tot->Fill(std::fabs(eta2),finalWeight);
	      p_eta_Rec->Fill(std::fabs(eta2),weight_indiv_SF_MU_Isol*weight_globalLeptonTriggerSF*weight_indiv_SF_MU_ID*weight_indiv_SF_MU_TTVA);
	      p_eta_Iso->Fill(std::fabs(eta2),weight_indiv_SF_MU_Isol);
	      p_eta_ID ->Fill(std::fabs(eta2),weight_indiv_SF_MU_ID);		     
	      p_eta_Tri->Fill(std::fabs(eta2),weight_globalLeptonTriggerSF);
	      p_eta_ttv->Fill(std::fabs(eta2),weight_indiv_SF_MU_TTVA);
	      p_eta_KF ->Fill(std::fabs(eta2),kFactor);
	      p_eta_PU ->Fill(std::fabs(eta2),pileup);
	      p_eta_MC ->Fill(std::fabs(eta2),mc_weight);
	    }

	    recoLeptons.clear();

	    //multijet fill histos
	    if(DoMultijetAnalysis && passFitRegion){
	      for(int iBin = 0; iBin<(int)eta_bin_name.size(); iBin++){
		FillHistosMultijet(iBin,eta_bin_low[iBin],eta_bin_high[iBin]);
	      }
	    }

	  }//if basic Z and W	    
	  
      }//WZSelection
      

	//fill histograms with no cuts************
	if(!DoMultijetAnalysis){
        h_truth_mass_nocuts->Fill(1.,pileup*mc_weight);
        h_nm->Fill(mu_pt->size(),pileup*mc_weight);     
        hPileUp->Fill(weight_pileup);}
      
      nb = fChain->GetEntry(jentry);   nbytes += nb; 
      // if (Cut(ientry) < 0) continue;   
    }//end reco loop
     

    //CUTFLOWS
    TString chanflowdir, cutflowdir, cutflowdir2;
    if(BosonZSelection)
      chanflowdir="mumu";
    if(BosonWplusSelection || BosonWminusSelection)
      chanflowdir="munu";    
    if(DataYears2015p2016)
      cutflowdir=chanflowdir+"_2016/cutflow_mc_pu";
    else if(DataYears2017)
      cutflowdir=chanflowdir+"_2017/cutflow_mc_pu";
    else if(DataYears2018)
      cutflowdir=chanflowdir+"_2018/cutflow_mc_pu";   
    
    //obtain the directories
    TSystemDirectory directory("",(config.InputFileDir+nameOfSample+"/").c_str());
    TList *list = directory.GetListOfFiles();  
    std::vector<TString> names;
    
    if(list){
      TSystemFile *file; TString fname; TIter next(list);    
      while((file=(TSystemFile*)next())){
	fname = file->GetName();
	if (!file->IsDirectory() && fname.EndsWith(".root"))
	  names.push_back(fname);
      }
      delete list;
    }
    
    TH1F *hTot; TH1F *hTot2; TFile *fIN; TFile *fIN2;
    int nFiles = (int)names.size();

    //get the cutflow histograms
    if(nFiles>1){
      for(int k=0; k<nFiles; k++){	
	fIN = new TFile((config.InputFileDir+nameOfSample+"/").c_str()+names[k]);
	TH1F *dummy;
	if(!fIN->IsOpen()){ std::cout<< "Error: something wrong with cutflow!!!"<<std::endl; break;}      
	if(k==0){
	  hTot = (TH1F*)fIN->Get(cutflowdir);
	}else if(k!=0){
	  dummy = (TH1F*)fIN->Get(cutflowdir);
	  hTot->Add(dummy);
	}               
      }//endloop
    }else if(nFiles==1){//ugly fix
      fIN2 = new TFile((config.InputFileDir+nameOfSample+"/").c_str()+names[0]);
      hTot2 = (TH1F*)fIN2->Get(cutflowdir);
    }

    //fill a new histo for cutflow
    unsigned int nx=0;
    if(nFiles==1){
      nx = hTot2->GetXaxis()->GetNbins();
      for(int k=1; k<nx+1;k++){
	h_cutflow     ->SetBinContent(k,hTot2->GetBinContent(k));
	h_cutflow_high->SetBinContent(k,hTot2->GetBinContent(k));	
      }
    }else{
      nx = hTot->GetXaxis()->GetNbins();
      for(int k=1; k<nx+1;k++){
	h_cutflow     ->SetBinContent(k,hTot->GetBinContent(k));
	h_cutflow_high->SetBinContent(k,hTot->GetBinContent(k));
      }
    }
      
    //Inclusive cutflow
    h_cutflow->SetBinContent(nx+1,cutFlowIncl1);
    h_cutflow->SetBinContent(nx+2,cutFlowIncl2);
    h_cutflow->SetBinContent(nx+3,cutFlowIncl3);
    h_cutflow->SetBinContent(nx+4,cutFlowIncl4);
    h_cutflow->SetBinContent(nx+5,cutFlowIncl5);
    h_cutflow->SetBinContent(nx+6,cutFlowIncl6);
    //High cutflow
    h_cutflow_high->SetBinContent(nx+1,cutFlowHigh1);
    h_cutflow_high->SetBinContent(nx+2,cutFlowHigh2);
    h_cutflow_high->SetBinContent(nx+3,cutFlowHigh3);
    h_cutflow_high->SetBinContent(nx+4,cutFlowHigh4);
    h_cutflow_high->SetBinContent(nx+5,cutFlowHigh5);
    h_cutflow_high->SetBinContent(nx+6,cutFlowHigh6);

    //normalise cutflow
    double cutFlowNorm=1.;
    if(!isData) cutFlowNorm=(lumi_sample/tt2p)*(weightsSum_7/weightsSum_8);
    h_cutflow     ->Scale(cutFlowNorm);
    h_cutflow_high->Scale(cutFlowNorm);

    //display cutflows
    std::cout<<"     Inclusive CutFlow        "<<std::endl;
    std::cout<<cutFlowIncl1*cutFlowNorm<<"    "<<cutFlowNorm<<std::endl;
    std::cout<<cutFlowIncl2*cutFlowNorm<<std::endl;
    std::cout<<cutFlowIncl3*cutFlowNorm<<std::endl;
    std::cout<<cutFlowIncl4*cutFlowNorm<<std::endl;
    std::cout<<cutFlowIncl5*cutFlowNorm<<std::endl;
    std::cout<<cutFlowIncl6*cutFlowNorm<<std::endl;
    // std::cout<<cutFlowIncl7<<std::endl;
    // std::cout<<cutFlowIncl8<<std::endl;

    std::cout<<"     HighMass CutFlow      "<<std::endl;
    std::cout<<cutFlowHigh1*cutFlowNorm<<"    "<<cutFlowNorm<<std::endl; 
    std::cout<<cutFlowHigh2*cutFlowNorm<<std::endl;
    std::cout<<cutFlowHigh3*cutFlowNorm<<std::endl;
    std::cout<<cutFlowHigh4*cutFlowNorm<<std::endl;
    std::cout<<cutFlowHigh5*cutFlowNorm<<std::endl; 
    std::cout<<cutFlowHigh6*cutFlowNorm<<std::endl;
    // std::cout<<cutFlowHigh7<<std::endl;
    // std::cout<<cutFlowHigh8<<std::endl;

    //normalise and write histograms!
    double xSecVar=1.;
    if(config.MultijetVariation=="xSecdown" && config.DoMultijet=="True") xSecVar=0.95;
    if(config.MultijetVariation=="xSecup"   && config.DoMultijet=="True") xSecVar=1.05;

    double d0sig_norm = 1.;
    if(Applyd0WeightAnalysis && !isData )
      d0sig_norm=weightsSum_9/weightsSum_10;// pending
    d0sig_norm=d0sig_norm*xSecVar;

    double pileup_norm=1;
    pileup_norm=weightsSum_7/weightsSum_8;

    double trig_norm=1;
    //trig_norm=weightsSum_1/weightsSum_2;

    std::cout<<lumi_sample*pileup_norm*d0sig_norm<<" precioso parrito"<<lumi_sample<<"   "<<pileup_norm<<"   "<<d0sig_norm<<std::endl;
    if(!DoMultijetAnalysis){
      for(int iSF_variations = 0; iSF_variations < SFloopSize; ++iSF_variations){
	LumiNormalisation(config, iSF_variations, isData, lumi_sample*pileup_norm*d0sig_norm*trig_norm, tt2p);
	WriteHist(iSF_variations);
      }
    }else if(DoMultijetAnalysis){
      for(int iBin = 0; iBin<(int)eta_bin_name.size();iBin++){
	LumiNormalisation(config,iBin,isData,lumi_sample*pileup_norm*d0sig_norm*trig_norm,tt2p);
	WriteHistMultijet(iBin);
      }
    }
    
    
    std::cout<<weightsSum_1<<"    "<<weightsSum_2<<"  "<<weightsSum_3<<"    "<<weightsSum_4<<"    "<<weightsSum_5<<"    "<<weightsSum_6<<"    "<<weightsSum_7<<"   "<<weightsSum_8<<std::endl;                  
    std::cout<<"    "<<weightsSumWZ<<"    "<<lumi_sample<<"    "<<tt2p<<"       Normalisation  Reco"<<std::endl;
    
    std::cout<<"   Factores de escala:     "<<lumi_sample<<"    "<<skimNorm<<"     "<<std::endl;
    std::cout<<"   Size of the map:  "<<m_myTruthMap.size()<<"      Size of Truth Tree:    "<<nTruth<<"    MapCounter:   "<<mapCounter<<"   Size of Reco Tree:    "<<"   "<<nReco<<std::endl;
    
    if(HasRecoilInfoAnalysis) m_recoilHistos->finalize();
    
    std::cout<<" Trigger scale factor norm   "<<weightsSum_1/weightsSum_2<<"     "<<weightsSum_1<<"       "<<weightsSum_2<<std::endl;
    std::cout<<" pile up Trigger scale factor norm   "<<weightsSum_7/weightsSum_8<<"     "<<weightsSum_7<<"       "<<weightsSum_8<<std::endl;
    
}//function MyWZAnalysis::Loop


void MyWZAnalysis::Getd0sigFit(Config config){
  
  std::string dirYear="", dataname="";
  if(DataYears2015p2016) {dirYear="2015p2016/"; dataname="data15p16_";}
  if(DataYears2017) {dirYear="2017/"; dataname="data17_";}
  if(DataYears2018) {dirYear="2018/"; dataname="data18_";}
  
  TFile *fdata       = new TFile((config.OutputFileDir + "Files/"+dirYear+"Multijet/"+dataname+config.WZSelection+"_multijet_signal.root").c_str());
  TFile *fmc_ori     = new TFile((config.OutputFileDir + "Files/"+dirYear+"Multijet/background_"+config.WZSelection+"_nominal_nonewPU_multijet_signal_original.root").c_str());
  TFile *fmc_shifted = new TFile((config.OutputFileDir + "Files/"+dirYear+"Multijet/background_"+config.WZSelection+"_nominal_nonewPU_multijet_signal_shifted.root").c_str());

  //if(fdata==NULL || fmc_ori==NULL || fmc_shifted==NULL) return;
  if( !fdata->IsOpen() || !fmc_ori->IsOpen() || !fmc_shifted->IsOpen()){std::cout<<"No d0sig shift files"<<std::endl; return;}

  TH1D *hData, *hBG_ori, *hBG;
  hData   = (TH1D*)fdata       ->GetObjectUnchecked("h_d0sig_0");//i is for d0sig<3, 0 is for no d0sig cut, 1 is for d0sig > 3
  hBG_ori = (TH1D*)fmc_ori     ->GetObjectUnchecked("h_d0sig_0");
  
  hBG     = (TH1D*)fmc_shifted ->GetObjectUnchecked("h_d0sig_0");
   
  std::cout<<hData->GetMean()<<"   data mean (average value)"<<std::endl;
  std::cout<<hBG  ->GetMean()<<"   MC   mean (average value) shifted "<<std::endl;
  std::cout<<hData->GetMean()-hBG_ori->GetMean()<<"   difference data - MCoriginal mean (average value)!!SHIFT"<<std::endl;

  shift= hData->GetMean()-hBG_ori->GetMean();
  
  double limit=0;
  if(config.MultijetVariation=="Nominal" || config.MultijetVariation=="xSecup"|| config.MultijetVariation=="xSecdown") limit = -2.5;
  else if(config.MultijetVariation=="20") limit = -2.0;
  else if(config.MultijetVariation=="30") limit = -3.0;

  hData->Fit("gaus","","",limit,std::fabs(limit));
  TF1 *fit_data = hData->GetFunction("gaus");
  TF1 *fit_bg;
  hBG  ->Fit("gaus","","",limit,std::fabs(limit));
  fit_bg  = hBG  ->GetFunction("gaus");

  TF1 *ExtraPolate;
  ExtraPolate = new TF1("ExtraPolate","([0]/[3])*exp(-0.5*((x-[1])/[2])**2)/exp(-0.5*((x-[4])/[5])**2)",-10,10);
  ExtraPolate->SetParameter(0, fit_data->GetParameter(0));
  ExtraPolate->SetParameter(1, fit_data->GetParameter(1));
  ExtraPolate->SetParameter(2, fit_data->GetParameter(2));
  ExtraPolate->SetParameter(3, fit_bg  ->GetParameter(0));
  ExtraPolate->SetParameter(4, fit_bg  ->GetParameter(1));
  ExtraPolate->SetParameter(5, fit_bg  ->GetParameter(2));
  
  fit_d0sig_gaussratio = ExtraPolate;
  
  fdata       ->Close(); 
  fmc_ori     ->Close();
  fmc_shifted ->Close();
    
  return;
}


void MyWZAnalysis::setHist(Config config, std::string sf_name){

  std::cout<<"  Histogram set for:  "<<sf_name<<"  systematic"<<std::endl;
  Int_t nbins = 7;
  Double_t x0 = TMath::Log10(50);
  Double_t x1 = TMath::Log10(100);
  Double_t x2 = 0.3;
  
  //Double_t bins[]={x0,x1,x1+1*(0.25),x1+2*(0.25),x1+3*(0.25),3.};
  double bins[]={x0,x1,x1+1*(0.25),x1+2*(0.25),x1+3*(0.25),3.,3.+x2,3.+2*x2};
  //Double_t bins[]={x0,2.5,2.5+x2,2.5+2*x2,3.};
  
  //****** Log binning *******
  float Ptmin = 55000.; float Ptmax = 3000000.;
  Ptmin = log10(Ptmin); Ptmax = log10(Ptmax);
  const int nPtEdge1 = 51;
  float PtBins1[nPtEdge1];
  float ptWidth = (Ptmax - Ptmin)/(nPtEdge1-1);
  for (int i = 0; i < nPtEdge1; i++){
    PtBins1[i  ] = Ptmin+(i  )*ptWidth;
    PtBins1[i  ] = pow(10,PtBins1[i  ]);
    //printf("bin  i=%3i  edge= %9.4f   width=%9.2f\n",i,PtBins1[i],ptWidth);
  }
  float Massmin = 110000.; float Massmax = 7000000.;
  Massmin = log10(Massmin); Massmax = log10(Massmax);
  const int nMassEdge = 51;
  float MassBins[nMassEdge];
  float binWidth = (Massmax - Massmin)/(nMassEdge-1);
  for (int i = 0; i < nMassEdge; i++){
    MassBins[i  ] = Massmin+(i  )*binWidth;
    MassBins[i  ] = pow(10,MassBins[i  ]);
    //printf("bin  i=%3i  lo=%9.4f  hi=%9.4f \n",i,MassBins[i-1],MassBins[i]);
  }
  
  float MassminT = 50000.; float MassmaxT = 7000000.;
  MassminT = log10(MassminT); MassmaxT = log10(MassmaxT);
  const int nMassEdgeT = 101;
  float MassBinsT[nMassEdgeT];
  float binWidthT = (MassmaxT - MassminT)/(nMassEdgeT-1);
  for (int i = 0; i < nMassEdgeT; i++){
    MassBinsT[i  ] = MassminT+(i  )*binWidthT;
    MassBinsT[i  ] = pow(10,MassBinsT[i  ]);
  }
      
  PtRmass.push_back( new TH2D(("PtRmass"+sf_name).c_str(),";   M_{W,T}  [GeV];   p_{T}^{#mu}",nbins,bins,nbins,bins) );	
  PtRmass_high.push_back( new TH2D(("PtRmass_high"+sf_name).c_str(),";   M_{W,T}  [GeV];   p_{T}^{#mu}",nbins,bins,nbins,bins) );
  
  h_rmassLogCw.push_back( new TH1D(("h_rmassLogCw"+sf_name).c_str(),"Cw;   M_{W,T}  [GeV]; Events",nbins,bins) );
  h_rmassLogCw_high.push_back( new TH1D(("h_rmassLogCw_high"+sf_name).c_str(),"Cw;   M_{W,T}  [GeV]; Events",nbins,bins) );
  
  h_mu.push_back( new TH1D(("h_mu"+sf_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_high.push_back( new TH1D(("h_mu_high"+sf_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  
  h_npvtx.push_back( new TH1I(("h_npvtx"+sf_name).c_str(),"nPVTX;  nPVTX;  Events",20,0,40) );
  h_npvtx_high.push_back( new TH1I(("h_npvtx_high"+sf_name).c_str(),"nPVTX;  nPVTX;  Events",20,0,40) );
  
  h_pt20Track.push_back( new TH1D(("h_pt20Track"+sf_name).c_str(),"Track Isolation    ;Track20^{#mu}; Events",100,0.0,0.5) );
  h_pt20Track_high.push_back( new TH1D(("h_pt20Track_high"+sf_name).c_str(),"Track Isolation    ;Track20^{#mu}; Events",100,0.0,0.5) );
  
  h_SumET.push_back( new TH1D(("h_SumET"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  h_SumET_high.push_back( new TH1D(("h_SumET_high"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  
  h_SumET_eta1.push_back( new TH1D(("h_SumET_eta1"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  h_SumET_eta2.push_back( new TH1D(("h_SumET_eta2"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  h_SumET_eta3.push_back( new TH1D(("h_SumET_eta3"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  h_SumET_eta4.push_back( new TH1D(("h_SumET_eta4"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  h_SumET_eta5.push_back( new TH1D(("h_SumET_eta5"+sf_name).c_str(),"SumET    ;#sum ET; Events",100,0.,2000000.) );
  
  h_pt20Calo.push_back( new TH1D(("h_pt20Calo"+sf_name).c_str(),"Calo Isolation    ;Calo20^{#mu}; Events",100,0.0,0.5) );
  h_pt20Calo_high.push_back( new TH1D(("h_pt20Calo_high"+sf_name).c_str(),"Calo Isolation    ;Calo20^{#mu}; Events",100,0.0,0.5) );
  
  h_pt.push_back( new TH1D(("h_pt"+sf_name).c_str(),"Muon Transverse Momentum    ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );//,50,0,120000.0
  h_pt_high.push_back( new TH1D(("h_pt_high"+sf_name).c_str(), "muon pt        ;p_{T}^{#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );//105 bins to nice slices
  
  double eta_bins[]={0.0,0.21,0.42,0.63,0.84,1.05,1.37,1.52,1.74,1.95,2.18,2.4};//8,7 TeV binning
  //double eta_bins[]={0.0,0.21,0.42,0.64,0.86,1.08,1.30,1.52,1.74,1.96,2.18,2.5}; //andres
  
  h_eta.push_back( new TH1D(("h_eta"+sf_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );//(52,-2.6,2.6);//13bins default
  h_eta_high.push_back( new TH1D(("h_eta_high"+sf_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  
  h_eta_d.push_back( new TH1D(("h_eta_d"+sf_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );//13bins default
  h_eta_high_d.push_back( new TH1D(("h_eta_d_high"+sf_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  
  h_phi.push_back( new TH1D(("h_phi"+sf_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_high.push_back( new TH1D(("h_phi_high"+sf_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  
  h_phi_met.push_back( new TH1D(("h_phi_met"+sf_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_high.push_back( new TH1D(("h_phi_met_high"+sf_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  
  h_delta_phi.push_back( new TH1D(("h_delta_phi"+sf_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_high.push_back( new TH1D(("h_delta_phi_high"+sf_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  
  h_met.push_back( new TH1D(("h_met"+sf_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_high.push_back( new TH1D(("h_met_high"+sf_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );
  
  h_met_recoil.push_back( new TH1D(("h_met_recoil"+sf_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_recoil_high.push_back( new TH1D(("h_met_recoil_high"+sf_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );
  
  h_rmass.push_back( new TH1D(("h_rmass"+sf_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_high.push_back( new TH1D(("h_rmass_high"+sf_name).c_str(),"mt high;  M_{T}^{#mu}  [MeV]; Events",nMassEdge-1,MassBins) );//55,40000.,1100000.
  
  h_rmass_recoil.push_back( new TH1D(("h_rmass_recoil"+sf_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",55,40000.,150000.) );
  h_rmass_recoil_high.push_back( new TH1D(("h_rmass_recoil_high"+sf_name).c_str(),"mt high;  M_{T}^{#mu}  [MeV]; Events",nMassEdge-1,MassBins) );//55,40000.,1100000.
    
  h_d0.push_back( new TH1D(("h_d0"+sf_name).c_str(),"MET Phi;   d0; Events",80,-10.,10.) );
  h_d0_high.push_back( new TH1D(("h_d0_high"+sf_name).c_str(),"MET Phi;   d0; Events",80,-10.,10.) );
  
  h_z0.push_back( new TH1D(("h_z0"+sf_name).c_str(),"MET Phi; z0;  Events",20,-.5,.5) );
  h_z0_high.push_back( new TH1D(("h_z0_high"+sf_name).c_str(),"MET Phi; z0;  Events",20,-.5,.5) );
   
  h_pt_boson_reco.push_back( new TH1D(("h_pt_boson_reco"+sf_name).c_str(),"Boson Transverse Momentum    ;p_{T}_{#mu#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );//100,0.,200000.);
  h_pt_boson_reco_p.push_back( new TH1D(("h_pt_boson_reco_p"+sf_name).c_str(),"Boson Transverse Momentum    ;p_{T}_{#mu#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );
  
  //h boson control plots
  h_pt_boson.push_back( new TH1D(("h_pt_boson"+sf_name).c_str(),"Boson Transverse Momentum    ;p_{T}_{#mu#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );
  
  //double GeV=1000.;
  h_mass_boson.push_back( new TH1D(("h_mass_boson"+sf_name).c_str(),"Boson Truth Mass ;m  [MeV]; Events",nMassEdgeT-1,MassBinsT) );
  h_eta_boson.push_back( new TH1D(("h_eta_boson"+sf_name).c_str(),"Boson Pseudorapidity    ;#eta^{#mu#mu}; Events",52,-2.6,2.6) );//13bins default
  
  h_pt_boson_high.push_back( new TH1D(("h_pt_boson_high"+sf_name).c_str(),"Boson Transverse Momentum    ;p_{T}_{#mu#mu}  [MeV]; Events",nPtEdge1-1, PtBins1) );
  
  h_mass_boson_high.push_back( new TH1D(("h_mass_boson_high"+sf_name).c_str(),"Boson Invariant mass ;m_{#mu#mu}  [MeV]; Events",nMassEdgeT-1,MassBinsT) );
  h_eta_boson_high.push_back( new TH1D(("h_eta_boson_high"+sf_name).c_str(),"Boson Pseudorapidity ;#eta^{#mu#mu}; Events",52,-2.6,2.6) );//13bins def 12,eta_bins);
  

  //simple histos with no cuts at all
  h_truth_mass_nocuts = new TH1D(("h_truth_mass_nocuts"+sf_name).c_str(),"mt high;  M_{T}^{#mu}  [MeV]; Events",nMassEdge-1,MassBins);
  h_nm = new TH1I(("h_nm"+sf_name).c_str(),"# of muons",10,0,10);
  hPileUp  = new  TH1D(("hPileup"+sf_name).c_str(), "; Pileup weight; Entries ", 100,0,5.);

  if(sf_name=="nominal"){

    //cutflow
    h_cutflow      = new TH1D(("h_cutflow"+sf_name).c_str(),      "Cutflow; Events", 25,0,25);
    h_cutflow_high = new TH1D(("h_cutflow_high"+sf_name).c_str(), "Cutflow; Events", 25,0,25);
       
    //TProfiles
    p_mwt_Tot = new TProfile(("p_mwt_Tot"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Rec = new TProfile(("p_mwt_Rec"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_KF  = new TProfile(("p_mwt_KF"+sf_name).c_str(),  ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_MC  = new TProfile(("p_mwt_MC"+sf_name).c_str(),  ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_Iso           = new TProfile(("p_mwt_Iso"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Iso_sys_UP    = new TProfile(("p_mwt_Iso_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Iso_sys_DOWN  = new TProfile(("p_mwt_Iso_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Iso_stat_UP   = new TProfile(("p_mwt_Iso_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Iso_stat_DOWN = new TProfile(("p_mwt_Iso_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_ID            = new TProfile(("p_mwt_ID"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ID_sys_UP     = new TProfile(("p_mwt_ID_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ID_sys_DOWN   = new TProfile(("p_mwt_ID_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ID_stat_UP    = new TProfile(("p_mwt_ID_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ID_stat_DOWN  = new TProfile(("p_mwt_ID_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_Tri           = new TProfile(("p_mwt_Tri"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Tri_sys_UP    = new TProfile(("p_mwt_Tri_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Tri_sys_DOWN  = new TProfile(("p_mwt_Tri_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Tri_stat_UP   = new TProfile(("p_mwt_Tri_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_Tri_stat_DOWN = new TProfile(("p_mwt_Tri_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_ttv           = new TProfile(("p_mwt_ttva"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ttv_sys_UP    = new TProfile(("p_mwt_ttva_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ttv_sys_DOWN  = new TProfile(("p_mwt_ttva_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ttv_stat_UP   = new TProfile(("p_mwt_ttva_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_ttv_stat_DOWN = new TProfile(("p_mwt_ttva_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_PU           = new TProfile(("p_mwt_PU"     +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_PU_UP        = new TProfile(("p_mwt_PU_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_PU_DOWN      = new TProfile(("p_mwt_PU_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_mwt_jvt          = new TProfile(("p_mwt_jvt_PU"     +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_mwt_jvt_UP       = new TProfile(("p_mwt_jvt_PU_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT); 
    p_mwt_jvt_DOWN     = new TProfile(("p_mwt_jvt_PU_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);   
    
    p_m_Tot = new TProfile(("p_m_Tot"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Rec = new TProfile(("p_m_Rec"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_KF  = new TProfile(("p_m_KF"+sf_name).c_str(),  ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_MC  = new TProfile(("p_m_MC"+sf_name).c_str(),  ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_Iso           = new TProfile(("p_m_Iso"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Iso_sys_UP    = new TProfile(("p_m_Iso_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Iso_sys_DOWN  = new TProfile(("p_m_Iso_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Iso_stat_UP   = new TProfile(("p_m_Iso_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Iso_stat_DOWN = new TProfile(("p_m_Iso_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_ID            = new TProfile(("p_m_ID"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ID_sys_UP     = new TProfile(("p_m_ID_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ID_sys_DOWN   = new TProfile(("p_m_ID_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ID_stat_UP    = new TProfile(("p_m_ID_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ID_stat_DOWN  = new TProfile(("p_m_ID_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_Tri           = new TProfile(("p_m_Tri"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Tri_sys_UP    = new TProfile(("p_m_Tri_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Tri_sys_DOWN  = new TProfile(("p_m_Tri_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Tri_stat_UP   = new TProfile(("p_m_Tri_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_Tri_stat_DOWN = new TProfile(("p_m_Tri_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_ttv           = new TProfile(("p_m_ttva"          +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ttv_sys_UP    = new TProfile(("p_m_ttva_sys_UP"   +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ttv_sys_DOWN  = new TProfile(("p_m_ttva_sys_DOWN" +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ttv_stat_UP   = new TProfile(("p_m_ttva_stat_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_ttv_stat_DOWN = new TProfile(("p_m_ttva_stat_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_PU           = new TProfile(("p_m_PU"     +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_PU_UP        = new TProfile(("p_m_PU_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_PU_DOWN      = new TProfile(("p_m_PU_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    
    p_m_jvt          = new TProfile(("p_m_jvt_PU"     +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);
    p_m_jvt_UP       = new TProfile(("p_m_jvt_PU_UP"  +sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT); 
    p_m_jvt_DOWN     = new TProfile(("p_m_jvt_PU_DOWN"+sf_name).c_str(), ";m_{T} [GeV]; Scale Factor",nMassEdgeT-1,MassBinsT);   
    
    
    p_pt_Tot = new TProfile(("p_pt_Tot"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Rec = new TProfile(("p_pt_Rec"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_KF  = new TProfile(("p_pt_KF"+sf_name).c_str(),  ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_MC  = new TProfile(("p_pt_MC"+sf_name).c_str(),  ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_Iso           = new TProfile(("p_pt_Iso"          +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Iso_sys_UP    = new TProfile(("p_pt_Iso_sys_UP"   +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Iso_sys_DOWN  = new TProfile(("p_pt_Iso_sys_DOWN" +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Iso_stat_UP   = new TProfile(("p_pt_Iso_stat_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Iso_stat_DOWN = new TProfile(("p_pt_Iso_stat_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_ID            = new TProfile(("p_pt_ID"          +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ID_sys_UP     = new TProfile(("p_pt_ID_sys_UP"   +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ID_sys_DOWN   = new TProfile(("p_pt_ID_sys_DOWN" +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ID_stat_UP    = new TProfile(("p_pt_ID_stat_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ID_stat_DOWN  = new TProfile(("p_pt_ID_stat_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_Tri           = new TProfile(("p_pt_Tri"          +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Tri_sys_UP    = new TProfile(("p_pt_Tri_sys_UP"   +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Tri_sys_DOWN  = new TProfile(("p_pt_Tri_sys_DOWN" +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Tri_stat_UP   = new TProfile(("p_pt_Tri_stat_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_Tri_stat_DOWN = new TProfile(("p_pt_Tri_stat_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_ttv           = new TProfile(("p_pt_ttva"          +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ttv_sys_UP    = new TProfile(("p_pt_ttva_sys_UP"   +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ttv_sys_DOWN  = new TProfile(("p_pt_ttva_sys_DOWN" +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ttv_stat_UP   = new TProfile(("p_pt_ttva_stat_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_ttv_stat_DOWN = new TProfile(("p_pt_ttva_stat_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_PU           = new TProfile(("p_pt_PU"     +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_PU_UP        = new TProfile(("p_pt_PU_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_PU_DOWN      = new TProfile(("p_pt_PU_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    
    p_pt_jvt          = new TProfile(("p_pt_jvt_PU"     +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);
    p_pt_jvt_UP       = new TProfile(("p_pt_jvt_PU_UP"  +sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1); 
    p_pt_jvt_DOWN     = new TProfile(("p_pt_jvt_PU_DOWN"+sf_name).c_str(), ";p_{T} [GeV]; Scale Factor",nPtEdge1-1, PtBins1);   
    
    p_eta_Tot = new TProfile(("p_eta_Tot"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Rec = new TProfile(("p_eta_Rec"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_KF  = new TProfile(("p_eta_KF"+sf_name).c_str(),  "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_MC  = new TProfile(("p_eta_MC"+sf_name).c_str(),  "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_Iso           = new TProfile(("p_eta_Iso"          +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Iso_sys_UP    = new TProfile(("p_eta_Iso_sys_UP"   +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Iso_sys_DOWN  = new TProfile(("p_eta_Iso_sys_DOWN" +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Iso_stat_UP   = new TProfile(("p_eta_Iso_stat_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Iso_stat_DOWN = new TProfile(("p_eta_Iso_stat_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_ID            = new TProfile(("p_eta_ID"          +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ID_sys_UP     = new TProfile(("p_eta_ID_sys_UP"   +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ID_sys_DOWN   = new TProfile(("p_eta_ID_sys_DOWN" +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ID_stat_UP    = new TProfile(("p_eta_ID_stat_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ID_stat_DOWN  = new TProfile(("p_eta_ID_stat_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_Tri           = new TProfile(("p_eta_Tri"          +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Tri_sys_UP    = new TProfile(("p_eta_Tri_sys_UP"   +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Tri_sys_DOWN  = new TProfile(("p_eta_Tri_sys_DOWN" +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Tri_stat_UP   = new TProfile(("p_eta_Tri_stat_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_Tri_stat_DOWN = new TProfile(("p_eta_Tri_stat_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_ttv           = new TProfile(("p_eta_ttva"          +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ttv_sys_UP    = new TProfile(("p_eta_ttva_sys_UP"   +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ttv_sys_DOWN  = new TProfile(("p_eta_ttva_sys_DOWN" +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ttv_stat_UP   = new TProfile(("p_eta_ttva_stat_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_ttv_stat_DOWN = new TProfile(("p_eta_ttva_stat_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_PU           = new TProfile(("p_eta_PU"     +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_PU_UP        = new TProfile(("p_eta_PU_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_PU_DOWN      = new TProfile(("p_eta_PU_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    
    p_eta_jvt          = new TProfile(("p_eta_jvt_PU"     +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
    p_eta_jvt_UP       = new TProfile(("p_eta_jvt_PU_UP"  +sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins); 
    p_eta_jvt_DOWN     = new TProfile(("p_eta_jvt_PU_DOWN"+sf_name).c_str(), "; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);   
  }

  //Cw factor calculation
  unsigned int nxbins = config.xBinsCw.size()-1;
  double xbins[nxbins+1];
  for(unsigned int i = 0; i<nxbins+1;i++) xbins[i]=config.xBinsCw.at(i);
  
  hMassCrec   = new  TH1D(("hMassCrec"+sf_name).c_str()  , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  hMassCgen   = new  TH1D(("hMassCgen"+sf_name).c_str()  , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  hMassCstay  = new  TH1D(("hMassCstay"+sf_name).c_str() , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  hMassCcome  = new  TH1D(("hMassCcome"+sf_name).c_str() , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  
  hEtaCrec   = new  TH1D(("hEtaCrec"+sf_name).c_str()  , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  hEtaCgen   = new  TH1D(("hEtaCgen"+sf_name).c_str()  , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  hEtaCstay  = new  TH1D(("hEtaCstay"+sf_name).c_str() , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  hEtaCcome  = new  TH1D(("hEtaCcome"+sf_name).c_str() , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  
  // TH1D* hMassCrec   = new  TH1D("hMassCrec"  , "mass; m_{T}^{reco} [GeV]; Entries ", 10, 40000.,150000.);
  // TH1D* hMassCgen   = new  TH1D("hMassCgen"  , "mass; m_{T}^{reco} [GeV]; Entries ", 10, 40000.,150000.);
  // TH1D* hMassCstay  = new  TH1D("hMassCstay" , "mass; m_{T}^{reco} [GeV]; Entries ", 10, 40000.,150000.);
  // TH1D* hMassCcome  = new  TH1D("hMassCcome" , "mass; m_{T}^{reco} [GeV]; Entries ", 10, 40000.,150000.);
  
  return;
}

void MyWZAnalysis::setHistMultijet(Config config, std::string bin_name){

  h_rmass_i.push_back( new TH1D(("h_rmass_i"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_0.push_back( new TH1D(("h_rmass_0"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_1.push_back( new TH1D(("h_rmass_1"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_2.push_back( new TH1D(("h_rmass_2"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_3.push_back( new TH1D(("h_rmass_3"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_4.push_back( new TH1D(("h_rmass_4"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  h_rmass_5.push_back( new TH1D(("h_rmass_5"+bin_name).c_str(),"Boson Reco Mass;  M_{T}^{#mu}  [MeV]; Events",60,0.,150000.) );
  
  h_pt_i.push_back( new TH1D(("h_pt_i"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_0.push_back( new TH1D(("h_pt_0"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_1.push_back( new TH1D(("h_pt_1"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_2.push_back( new TH1D(("h_pt_2"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_3.push_back( new TH1D(("h_pt_3"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_4.push_back( new TH1D(("h_pt_4"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  h_pt_5.push_back( new TH1D(("h_pt_5"+bin_name).c_str(),"Muon Transverse Momentum  ;p_{T}^{#mu}  [MeV]; Events",100,20000.,200000.) );
  
  double eta_bins[]={0.0,0.21,0.42,0.63,0.84,1.05,1.37,1.52,1.74,1.95,2.18,2.4};//8,7 TeV binning
  h_eta_i.push_back( new TH1D(("h_eta_i"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_0.push_back( new TH1D(("h_eta_0"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_1.push_back( new TH1D(("h_eta_1"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_2.push_back( new TH1D(("h_eta_2"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_3.push_back( new TH1D(("h_eta_3"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_4.push_back( new TH1D(("h_eta_4"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );
  h_eta_5.push_back( new TH1D(("h_eta_5"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins) );

  //display eta
  h_eta_d_i.push_back( new TH1D(("h_eta_d_i"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_0.push_back( new TH1D(("h_eta_d_0"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_1.push_back( new TH1D(("h_eta_d_1"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_2.push_back( new TH1D(("h_eta_d_2"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_3.push_back( new TH1D(("h_eta_d_3"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_4.push_back( new TH1D(("h_eta_d_4"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  h_eta_d_5.push_back( new TH1D(("h_eta_d_5"+bin_name).c_str(),"Muon Pseudorapidity    ;#eta^{#mu}; Events",52,-2.6,2.6) );
  
  h_met_i.push_back( new TH1D(("h_met_i"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_0.push_back( new TH1D(("h_met_0"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_1.push_back( new TH1D(("h_met_1"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_2.push_back( new TH1D(("h_met_2"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_3.push_back( new TH1D(("h_met_3"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_4.push_back( new TH1D(("h_met_4"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  h_met_5.push_back( new TH1D(("h_met_5"+bin_name).c_str(),"Missing Transverse Energy;  E_{T}^{#mu}  [MeV]; Events",80,0,200000.) );
  
  h_iso30_i.push_back( new TH1D(("h_iso30_i"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_0.push_back( new TH1D(("h_iso30_0"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_1.push_back( new TH1D(("h_iso30_1"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_2.push_back( new TH1D(("h_iso30_2"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_3.push_back( new TH1D(("h_iso30_3"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_4.push_back( new TH1D(("h_iso30_4"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );
  h_iso30_5.push_back( new TH1D(("h_iso30_5"+bin_name).c_str(),"Isolation;  p_{T}var30/p_{T}^{#mu}; Events",80,0,1.) );

  h_d0sig_i.push_back( new TH1D(("h_d0sig_i"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_0.push_back( new TH1D(("h_d0sig_0"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_1.push_back( new TH1D(("h_d0sig_1"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_2.push_back( new TH1D(("h_d0sig_2"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_3.push_back( new TH1D(("h_d0sig_3"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_4.push_back( new TH1D(("h_d0sig_4"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );
  h_d0sig_5.push_back( new TH1D(("h_d0sig_5"+bin_name).c_str(),"Isolation;  d0sig; Events",80,-10.,10.) );

  h_ptDmt_i.push_back( new TH1D(("h_ptDmt_i"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_0.push_back( new TH1D(("h_ptDmt_0"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_1.push_back( new TH1D(("h_ptDmt_1"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_2.push_back( new TH1D(("h_ptDmt_2"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_3.push_back( new TH1D(("h_ptDmt_3"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_4.push_back( new TH1D(("h_ptDmt_4"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );
  h_ptDmt_5.push_back( new TH1D(("h_ptDmt_5"+bin_name).c_str(),"Isolation;  p_{T}^{#mu}/m_{T}; Events",40,0,3) );

  h_phi_i.push_back( new TH1D(("h_phi_i"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_0.push_back( new TH1D(("h_phi_0"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_1.push_back( new TH1D(("h_phi_1"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_2.push_back( new TH1D(("h_phi_2"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_3.push_back( new TH1D(("h_phi_3"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_4.push_back( new TH1D(("h_phi_4"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );
  h_phi_5.push_back( new TH1D(("h_phi_5"+bin_name).c_str(),"Muon Phi    ;#Phi ^{#mu} ;Events",50,-3.1416,3.1416) );

  h_phi_met_i.push_back( new TH1D(("h_phi_met_i"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_0.push_back( new TH1D(("h_phi_met_0"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_1.push_back( new TH1D(("h_phi_met_1"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_2.push_back( new TH1D(("h_phi_met_2"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_3.push_back( new TH1D(("h_phi_met_3"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_4.push_back( new TH1D(("h_phi_met_4"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );
  h_phi_met_5.push_back( new TH1D(("h_phi_met_5"+bin_name).c_str(),"MET Phi   ;#Phi ^{Met}; Events",50,-3.1416,3.1416) );

  h_delta_phi_i.push_back( new TH1D(("h_delta_phi_i"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_0.push_back( new TH1D(("h_delta_phi_0"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_1.push_back( new TH1D(("h_delta_phi_1"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_2.push_back( new TH1D(("h_delta_phi_2"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_3.push_back( new TH1D(("h_delta_phi_3"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_4.push_back( new TH1D(("h_delta_phi_4"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );
  h_delta_phi_5.push_back( new TH1D(("h_delta_phi_5"+bin_name).c_str(),"MET Phi;   #Detla ^{#Phi}; Events",50,0.,3.1416) );

  h_mu_i.push_back( new TH1D(("h_mu_i"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_0.push_back( new TH1D(("h_mu_0"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_1.push_back( new TH1D(("h_mu_1"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_2.push_back( new TH1D(("h_mu_2"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_3.push_back( new TH1D(("h_mu_3"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_4.push_back( new TH1D(("h_mu_4"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );
  h_mu_5.push_back( new TH1D(("h_mu_5"+bin_name).c_str(),"Muon Transverse Momentum;  <#mu>;  Events",100,0.,100.) );


  return;
}


void MyWZAnalysis::FillHistos(int iSF_variations, double var){

  
  //hadronic recoil histos
  // if(goodWEventRecoil && !DoMultijetAnalysis){
  //   v_recoil_met+=lep1;
  //   v_recoil_met= -1* v_recoil_met;
  //   met_recoil=v_recoil_met.Pt();
  //   double dPhiRec=std::fabs(v_recoil_met.Phi()-mu_phi->at(0));	      
  //   if(dPhiRec>3.14159265) dPhiRec=2.*3.14159265-dPhiRec;
  //   recoMass_recoil = std::sqrt(2.*((mu_pt->at(0))*met_recoil)*(1.-std::cos(dPhiRec))); 
  //   h_met_recoil->Fill(met_recoil,finalWeight);
  //   h_rmass_recoil->Fill(recoMass_recoil,finalWeight);
  // }
  
  // if(goodWEventRecoilHigh && !DoMultijetAnalysis){
  //   v_recoil_met+=lep1;
  //   v_recoil_met= -1* v_recoil_met;
  //   met_recoil=v_recoil_met.Pt();
  //   double dPhiRec=std::fabs(v_recoil_met.Phi()-mu_phi->at(0));	      
  //   if(dPhiRec>3.14159265) dPhiRec=2.*3.14159265-dPhiRec;
  //   recoMass_recoil = std::sqrt(2.*((mu_pt->at(0))*met_recoil)*(1.-std::cos(dPhiRec))); 
  //   h_met_recoil_high->Fill(met_recoil,finalWeight);
  //   h_rmass_recoil_high->Fill(recoMass_recoil,finalWeight);
  // }

  
  if(goodWEvent || goodZEvent){

    //filling histos: event variable
    h_mu[iSF_variations]->Fill(mu_scaled,finalWeight*var);
    h_npvtx[iSF_variations]->Fill(primaryVertices,finalWeight*var);
    h_phi_met[iSF_variations]->Fill(phimet,finalWeight*var);
    h_SumET[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    h_met[iSF_variations]->Fill(met,finalWeight*var);   
    h_rmass[iSF_variations]->Fill(recoMass,finalWeight*var);
    h_rmassLogCw[iSF_variations]->Fill(lrmass,finalWeight*var);
    PtRmass[iSF_variations]->Fill(lrmass,lpt,finalWeight*var);
   
    if(std::fabs(hadronic_eta)<=0.8                               ) h_SumET_eta1[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    if(std::fabs(hadronic_eta)>0.8 && std::fabs(hadronic_eta)<=1.4) h_SumET_eta2[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    if(std::fabs(hadronic_eta)>1.4 && std::fabs(hadronic_eta)<=2.5) h_SumET_eta3[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    if(std::fabs(hadronic_eta)>2.5 && std::fabs(hadronic_eta)<=3.2) h_SumET_eta4[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    if(std::fabs(hadronic_eta)>3.2                                ) h_SumET_eta5[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
        
    //h_pt_boson_reco[iSF_variations]->Fill(pT_tot,finalWeight*var);
    if(goodWEvent) h_pt_boson_reco_p[iSF_variations]->Fill(Wboson.Pt(),finalWeight*var);

    //filling histos: muon variables
    h_pt[iSF_variations]->Fill(pt  ,finalWeight*var);
    h_eta[iSF_variations]->Fill(std::fabs(eta),finalWeight*var);
    h_eta_d[iSF_variations]->Fill(eta,finalWeight*var);
    h_phi[iSF_variations]->Fill(phi,finalWeight*var);
    h_pt20Calo[iSF_variations]->Fill(calo20/pt,finalWeight*var);	    
    h_pt20Track[iSF_variations]->Fill(track20,finalWeight*var);
    h_delta_phi[iSF_variations]->Fill(dPhi,finalWeight*var);
    h_z0[iSF_variations]->Fill(z0,finalWeight*var);
    h_d0[iSF_variations]->Fill(d0sig,finalWeight*var);
    
    if(goodZEvent){
      h_pt[iSF_variations]->Fill(pt2  ,finalWeight*var);
      h_eta[iSF_variations]->Fill(std::fabs(eta2),finalWeight*var);
      h_eta_d[iSF_variations]->Fill(eta2,finalWeight*var);
      h_phi[iSF_variations]->Fill(phi2,finalWeight*var);
      h_pt20Calo[iSF_variations]->Fill(calo202/pt2,finalWeight*var);	    
      h_pt20Track[iSF_variations]->Fill(track202,finalWeight*var);
      h_delta_phi[iSF_variations]->Fill(dPhi2,finalWeight*var);
      h_z0[iSF_variations]->Fill(z02,finalWeight*var);
      h_d0[iSF_variations]->Fill(d0sig2,finalWeight*var);
      h_pt_boson_reco_p[iSF_variations]->Fill(Zboson.Pt(),finalWeight*var);
    }
    
    h_pt_boson[iSF_variations]->Fill(tboson.Pt(),finalWeight*var);
    h_mass_boson[iSF_variations]->Fill(tboson.M(),finalWeight*var);
    h_eta_boson[iSF_variations]->Fill(tboson.Eta(),finalWeight*var);
  }//good Z or W event
  

  if(goodWEventHigh || goodZEventHigh){
    h_mu_high[iSF_variations]->Fill(mu_scaled,finalWeight*var);
    h_npvtx_high[iSF_variations]->Fill(primaryVertices,finalWeight*var);
    h_phi_met_high[iSF_variations]->Fill(phimet,finalWeight*var);
    h_SumET_high[iSF_variations]->Fill(sumET_PFO,finalWeight*var);
    h_met_high[iSF_variations]->Fill(met,finalWeight*var);
    h_rmass_high[iSF_variations]->Fill(recoMass,finalWeight*var);
    h_rmassLogCw_high[iSF_variations]->Fill(lrmass,finalWeight*var);
    PtRmass_high[iSF_variations]->Fill(lrmass,lpt,SF*pileup*kFactor*mc_weight);
    
    //filling histos: muon variables
    h_pt_high[iSF_variations]->Fill(pt  ,finalWeight*var);
    h_eta_high[iSF_variations]->Fill(std::fabs(eta),finalWeight*var);
    h_eta_high_d[iSF_variations]->Fill(eta,finalWeight*var);
    h_phi_high[iSF_variations]->Fill(phi,finalWeight*var);
    h_pt20Calo_high[iSF_variations]->Fill(calo20/pt,finalWeight*var);	    
    h_pt20Track_high[iSF_variations]->Fill(track20,finalWeight*var);
    h_delta_phi_high[iSF_variations]->Fill(dPhi,finalWeight*var);
    h_z0_high[iSF_variations]->Fill(z0,finalWeight*var);
    h_d0_high[iSF_variations]->Fill(d0sig,finalWeight*var);
    
    if(goodZEventHigh){
      h_pt_high[iSF_variations]->Fill(pt2  ,finalWeight*var);
      h_eta_high[iSF_variations]->Fill(std::abs(eta2),finalWeight*var);
      h_eta_high_d[iSF_variations]->Fill(eta2,finalWeight*var);
      h_phi_high[iSF_variations]->Fill(phi2,finalWeight*var);
      h_pt20Calo_high[iSF_variations]->Fill(calo202/pt2,finalWeight*var);	    
      h_pt20Track_high[iSF_variations]->Fill(track202,finalWeight*var);
      h_delta_phi_high[iSF_variations]->Fill(dPhi2,finalWeight*var);
      h_z0_high[iSF_variations]->Fill(z02,finalWeight*var);
      h_d0_high[iSF_variations]->Fill(d0sig2,finalWeight*var);
    }
    
    h_pt_boson_high[iSF_variations]->Fill(tboson.Pt(),finalWeight*var);
    h_mass_boson_high[iSF_variations]->Fill(tboson.M(),finalWeight*var);
    h_eta_boson_high[iSF_variations]->Fill(tboson.Eta(),finalWeight*var);
  }//high W or Z selection
 
  return; 
}


void MyWZAnalysis::FillHistosMultijet(int iEta_bins, double eta_low, double eta_high){


  if(std::fabs(eta) >= eta_low && std::fabs(eta)< eta_high || iEta_bins==0){

    //if(iso30<0.1){
    if(std::fabs(d0sig)<3.){
      h_rmass_i[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_i   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_i  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_met_i  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_i[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_i[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_i[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_i    [iEta_bins]->Fill(eta,finalWeight);
      h_phi_i      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_i  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_i[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_i       [iEta_bins]->Fill(mu_scaled,finalWeight);

    }
    
    h_rmass_0[iEta_bins]->Fill(wtm,finalWeight);
    h_pt_0   [iEta_bins]->Fill(pt,finalWeight);
    h_eta_0  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
    h_met_0  [iEta_bins]->Fill(met,finalWeight);
    h_iso30_0[iEta_bins]->Fill(iso30,finalWeight);
    h_d0sig_0[iEta_bins]->Fill(d0sig,finalWeight);
    h_ptDmt_0[iEta_bins]->Fill(pt/wtm,finalWeight);	   

    h_eta_d_0    [iEta_bins]->Fill(eta,finalWeight);
    h_phi_0      [iEta_bins]->Fill(phi,finalWeight);
    h_phi_met_0  [iEta_bins]->Fill(phimet,finalWeight);
    h_delta_phi_0[iEta_bins]->Fill(dPhi,finalWeight);
    h_mu_0       [iEta_bins]->Fill(mu_scaled,finalWeight);

    
    //if(iso30>0.1 && iso30<0.16){
    if(std::fabs(d0sig)>3.){
      h_rmass_1[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_1   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_1  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_met_1  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_1[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_1[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_1[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_1    [iEta_bins]->Fill(eta,finalWeight);
      h_phi_1      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_1  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_1[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_1       [iEta_bins]->Fill(mu_scaled,finalWeight);
    }
    
    if(std::fabs(d0sig)>4.){
      h_rmass_2[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_2   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_2  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_met_2  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_2[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_2[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_2[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_2    [iEta_bins]->Fill(eta,finalWeight);
      h_phi_2      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_2  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_2[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_2       [iEta_bins]->Fill(mu_scaled,finalWeight);

    }
    
    //if(iso30>0.22 && iso30<0.28){
    if(std::fabs(d0sig)>5.){
      h_rmass_3[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_3   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_3  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_met_3  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_3[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_3[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_3[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_3    [iEta_bins]->Fill(eta,finalWeight);
      h_phi_3      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_3  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_3[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_3       [iEta_bins]->Fill(mu_scaled,finalWeight);

    }
    
    //if(iso30>0.28 && iso30<0.34){
    if(std::fabs(d0sig)>5.){
      h_rmass_4[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_4   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_4  [iEta_bins]->Fill(eta,finalWeight);
      h_met_4  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_4[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_4[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_4[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_4    [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_phi_4      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_4  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_4[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_4       [iEta_bins]->Fill(mu_scaled,finalWeight);

    }
    
    //if(iso30>0.34 && iso30<0.4){
    if(std::fabs(d0sig)>6.){
      h_rmass_5[iEta_bins]->Fill(wtm,finalWeight);
      h_pt_5   [iEta_bins]->Fill(pt,finalWeight);
      h_eta_5  [iEta_bins]->Fill(std::fabs(eta),finalWeight);
      h_met_5  [iEta_bins]->Fill(met,finalWeight);
      h_iso30_5[iEta_bins]->Fill(iso30,finalWeight);
      h_d0sig_5[iEta_bins]->Fill(d0sig,finalWeight);
      h_ptDmt_5[iEta_bins]->Fill(pt/wtm,finalWeight);

      h_eta_d_5    [iEta_bins]->Fill(eta,finalWeight);
      h_phi_5      [iEta_bins]->Fill(phi,finalWeight);
      h_phi_met_5  [iEta_bins]->Fill(phimet,finalWeight);
      h_delta_phi_5[iEta_bins]->Fill(dPhi,finalWeight);
      h_mu_5       [iEta_bins]->Fill(mu_scaled,finalWeight);
    }
  
  }
      
  return;
}




void MyWZAnalysis::LumiNormalisation(Config config, int it_sf_eta, bool isData, double lumi_sample, double tt2p){

  
  double scale=lumi_sample/tt2p;

  //Luminosity normalisation
  if(!isData && !DoMultijetAnalysis){
    
    //Inclusive normalisation
    //double scale=(weightsSumWZ*lumi_sample)/tt2p;
    
    //double scaleRecoil=(weightsSumWZRecoil*lumi_sample)/tt2p;
    
    // if(BosonZSelection)  scaleRecoil=scale;
    // if(HasRecoilInfoAnalysis) m_recoilHistos->scalehists(scaleRecoil);
    
    if(h_pt[it_sf_eta]->Integral()!=0){
      
      // h_met_recoil[it_sf_eta]    ->Scale(1./(h_met_recoil->Integral()));    h_met_recoil->Scale(scaleRecoil);
      // h_rmass_recoil[it_sf_eta]  ->Scale(1./(h_rmass_recoil->Integral()));  h_rmass_recoil->Scale(scaleRecoil);
      
      // h_mu       ->Scale(1./(h_mu->Integral()));           
      // h_npvtx    ->Scale(1./(h_npvtx->Integral()));        
      // h_phi_met  ->Scale(1./(h_phi_met->Integral()));      	  
      // h_SumET    ->Scale(1./(h_SumET->Integral()));        
      // h_met      ->Scale(1./(h_met->Integral()));          
      // h_rmass      ->Scale(1./(h_rmass->Integral()));      
      // h_rmassLogCw ->Scale(1./(h_rmassLogCw->Integral())); 
      // PtRmass      ->Scale(1./(PtRmass->Integral()));
      
      h_mu[it_sf_eta]->Scale(scale);	    
      h_npvtx[it_sf_eta]->Scale(scale);	    
      h_phi_met[it_sf_eta]->Scale(scale);    
      h_SumET[it_sf_eta]->Scale(scale);	    
      h_met[it_sf_eta]->Scale(scale);   	    
      h_rmass[it_sf_eta]->Scale(scale);	    
      h_rmassLogCw[it_sf_eta]->Scale(scale); 
      PtRmass[it_sf_eta]->Scale(scale);      


      // h_pt        [it_sf_eta]->Scale(1./(h_pt->Integral()));               
      // h_eta       [it_sf_eta]->Scale(1./(h_eta->Integral()));              
      // h_eta_d     [it_sf_eta]->Scale(1./(h_eta_d->Integral()));
      // h_phi       [it_sf_eta]->Scale(1./(h_phi->Integral()));              
      // h_pt20Calo  [it_sf_eta]->Scale(1./(h_pt20Calo->Integral()));       	    
      // h_pt20Track [it_sf_eta]->Scale(1./(h_pt20Track->Integral()));      
      // h_delta_phi [it_sf_eta]->Scale(1./(h_delta_phi->Integral()));      
      // h_z0        [it_sf_eta]->Scale(1./(h_z0->Integral()));               
      // h_d0        [it_sf_eta]->Scale(1./(h_d0->Integral()));           
      
      h_pt[it_sf_eta]->Scale(scale);	  
      h_eta[it_sf_eta]->Scale(scale);	  
      h_eta_d[it_sf_eta]->Scale(scale);	  
      h_phi[it_sf_eta]->Scale(scale);	  
      h_pt20Calo[it_sf_eta]->Scale(scale); 
      h_pt20Track[it_sf_eta]->Scale(scale);
      h_delta_phi[it_sf_eta]->Scale(scale);
      h_z0[it_sf_eta]->Scale(scale);	  
      h_d0[it_sf_eta]->Scale(scale);       
      
      // h_SumET_eta1[it_sf_eta]->Scale(1./(h_SumET_eta1 ->Integral()));   
      // h_SumET_eta2[it_sf_eta]->Scale(1./(h_SumET_eta2 ->Integral()));   
      // h_SumET_eta3[it_sf_eta]->Scale(1./(h_SumET_eta3 ->Integral()));   
      // h_SumET_eta4[it_sf_eta]->Scale(1./(h_SumET_eta4 ->Integral()));   
      // h_SumET_eta5[it_sf_eta]->Scale(1./(h_SumET_eta5 ->Integral()));   

      h_SumET_eta1[it_sf_eta]->Scale(scale);
      h_SumET_eta2[it_sf_eta]->Scale(scale);
      h_SumET_eta3[it_sf_eta]->Scale(scale);
      h_SumET_eta4[it_sf_eta]->Scale(scale);
      h_SumET_eta5[it_sf_eta]->Scale(scale);
      
      //cw calculation
      // hMassCrec  [it_sf_eta]->Scale(1./hMassCrec ->Integral()); 
      // hMassCgen  [it_sf_eta]->Scale(1./hMassCgen ->Integral()); 
      // hMassCstay [it_sf_eta]->Scale(1./hMassCstay->Integral()); 
      // hMassCcome [it_sf_eta]->Scale(1./hMassCcome->Integral()); 
      
      hMassCrec  ->Scale(scale);
      hMassCgen  ->Scale(scale);
      hMassCstay ->Scale(scale);
      hMassCcome ->Scale(scale);
      
      // hEtaCrec  ->Scale(1./hEtaCrec ->Integral());  
      // hEtaCgen  ->Scale(1./hEtaCgen ->Integral());  
      // hEtaCstay ->Scale(1./hEtaCstay->Integral());  
      // hEtaCcome ->Scale(1./hEtaCcome->Integral());  
      
      hEtaCrec  ->Scale(scale);
      hEtaCgen  ->Scale(scale);
      hEtaCstay ->Scale(scale);
      hEtaCcome ->Scale(scale);
      
      // h_pt_boson_reco[it_sf_eta]->Scale(1./h_pt_boson_reco->Integral());     
      // h_pt_boson_reco_p[it_sf_eta]->Scale(1./h_pt_boson_reco_p->Integral()); 
      
      h_pt_boson_reco[it_sf_eta]->Scale(scale);  
      h_pt_boson_reco_p[it_sf_eta]->Scale(scale);
      
      // h_pt_boson   [it_sf_eta]->Scale(1./(h_pt_boson->Integral()));    	  
      // h_mass_boson [it_sf_eta]->Scale(1./(h_mass_boson->Integral()));   
      // h_eta_boson  [it_sf_eta]->Scale(1./(h_eta_boson->Integral()));      
      
      h_pt_boson[it_sf_eta]->Scale(scale);  
      h_mass_boson[it_sf_eta]->Scale(scale);
      h_eta_boson[it_sf_eta]->Scale(scale); 
      
    }
    //High mass normalisation
    //double scaleHigh=(weightsSumWZHigh*lumi_sample)/tt2p;
    double scaleHigh=lumi_sample/tt2p;
    // double scaleRecoilHigh=(weightsSumWZRecoilHigh*lumi_sample)/tt2p;
    
    // if(BosonZSelection)  scaleRecoilHigh=scale;
    // if(HasRecoilInfoAnalysis) m_recoilHistos->scalehists(scaleRecoilHigh);
    // if(h_pt_high->Integral()!=0){
      
    //   h_met_recoil_high   [it_sf_eta]->Scale(1./(h_met_recoil_high->Integral()));   h_met_recoil_high[it_sf_eta]->Scale(scaleRecoil);
    //   h_rmass_recoil_high [it_sf_eta]->Scale(1./(h_rmass_recoil_high->Integral()));       h_rmass_recoil_high[it_sf_eta]->Scale(scaleRecoil);
      
      // h_mu_high       [it_sf_eta]->Scale(1./(h_mu_high->Integral()));             
      // h_npvtx_high    [it_sf_eta]->Scale(1./(h_npvtx_high->Integral()));          
      // h_phi_met_high  [it_sf_eta]->Scale(1./(h_phi_met_high->Integral()));        
      // h_SumET_high    [it_sf_eta]->Scale(1./(h_SumET_high->Integral()));          
      // h_met_high      [it_sf_eta]->Scale(1./(h_met_high->Integral()));            
      // h_rmass_high      [it_sf_eta]->Scale(1./(h_rmass_high->Integral()));        
      // h_rmassLogCw_high [it_sf_eta]->Scale(1./(h_rmassLogCw_high->Integral()));   
      // PtRmass_high      [it_sf_eta]->Scale(1./(PtRmass_high->Integral()));        
      
      h_mu_high[it_sf_eta]->Scale(scaleHigh);		
      h_npvtx_high[it_sf_eta]->Scale(scaleHigh);	
      h_phi_met_high[it_sf_eta]->Scale(scaleHigh);	
      h_SumET_high[it_sf_eta]->Scale(scaleHigh);	
      h_met_high[it_sf_eta]->Scale(scaleHigh);   	
      h_rmass_high[it_sf_eta]->Scale(scaleHigh);	
      h_rmassLogCw_high[it_sf_eta]->Scale(scaleHigh);	
      PtRmass_high[it_sf_eta]->Scale(scaleHigh);       
      
      // h_pt_high        [it_sf_eta]->Scale(1./(h_pt_high->Integral()));        
      // h_eta_high       [it_sf_eta]->Scale(1./(h_eta_high->Integral()));       
      // h_eta_high_d     [it_sf_eta]->Scale(1./(h_eta_high_d->Integral()));
      // h_phi_high       [it_sf_eta]->Scale(1./(h_phi_high->Integral()));       
      // h_pt20Calo_high  [it_sf_eta]->Scale(1./(h_pt20Calo_high->Integral()));  
      // h_pt20Track_high [it_sf_eta]->Scale(1./(h_pt20Track_high->Integral())); 
      // h_delta_phi_high [it_sf_eta]->Scale(1./(h_delta_phi_high->Integral())); 
      // h_z0_high        [it_sf_eta]->Scale(1./(h_z0_high->Integral()));        
      // h_d0_high        [it_sf_eta]->Scale(1./(h_d0_high->Integral()));        
      
      h_pt_high[it_sf_eta]->Scale(scaleHigh);	   
      h_eta_high[it_sf_eta]->Scale(scaleHigh);	   
      h_eta_high_d[it_sf_eta]->Scale(scaleHigh);	   
      h_phi_high[it_sf_eta]->Scale(scaleHigh);	   
      h_pt20Calo_high[it_sf_eta]->Scale(scaleHigh);  
      h_pt20Track_high[it_sf_eta]->Scale(scaleHigh); 
      h_delta_phi_high[it_sf_eta]->Scale(scaleHigh); 
      h_z0_high[it_sf_eta]->Scale(scaleHigh);	   
      h_d0_high[it_sf_eta]->Scale(scaleHigh);        
      
      // h_pt_boson_high   [it_sf_eta]->Scale(1./(h_pt_boson_high->Integral()));   
      // h_mass_boson_high [it_sf_eta]->Scale(1./(h_mass_boson_high->Integral())); 
      // h_eta_boson_high  [it_sf_eta]->Scale(1./(h_eta_boson_high->Integral()));  
      
      h_pt_boson_high[it_sf_eta]->Scale(scaleHigh);  
      h_mass_boson_high[it_sf_eta]->Scale(scaleHigh);
      h_eta_boson_high[it_sf_eta]->Scale(scaleHigh);       
  }
    
  
  // weightsSumWZ=0.;      weightsSumWZRecoil=0.;
  // weightsSumWZHigh=0.;  weightsSumWZRecoilHigh=0.;
  
  //Luminosity normalisation
  if(!isData && DoMultijetAnalysis){
    
    //double scale=lumi_sample/tt2p;	
    
    h_rmass_i[it_sf_eta]->Scale(scale);
    h_rmass_0[it_sf_eta]->Scale(scale);
    h_rmass_1[it_sf_eta]->Scale(scale);
    h_rmass_2[it_sf_eta]->Scale(scale);
    h_rmass_3[it_sf_eta]->Scale(scale);
    h_rmass_4[it_sf_eta]->Scale(scale);
    h_rmass_5[it_sf_eta]->Scale(scale);
    
    h_pt_i[it_sf_eta]->Scale(scale);
    h_pt_0[it_sf_eta]->Scale(scale);
    h_pt_1[it_sf_eta]->Scale(scale);
    h_pt_2[it_sf_eta]->Scale(scale);
    h_pt_3[it_sf_eta]->Scale(scale);
    h_pt_4[it_sf_eta]->Scale(scale);
    h_pt_5[it_sf_eta]->Scale(scale);
    
    h_eta_i[it_sf_eta]->Scale(scale);	    
    h_eta_0[it_sf_eta]->Scale(scale);	    
    h_eta_1[it_sf_eta]->Scale(scale);	    
    h_eta_2[it_sf_eta]->Scale(scale);	    
    h_eta_3[it_sf_eta]->Scale(scale);	    
    h_eta_4[it_sf_eta]->Scale(scale);	    
    h_eta_5[it_sf_eta]->Scale(scale);	    
    
    h_met_i[it_sf_eta]->Scale(scale);
    h_met_0[it_sf_eta]->Scale(scale);
    h_met_1[it_sf_eta]->Scale(scale);
    h_met_2[it_sf_eta]->Scale(scale);
    h_met_3[it_sf_eta]->Scale(scale);
    h_met_4[it_sf_eta]->Scale(scale);
    h_met_5[it_sf_eta]->Scale(scale);
    
    h_iso30_i[it_sf_eta]->Scale(scale);
    h_iso30_0[it_sf_eta]->Scale(scale);
    h_iso30_1[it_sf_eta]->Scale(scale);
    h_iso30_2[it_sf_eta]->Scale(scale);
    h_iso30_3[it_sf_eta]->Scale(scale);
    h_iso30_4[it_sf_eta]->Scale(scale);
    h_iso30_5[it_sf_eta]->Scale(scale);
    
    h_d0sig_i[it_sf_eta]->Scale(scale);
    h_d0sig_0[it_sf_eta]->Scale(scale);
    h_d0sig_1[it_sf_eta]->Scale(scale);
    h_d0sig_2[it_sf_eta]->Scale(scale);
    h_d0sig_3[it_sf_eta]->Scale(scale);
    h_d0sig_4[it_sf_eta]->Scale(scale);
    h_d0sig_5[it_sf_eta]->Scale(scale);
    
    h_ptDmt_i[it_sf_eta]->Scale(scale);
    h_ptDmt_0[it_sf_eta]->Scale(scale);
    h_ptDmt_1[it_sf_eta]->Scale(scale);
    h_ptDmt_2[it_sf_eta]->Scale(scale);
    h_ptDmt_3[it_sf_eta]->Scale(scale);
    h_ptDmt_4[it_sf_eta]->Scale(scale);
    h_ptDmt_5[it_sf_eta]->Scale(scale);

    h_eta_d_i[it_sf_eta]->Scale(scale);	    
    h_eta_d_0[it_sf_eta]->Scale(scale);	    
    h_eta_d_1[it_sf_eta]->Scale(scale);	    
    h_eta_d_2[it_sf_eta]->Scale(scale);	    
    h_eta_d_3[it_sf_eta]->Scale(scale);	    
    h_eta_d_4[it_sf_eta]->Scale(scale);	    
    h_eta_d_5[it_sf_eta]->Scale(scale);	    

    h_phi_i[it_sf_eta]->Scale(scale);
    h_phi_0[it_sf_eta]->Scale(scale);
    h_phi_1[it_sf_eta]->Scale(scale);
    h_phi_2[it_sf_eta]->Scale(scale);
    h_phi_3[it_sf_eta]->Scale(scale);
    h_phi_4[it_sf_eta]->Scale(scale);
    h_phi_5[it_sf_eta]->Scale(scale);
    
    h_phi_met_i[it_sf_eta]->Scale(scale);
    h_phi_met_0[it_sf_eta]->Scale(scale);
    h_phi_met_1[it_sf_eta]->Scale(scale);
    h_phi_met_2[it_sf_eta]->Scale(scale);
    h_phi_met_3[it_sf_eta]->Scale(scale);
    h_phi_met_4[it_sf_eta]->Scale(scale);
    h_phi_met_5[it_sf_eta]->Scale(scale);
    
    h_delta_phi_i[it_sf_eta]->Scale(scale);
    h_delta_phi_0[it_sf_eta]->Scale(scale);
    h_delta_phi_1[it_sf_eta]->Scale(scale);
    h_delta_phi_2[it_sf_eta]->Scale(scale);
    h_delta_phi_3[it_sf_eta]->Scale(scale);
    h_delta_phi_4[it_sf_eta]->Scale(scale);
    h_delta_phi_5[it_sf_eta]->Scale(scale);
    
    h_mu_i[it_sf_eta]->Scale(scale);
    h_mu_0[it_sf_eta]->Scale(scale);
    h_mu_1[it_sf_eta]->Scale(scale);
    h_mu_2[it_sf_eta]->Scale(scale);
    h_mu_3[it_sf_eta]->Scale(scale);
    h_mu_4[it_sf_eta]->Scale(scale);
    h_mu_5[it_sf_eta]->Scale(scale);
    
  }//end norm MultiJetMaker
  
  return;
}

void MyWZAnalysis::WriteHist(int iSF_variations){

  newfile[iSF_variations]->cd();
  h_mu[iSF_variations]->Write("h_mu");  
  h_npvtx[iSF_variations]->Write("h_npvtx");	  
  h_phi_met[iSF_variations]->Write("h_phi_met");	  
  h_SumET[iSF_variations]->Write("h_SumET");	  
  h_met[iSF_variations]->Write("h_met");   	  
  h_rmass[iSF_variations]->Write("h_rmass");	  
  h_rmassLogCw[iSF_variations]->Write("h_rmassLogCw"); 
  PtRmass[iSF_variations]->Write("PtRmass");        
  
  h_met_recoil[iSF_variations]->Write("h_met_recoil");  	  
  h_rmass_recoil[iSF_variations]->Write("h_rmass_recoil");
  
  h_pt[iSF_variations]->Write("h_pt");	      
  h_eta[iSF_variations]->Write("h_eta");	      
  h_eta_d[iSF_variations]->Write("h_eta_d");	      
  h_phi[iSF_variations]->Write("h_phi");	      
  h_pt20Calo[iSF_variations]->Write("h_pt20Calo");  
  h_pt20Track[iSF_variations]->Write("h_pt20Track"); 
  h_delta_phi[iSF_variations]->Write("h_delta_phi"); 
  h_z0[iSF_variations]->Write("h_z0");	      
  h_d0[iSF_variations]->Write("h_d0");
  h_pt_boson_reco[iSF_variations]->Write("h_pt_boson_reco");
  h_pt_boson_reco_p[iSF_variations]->Write("h_pt_boson_reco_p");
  
  h_pt_boson[iSF_variations]->Write("h_pt_boson_truth");
  h_mass_boson[iSF_variations]->Write("h_mass_boson_truth");
  h_eta_boson[iSF_variations]->Write("h_eta_boson_truth"); 
  
  h_mu_high[iSF_variations]->Write("h_mu_high");	  
  h_npvtx_high[iSF_variations]->Write("h_npvtx_high");	  
  h_phi_met_high[iSF_variations]->Write("h_phi_met_high");	  
  h_SumET_high[iSF_variations]->Write("h_SumET_high");	  
  h_met_high[iSF_variations]->Write("h_met_high");   	  
  h_rmass_high[iSF_variations]->Write("h_rmass_high");	  
  h_rmassLogCw_high[iSF_variations]->Write("h_rmassLogCw_high"); 
  PtRmass_high[iSF_variations]->Write("PtRmass_high");        
  
  h_met_recoil_high[iSF_variations]->Write("h_met_recoil_high");  	  
  h_rmass_recoil_high[iSF_variations]->Write("h_rmass_recoil_high");
  
  h_pt_high[iSF_variations]->Write("h_pt_high");	      
  h_eta_high[iSF_variations]->Write("h_eta_high");	      
  h_eta_high_d[iSF_variations]->Write("h_eta_d_high");	      
  h_phi_high[iSF_variations]->Write("h_phi_high");	      
  h_pt20Calo_high[iSF_variations]->Write("h_pt20Calo_high");  
  h_pt20Track_high[iSF_variations]->Write("h_pt20Track_high"); 
  h_delta_phi_high[iSF_variations]->Write("h_delta_phi_high"); 
  h_z0_high[iSF_variations]->Write("h_z0_high");	      
  h_d0_high[iSF_variations]->Write("h_d0_high"); 
  
  h_pt_boson_high[iSF_variations]->Write("h_pt_boson_truth_high");	
  h_mass_boson_high[iSF_variations]->Write("h_mass_boson_truth_high");
  h_eta_boson_high[iSF_variations]->Write("h_eta_boson_truth_high"); 
    
  h_SumET_eta1[iSF_variations]->Write("h_SumET_eta1");
  h_SumET_eta2[iSF_variations]->Write("h_SumET_eta2");
  h_SumET_eta3[iSF_variations]->Write("h_SumET_eta3");
  h_SumET_eta4[iSF_variations]->Write("h_SumET_eta4");
  h_SumET_eta5[iSF_variations]->Write("h_SumET_eta5");
  
  hPileUp ->Write("hPileUp");
  h_truth_mass_nocuts->Write("h_truth_mass_nocuts");
  h_nm->Write("h_nm");

  hMassCrec  ->Write("hMassCrec");
  hMassCgen  ->Write("hMassCgen");
  hMassCstay ->Write("hMassCsta");
  hMassCcome ->Write("hMassCcome");
  
  hEtaCrec  ->Write("hEtaCrec");
  hEtaCgen  ->Write("hEtaCgen");
  hEtaCstay ->Write("hEtaCsta");
  hEtaCcome ->Write("hEtaCcome");
  

  //TProfiles and cutflow
  if(iSF_variations==0){
    h_cutflow->Write("h_cutflow");
    h_cutflow_high->Write("h_cutflow_high");

    p_mwt_Tot ->Write("p_mwt_Tot");
    p_mwt_Rec ->Write("p_mwt_Rec");
    p_mwt_KF  ->Write("p_mwt_KF" );
    p_mwt_MC  ->Write("p_mwt_MC" );
    
    p_mwt_Iso           ->Write("p_mwt_Iso"          );
    p_mwt_Iso_sys_UP    ->Write("p_mwt_Iso_sys_UP"   );
    p_mwt_Iso_sys_DOWN  ->Write("p_mwt_Iso_sys_DOWN" );
    p_mwt_Iso_stat_UP   ->Write("p_mwt_Iso_stat_UP"  );
    p_mwt_Iso_stat_DOWN ->Write("p_mwt_Iso_stat_DOWN");
    
    p_mwt_ID            ->Write("p_mwt_ID"          );
    p_mwt_ID_sys_UP     ->Write("p_mwt_ID_sys_UP"   );
    p_mwt_ID_sys_DOWN   ->Write("p_mwt_ID_sys_DOWN" );
    p_mwt_ID_stat_UP    ->Write("p_mwt_ID_stat_UP"  );
    p_mwt_ID_stat_DOWN  ->Write("p_mwt_ID_stat_DOWN");
    
    p_mwt_Tri           ->Write("p_mwt_Tri"          );
    p_mwt_Tri_sys_UP    ->Write("p_mwt_Tri_sys_UP"   );
    p_mwt_Tri_sys_DOWN  ->Write("p_mwt_Tri_sys_DOWN" );
    p_mwt_Tri_stat_UP   ->Write("p_mwt_Tri_stat_UP"  );
    p_mwt_Tri_stat_DOWN ->Write("p_mwt_Tri_stat_DOWN");
    
    p_mwt_ttv           ->Write("p_mwt_ttva"          );
    p_mwt_ttv_sys_UP    ->Write("p_mwt_ttva_sys_UP"   );
    p_mwt_ttv_sys_DOWN  ->Write("p_mwt_ttva_sys_DOWN" );
    p_mwt_ttv_stat_UP   ->Write("p_mwt_ttva_stat_UP"  );
    p_mwt_ttv_stat_DOWN ->Write("p_mwt_ttva_stat_DOWN");
    
    p_mwt_PU           ->Write("p_mwt_PU"     );
    p_mwt_PU_UP        ->Write("p_mwt_PU_UP"  );
    p_mwt_PU_DOWN      ->Write("p_mwt_PU_DOWN");
    
    p_mwt_jvt          ->Write("p_mwt_jvt_PU"     );
    p_mwt_jvt_UP       ->Write("p_mwt_jvt_PU_UP"  ); 
    p_mwt_jvt_DOWN     ->Write("p_mwt_jvt_PU_DOWN");   
    
    p_m_Tot ->Write("p_m_Tot");
    p_m_Rec ->Write("p_m_Rec");
    p_m_KF  ->Write("p_m_KF" );
    p_m_MC  ->Write("p_m_MC" );
    
    p_m_Iso           ->Write("p_m_Iso"          );
    p_m_Iso_sys_UP    ->Write("p_m_Iso_sys_UP"   );
    p_m_Iso_sys_DOWN  ->Write("p_m_Iso_sys_DOWN" );
    p_m_Iso_stat_UP   ->Write("p_m_Iso_stat_UP"  );
    p_m_Iso_stat_DOWN ->Write("p_m_Iso_stat_DOWN");
    
    p_m_ID            ->Write("p_m_ID"          );
    p_m_ID_sys_UP     ->Write("p_m_ID_sys_UP"   );
    p_m_ID_sys_DOWN   ->Write("p_m_ID_sys_DOWN" );
    p_m_ID_stat_UP    ->Write("p_m_ID_stat_UP"  );
    p_m_ID_stat_DOWN  ->Write("p_m_ID_stat_DOWN");
    
    p_m_Tri           ->Write("p_m_Tri"          );
    p_m_Tri_sys_UP    ->Write("p_m_Tri_sys_UP"   );
    p_m_Tri_sys_DOWN  ->Write("p_m_Tri_sys_DOWN" );
    p_m_Tri_stat_UP   ->Write("p_m_Tri_stat_UP"  );
    p_m_Tri_stat_DOWN ->Write("p_m_Tri_stat_DOWN");
    
    p_m_ttv           ->Write("p_m_ttva"          );
    p_m_ttv_sys_UP    ->Write("p_m_ttva_sys_UP"   );
    p_m_ttv_sys_DOWN  ->Write("p_m_ttva_sys_DOWN" );
    p_m_ttv_stat_UP   ->Write("p_m_ttva_stat_UP"  );
    p_m_ttv_stat_DOWN ->Write("p_m_ttva_stat_DOWN");
    
    p_m_PU           ->Write("p_m_PU"     );
    p_m_PU_UP        ->Write("p_m_PU_UP"  );
    p_m_PU_DOWN      ->Write("p_m_PU_DOWN");
    
    p_m_jvt          ->Write("p_m_jvt_PU"     );
    p_m_jvt_UP       ->Write("p_m_jvt_PU_UP"  ); 
    p_m_jvt_DOWN     ->Write("p_m_jvt_PU_DOWN");   
    
    
    p_pt_Tot ->Write("p_pt_Tot");
    p_pt_Rec ->Write("p_pt_Rec");
    p_pt_KF  ->Write("p_pt_KF" );
    p_pt_MC  ->Write("p_pt_MC" );
    
    p_pt_Iso           ->Write("p_pt_Iso"          );
    p_pt_Iso_sys_UP    ->Write("p_pt_Iso_sys_UP"   );
    p_pt_Iso_sys_DOWN  ->Write("p_pt_Iso_sys_DOWN" );
    p_pt_Iso_stat_UP   ->Write("p_pt_Iso_stat_UP"  );
    p_pt_Iso_stat_DOWN ->Write("p_pt_Iso_stat_DOWN");
    
    p_pt_ID            ->Write("p_pt_ID"          );
    p_pt_ID_sys_UP     ->Write("p_pt_ID_sys_UP"   );
    p_pt_ID_sys_DOWN   ->Write("p_pt_ID_sys_DOWN" );
    p_pt_ID_stat_UP    ->Write("p_pt_ID_stat_UP"  );
    p_pt_ID_stat_DOWN  ->Write("p_pt_ID_stat_DOWN");
    
    p_pt_Tri           ->Write("p_pt_Tri"          );
    p_pt_Tri_sys_UP    ->Write("p_pt_Tri_sys_UP"   );
    p_pt_Tri_sys_DOWN  ->Write("p_pt_Tri_sys_DOWN" );
    p_pt_Tri_stat_UP   ->Write("p_pt_Tri_stat_UP"  );
    p_pt_Tri_stat_DOWN ->Write("p_pt_Tri_stat_DOWN");
    
    p_pt_ttv           ->Write("p_pt_ttva"          );
    p_pt_ttv_sys_UP    ->Write("p_pt_ttva_sys_UP"   );
    p_pt_ttv_sys_DOWN  ->Write("p_pt_ttva_sys_DOWN" );
    p_pt_ttv_stat_UP   ->Write("p_pt_ttva_stat_UP"  );
    p_pt_ttv_stat_DOWN ->Write("p_pt_ttva_stat_DOWN");
    
    p_pt_PU           ->Write("p_pt_PU"     );
    p_pt_PU_UP        ->Write("p_pt_PU_UP"  );
    p_pt_PU_DOWN      ->Write("p_pt_PU_DOWN");
    
    p_pt_jvt          ->Write("p_pt_jvt_PU"     );
    p_pt_jvt_UP       ->Write("p_pt_jvt_PU_UP"  ); 
    p_pt_jvt_DOWN     ->Write("p_pt_jvt_PU_DOWN");   
    
    p_eta_Tot ->Write("p_eta_Tot");
    p_eta_Rec ->Write("p_eta_Rec");
    p_eta_KF  ->Write("p_eta_KF" );
    p_eta_MC  ->Write("p_eta_MC" );
    
    p_eta_Iso           ->Write("p_eta_Iso"          );
    p_eta_Iso_sys_UP    ->Write("p_eta_Iso_sys_UP"   );
    p_eta_Iso_sys_DOWN  ->Write("p_eta_Iso_sys_DOWN" );
    p_eta_Iso_stat_UP   ->Write("p_eta_Iso_stat_UP"  );
    p_eta_Iso_stat_DOWN ->Write("p_eta_Iso_stat_DOWN");
    
    p_eta_ID            ->Write("p_eta_ID"          );
    p_eta_ID_sys_UP     ->Write("p_eta_ID_sys_UP"   );
    p_eta_ID_sys_DOWN   ->Write("p_eta_ID_sys_DOWN" );
    p_eta_ID_stat_UP    ->Write("p_eta_ID_stat_UP"  );
    p_eta_ID_stat_DOWN  ->Write("p_eta_ID_stat_DOWN");
    
    p_eta_Tri           ->Write("p_eta_Tri"          );
    p_eta_Tri_sys_UP    ->Write("p_eta_Tri_sys_UP"   );
    p_eta_Tri_sys_DOWN  ->Write("p_eta_Tri_sys_DOWN" );
    p_eta_Tri_stat_UP   ->Write("p_eta_Tri_stat_UP"  );
    p_eta_Tri_stat_DOWN ->Write("p_eta_Tri_stat_DOWN");
    
    p_eta_ttv           ->Write("p_eta_ttva"          );
    p_eta_ttv_sys_UP    ->Write("p_eta_ttva_sys_UP"   );
    p_eta_ttv_sys_DOWN  ->Write("p_eta_ttva_sys_DOWN" );
    p_eta_ttv_stat_UP   ->Write("p_eta_ttva_stat_UP"  );
    p_eta_ttv_stat_DOWN ->Write("p_eta_ttva_stat_DOWN");
    
    p_eta_PU           ->Write("p_eta_PU"     );
    p_eta_PU_UP        ->Write("p_eta_PU_UP"  );
    p_eta_PU_DOWN      ->Write("p_eta_PU_DOWN");
    
    p_eta_jvt          ->Write("p_eta_jvt_PU"     );
    p_eta_jvt_UP       ->Write("p_eta_jvt_PU_UP"  ); 
    p_eta_jvt_DOWN     ->Write("p_eta_jvt_PU_DOWN");
  }

     
  delete h_mu[iSF_variations];	  
  delete h_npvtx[iSF_variations];	  
  delete h_phi_met[iSF_variations];	  
  delete h_SumET[iSF_variations];	  
  delete h_met[iSF_variations];   	  
  delete h_rmass[iSF_variations];	  
  delete h_rmassLogCw[iSF_variations]; 
  delete PtRmass[iSF_variations];
  
  delete h_met_recoil[iSF_variations];
  delete h_rmass_recoil[iSF_variations];
  
  delete h_pt[iSF_variations];	  
  delete h_eta[iSF_variations];	  
  delete h_eta_d[iSF_variations];	  
  delete h_phi[iSF_variations];	  
  delete h_pt20Calo[iSF_variations];   
  delete h_pt20Track[iSF_variations];  
  delete h_delta_phi[iSF_variations];  
  delete h_z0[iSF_variations];	  
  delete h_d0[iSF_variations];
  delete h_pt_boson_reco[iSF_variations];      
  delete h_pt_boson_reco_p[iSF_variations];
  
  delete h_pt_boson[iSF_variations];	  
  delete h_mass_boson[iSF_variations]; 
  delete h_eta_boson[iSF_variations];

  //high mass
  delete h_mu_high[iSF_variations];	  
  delete h_npvtx_high[iSF_variations];	  
  delete h_phi_met_high[iSF_variations];	  
  delete h_SumET_high[iSF_variations];	  
  delete h_met_high[iSF_variations];   	  
  delete h_rmass_high[iSF_variations];	  
  delete h_rmassLogCw_high[iSF_variations]; 
  delete PtRmass_high[iSF_variations];  
  
  delete h_met_recoil_high[iSF_variations];
  delete h_rmass_recoil_high[iSF_variations];      
      
  delete h_pt_high[iSF_variations];	  
  delete h_eta_high[iSF_variations];	  
  delete h_eta_high_d[iSF_variations];	  
  delete h_phi_high[iSF_variations];	  
  delete h_pt20Calo_high[iSF_variations];   
  delete h_pt20Track_high[iSF_variations];  
  delete h_delta_phi_high[iSF_variations];  
  delete h_z0_high[iSF_variations];	  
  delete h_d0_high[iSF_variations];         
  
  delete h_pt_boson_high[iSF_variations];	  
  delete h_mass_boson_high[iSF_variations]; 
  delete h_eta_boson_high[iSF_variations]; 
  
  delete h_SumET_eta1[iSF_variations];
  delete h_SumET_eta2[iSF_variations];
  delete h_SumET_eta3[iSF_variations];
  delete h_SumET_eta4[iSF_variations];
  delete h_SumET_eta5[iSF_variations];


  //no cuts
  // delete hPileUp;
  // delete h_truth_mass_nocuts;
  // delete h_nm;

  //Cutflow
  // delete h_cutflow;
  // delete h_cutflow_high;

  //cw factor
  // delete hMassCrec;
  // delete hMassCgen;
  // delete hMassCstay;
  // delete hMassCcome;
  
  // delete hEtaCrec;
  // delete hEtaCgen;
  // delete hEtaCstay;
  // delete hEtaCcome;

  // //profiles 
  // delete p_mwt_Tot;
  // delete p_mwt_Iso;
  // delete p_mwt_ID ;
  // delete p_mwt_Tri;
  // delete p_mwt_Rec;
  // delete p_mwt_ttv;
  // delete p_mwt_KF;
  // delete p_mwt_PU;
  // delete p_mwt_MC;
  
  // delete p_pt_Tot;
  // delete p_pt_Iso;
  // delete p_pt_ID;
  // delete p_pt_Tri;
  // delete p_pt_Rec;
  // delete p_pt_ttv;
  // delete p_pt_KF;
  // delete p_pt_PU;
  // delete p_pt_MC;
  
  // delete p_eta_Tot;
  // delete p_eta_Iso;
  // delete p_eta_ID;
  // delete p_eta_Tri;
  // delete p_eta_Rec;
  // delete p_eta_ttv;
  // delete p_eta_KF;
  // delete p_eta_PU;
  // delete p_eta_MC;
  
  newfile[iSF_variations]->Close();
  
  return;
}

void MyWZAnalysis::WriteHistMultijet(int iEta_bin){

  newfile[0]->cd();

  h_rmass_i[iEta_bin]->Write();
  h_rmass_0[iEta_bin]->Write();
  h_rmass_1[iEta_bin]->Write();
  h_rmass_2[iEta_bin]->Write();
  h_rmass_3[iEta_bin]->Write();
  h_rmass_4[iEta_bin]->Write();
  h_rmass_5[iEta_bin]->Write();
                          
  h_pt_i[iEta_bin]->Write();  
  h_pt_0[iEta_bin]->Write();
  h_pt_1[iEta_bin]->Write();
  h_pt_2[iEta_bin]->Write();
  h_pt_3[iEta_bin]->Write();
  h_pt_4[iEta_bin]->Write();
  h_pt_5[iEta_bin]->Write();
  
  h_eta_i[iEta_bin]->Write();	    
  h_eta_0[iEta_bin]->Write();	    
  h_eta_1[iEta_bin]->Write();	    
  h_eta_2[iEta_bin]->Write();	    
  h_eta_3[iEta_bin]->Write();	    
  h_eta_4[iEta_bin]->Write();	    
  h_eta_5[iEta_bin]->Write();	    
  
  h_met_i[iEta_bin]->Write();
  h_met_0[iEta_bin]->Write();
  h_met_1[iEta_bin]->Write();
  h_met_2[iEta_bin]->Write();
  h_met_3[iEta_bin]->Write();
  h_met_4[iEta_bin]->Write();
  h_met_5[iEta_bin]->Write();
  
  h_iso30_i[iEta_bin]->Write();
  h_iso30_0[iEta_bin]->Write();
  h_iso30_1[iEta_bin]->Write();
  h_iso30_2[iEta_bin]->Write();
  h_iso30_3[iEta_bin]->Write();
  h_iso30_4[iEta_bin]->Write();
  h_iso30_5[iEta_bin]->Write();

  h_d0sig_i[iEta_bin]->Write();
  h_d0sig_0[iEta_bin]->Write();
  h_d0sig_1[iEta_bin]->Write();
  h_d0sig_2[iEta_bin]->Write();
  h_d0sig_3[iEta_bin]->Write();
  h_d0sig_4[iEta_bin]->Write();
  h_d0sig_5[iEta_bin]->Write();

  h_ptDmt_i[iEta_bin]->Write();
  h_ptDmt_0[iEta_bin]->Write();
  h_ptDmt_1[iEta_bin]->Write();
  h_ptDmt_2[iEta_bin]->Write();
  h_ptDmt_3[iEta_bin]->Write();
  h_ptDmt_4[iEta_bin]->Write();
  h_ptDmt_5[iEta_bin]->Write();

  h_eta_d_i[iEta_bin]->Write();	    
  h_eta_d_0[iEta_bin]->Write();	    
  h_eta_d_1[iEta_bin]->Write();	    
  h_eta_d_2[iEta_bin]->Write();	    
  h_eta_d_3[iEta_bin]->Write();	    
  h_eta_d_4[iEta_bin]->Write();	    
  h_eta_d_5[iEta_bin]->Write();	    
  
  h_phi_i[iEta_bin]->Write();
  h_phi_0[iEta_bin]->Write();
  h_phi_1[iEta_bin]->Write();
  h_phi_2[iEta_bin]->Write();
  h_phi_3[iEta_bin]->Write();
  h_phi_4[iEta_bin]->Write();
  h_phi_5[iEta_bin]->Write();
  
  h_phi_met_i[iEta_bin]->Write();
  h_phi_met_0[iEta_bin]->Write();
  h_phi_met_1[iEta_bin]->Write();
  h_phi_met_2[iEta_bin]->Write();
  h_phi_met_3[iEta_bin]->Write();
  h_phi_met_4[iEta_bin]->Write();
  h_phi_met_5[iEta_bin]->Write();
  
  h_delta_phi_i[iEta_bin]->Write();
  h_delta_phi_0[iEta_bin]->Write();
  h_delta_phi_1[iEta_bin]->Write();
  h_delta_phi_2[iEta_bin]->Write();
  h_delta_phi_3[iEta_bin]->Write();
  h_delta_phi_4[iEta_bin]->Write();
  h_delta_phi_5[iEta_bin]->Write();
  
  h_mu_i[iEta_bin]->Write();
  h_mu_0[iEta_bin]->Write();
  h_mu_1[iEta_bin]->Write();
  h_mu_2[iEta_bin]->Write();
  h_mu_3[iEta_bin]->Write();
  h_mu_4[iEta_bin]->Write();
  h_mu_5[iEta_bin]->Write();
  
  delete h_rmass_i[iEta_bin];
  delete h_rmass_0[iEta_bin];
  delete h_rmass_1[iEta_bin];
  delete h_rmass_2[iEta_bin];
  delete h_rmass_3[iEta_bin];
  delete h_rmass_4[iEta_bin];
  delete h_rmass_5[iEta_bin];

  delete h_pt_i[iEta_bin];
  delete h_pt_0[iEta_bin];
  delete h_pt_1[iEta_bin];
  delete h_pt_2[iEta_bin];
  delete h_pt_3[iEta_bin];
  delete h_pt_4[iEta_bin];
  delete h_pt_5[iEta_bin];

  delete h_eta_i[iEta_bin];
  delete h_eta_0[iEta_bin];	    
  delete h_eta_1[iEta_bin];	    
  delete h_eta_2[iEta_bin];	    
  delete h_eta_3[iEta_bin];	    
  delete h_eta_4[iEta_bin];	    
  delete h_eta_5[iEta_bin];	    

  delete h_met_i[iEta_bin];
  delete h_met_0[iEta_bin];
  delete h_met_1[iEta_bin];
  delete h_met_2[iEta_bin];
  delete h_met_3[iEta_bin];
  delete h_met_4[iEta_bin];
  delete h_met_5[iEta_bin];

  delete h_iso30_i[iEta_bin];
  delete h_iso30_0[iEta_bin];
  delete h_iso30_1[iEta_bin];
  delete h_iso30_2[iEta_bin];
  delete h_iso30_3[iEta_bin];
  delete h_iso30_4[iEta_bin];
  delete h_iso30_5[iEta_bin];

  delete h_d0sig_i[iEta_bin];
  delete h_d0sig_0[iEta_bin];
  delete h_d0sig_1[iEta_bin];
  delete h_d0sig_2[iEta_bin];
  delete h_d0sig_3[iEta_bin];
  delete h_d0sig_4[iEta_bin];
  delete h_d0sig_5[iEta_bin];

  delete h_ptDmt_i[iEta_bin];
  delete h_ptDmt_0[iEta_bin];
  delete h_ptDmt_1[iEta_bin];
  delete h_ptDmt_2[iEta_bin];
  delete h_ptDmt_3[iEta_bin];
  delete h_ptDmt_4[iEta_bin];
  delete h_ptDmt_5[iEta_bin];

  delete h_eta_d_i[iEta_bin];	    
  delete h_eta_d_0[iEta_bin];	    
  delete h_eta_d_1[iEta_bin];	    
  delete h_eta_d_2[iEta_bin];	    
  delete h_eta_d_3[iEta_bin];	    
  delete h_eta_d_4[iEta_bin];	    
  delete h_eta_d_5[iEta_bin];	    
  
  delete h_phi_i[iEta_bin];
  delete h_phi_0[iEta_bin];
  delete h_phi_1[iEta_bin];
  delete h_phi_2[iEta_bin];
  delete h_phi_3[iEta_bin];
  delete h_phi_4[iEta_bin];
  delete h_phi_5[iEta_bin];
  
  delete h_phi_met_i[iEta_bin];
  delete h_phi_met_0[iEta_bin];
  delete h_phi_met_1[iEta_bin];
  delete h_phi_met_2[iEta_bin];
  delete h_phi_met_3[iEta_bin];
  delete h_phi_met_4[iEta_bin];
  delete h_phi_met_5[iEta_bin];
  
  delete h_delta_phi_i[iEta_bin];
  delete h_delta_phi_0[iEta_bin];
  delete h_delta_phi_1[iEta_bin];
  delete h_delta_phi_2[iEta_bin];
  delete h_delta_phi_3[iEta_bin];
  delete h_delta_phi_4[iEta_bin];
  delete h_delta_phi_5[iEta_bin];
    
  delete h_mu_i[iEta_bin];
  delete h_mu_0[iEta_bin];
  delete h_mu_1[iEta_bin];
  delete h_mu_2[iEta_bin];
  delete h_mu_3[iEta_bin];
  delete h_mu_4[iEta_bin];
  delete h_mu_5[iEta_bin];
  
  if(iEta_bin==11) newfile[0]->Close();
  
  return;
}


void MyWZAnalysis::SetConfiguration(Config config){

  if(config.DataYears=="2015+2016"){
    DataYears2015p2016 = true;
    DataYears2017 = false;
    DataYears2018 = false;
  }else if(config.DataYears=="2017"){
    DataYears2015p2016 = false;
    DataYears2017 = true;
    DataYears2018 = false;
  }else if(config.DataYears=="2018"){
    DataYears2015p2016 = false;
    DataYears2017 = false;
    DataYears2018 = true;
  }

  if(config.OnlyInclusive=="True")      OnlyInclusiveAnalysis = true;
  else if(config.OnlyInclusive!="True") OnlyInclusiveAnalysis = false;
   
  if(config.DoMultijet=="True")       DoMultijetAnalysis = true;
  else if(config.DoMultijet!="True")  DoMultijetAnalysis = false;
       
  if(config.MultiFitRegion=="Signal"){
    MJFitRegionSignal = true;
    MJFitRegionFitR1  = false;
    MJFitRegionFitR2  = false;
  }else if(config.MultiFitRegion=="FitR1"){
    MJFitRegionSignal = false;
    MJFitRegionFitR1  = true;
    MJFitRegionFitR2  = false;    
  }else if(config.MultiFitRegion=="FitR2"){
    MJFitRegionSignal = false;
    MJFitRegionFitR1  = false;
    MJFitRegionFitR2  = true;
  }
  
  if(config.Dod0Cut=="True")        Dod0CutAnalysis = true;
  else if(config.Dod0Cut!="True")   Dod0CutAnalysis = false;

  if(config.Applyd0Weight=="True")  Applyd0WeightAnalysis = true;
  else if(config.Dod0Cut!="True")   Applyd0WeightAnalysis = false;

  if(config.WZSelection=="zmumu"){     
    BosonZSelection = true;	  
    BosonWplusSelection = false; 
    BosonWminusSelection = false;
  }else if(config.WZSelection=="wplus"){
    BosonZSelection = false;	  
    BosonWplusSelection = true; 
    BosonWminusSelection = false;
  }else if(config.WZSelection=="wminus"){
    BosonZSelection = false;	  
    BosonWplusSelection = false; 
    BosonWminusSelection = true;
  }

  if(config.HasRecoilInfo=="True")       HasRecoilInfoAnalysis = true;
  else if(config.HasRecoilInfo!="True")  HasRecoilInfoAnalysis = false;

  if(config.OnTheFlyPileUp=="True")      OnTheFlyPileUpAnalysis = true;
  else if(config.OnTheFlyPileUp!="True") OnTheFlyPileUpAnalysis = false;

  if(config.SETCalibration=="True")       SETCalibrationAnalysis = true;
  else if(config.SETCalibration!="True")  SETCalibrationAnalysis = false;

  if(config.InsituCorrection=="True")      InsituCorrectionAnalysis = true;
  else if(config.InsituCorrection!="True") InsituCorrectionAnalysis = false;

  if(config.ResolResponse=="True")      ResolResponseAnalysis = true;
  else if(config.ResolResponse!="True") ResolResponseAnalysis = false;

  if(config.TruthMatching=="True")      TruthMatchingAnalysis = true;
  else if(config.TruthMatching!="True") TruthMatchingAnalysis = false;

  if(config.SFVariations=="True")       SFVariationsAnalysis = true;
  else if(config.SFVariations!="True")  SFVariationsAnalysis = false;       

  return;
}
