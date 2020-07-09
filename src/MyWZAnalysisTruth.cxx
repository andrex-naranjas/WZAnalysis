#define MyWZAnalysisTruth_cxx
#include "MyWZAnalysisTruth.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include "TProfile.h"
#include "TProfile.h"
#include "TProfile.h"
#include "TMath.h"
#include <iostream>

#include <vector>
#include <vector>
#include <vector>

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

#include "TLorentzVector.h"
#include <vector>

#include "BosonReco.h"

void MyWZAnalysisTruth::Loop()
{
   if (fChain == 0) return;
   
   std::string wzchannel;
   if(config.WZSelection=="zmumu"){
     wzchannel="z";
   }else if(config.WZSelection=="wplus"){
     wzchannel="wplus";
   }else if(config.WZSelection=="wminus"){
     wzchannel="wminus";}

   //ouput directory according year
   std::string dirYear="";
   if(config.DataYears=="2015+2016") dirYear="2015p2016/";
   if(config.DataYears=="2017") dirYear="2017/";

   //Cutflow
   TString chanflowdir, cutflowdir, cutflowdir2;
   if(config.WZSelection=="zmumu") chanflowdir="mumu";
   if(config.WZSelection=="wplus" || config.WZSelection=="wminus") chanflowdir="munu";
   if((nameOfSample=="wminmunu_15"|| nameOfSample=="zmumu_15")){
     cutflowdir=chanflowdir+"/cutflow";
   }else if(config.DataYears=="2015+2016"){
     cutflowdir=chanflowdir+"_2016/cutflow";
     cutflowdir2=chanflowdir+"_2015/cutflow";
   }else if(config.DataYears=="2017"){
     cutflowdir=chanflowdir+"_2017/cutflow";
     cutflowdir2=chanflowdir+"_2017/cutflow";
   }

   //skim normalisation
   float tt=0.; double tt2=0.;   
   fChain_sumWeights->SetBranchAddress("totalEventsWeighted",&tt);
   int ntree = fChain_sumWeights->GetEntries();
   for(Int_t i=0;i<ntree;i++){
     fChain_sumWeights->GetEntry(i);
     tt2=tt2+tt;
   }
   
   float tt2p=tt2/1.;

   //define lepton charge
   float chargeSelection = 0.;
   if(config.WZSelection == "wplus"){
     chargeSelection = 1.;
   }else if(config.WZSelection == "wminus"){
     chargeSelection = -1.;
   }
   
   //reco matching
   std::unordered_map<int,std::shared_ptr<BosonReco>> m_myRecoMap;
   typedef std::unordered_map<int,std::shared_ptr<BosonReco>>::const_iterator Itr;

   Long64_t nReco=0, nTruth = 0, mapCounter=0;

   if(config.RecoMatching=="True"){
     
     ULong64_t R_eventNumber;
     UInt_t    R_runNumber, R_mcChannelNumber;
     
     Double_t  R_KFactor_weight;
     Float_t   R_weight_mc, R_weight_pileup;
     
     //references to pointers, thanks a lot c++ :dd
     std::vector<float>  ref_R_mu_pt;
     std::vector<float>  *R_mu_pt = &ref_R_mu_pt;

     std::vector<float>  ref_R_mu_eta;
     std::vector<float>  *R_mu_eta = &ref_R_mu_eta;

     std::vector<float>  ref_R_mu_phi;
     std::vector<float>  *R_mu_phi = &ref_R_mu_phi;

     std::vector<float>  ref_R_mu_charge;
     std::vector<float>  *R_mu_charge = &ref_R_mu_charge;
     
     Float_t R_met_met;
     Float_t R_met_phi;
     
     nReco  = fChain_reco->GetEntries();
     
     fChain_reco->SetBranchAddress("eventNumber", &R_eventNumber);
     fChain_reco->SetBranchAddress("runNumber",   &R_runNumber);
     fChain_reco->SetBranchAddress("mcChannelNumber", &R_mcChannelNumber);
     
     fChain_reco->SetBranchAddress("weight_mc", &R_weight_mc);
     fChain_reco->SetBranchAddress("weight_pileup", &R_weight_pileup);
     fChain_reco->SetBranchAddress("weight_KFactor", &R_KFactor_weight);

     fChain_reco->SetBranchAddress("mu_pt",  &R_mu_pt);
     fChain_reco->SetBranchAddress("mu_eta", &R_mu_eta);
     fChain_reco->SetBranchAddress("mu_phi", &R_mu_phi);
     fChain_reco->SetBranchAddress("mu_charge",&R_mu_charge);

     fChain_reco->SetBranchAddress("met_met", &R_met_met);
     fChain_reco->SetBranchAddress("met_phi", &R_met_phi);
     
     if(config.NumberOfEvents!=-1) nReco=config.NumberOfEvents;//custom # events

     for(Long64_t reco_entry=0; reco_entry <nReco; reco_entry++){//loop over truth tree
       
       fChain_reco->GetEntry(reco_entry);
       if(reco_entry % 100000 == 0) std::cout<<"Processing reco event  "<<reco_entry<<"/"<<nReco<<std::endl;
       
       std::shared_ptr<BosonReco> recoEvent = std::shared_ptr<BosonReco>(new BosonReco());
       
       recoEvent->eventNumberTruth = R_eventNumber;
       recoEvent->runNumberTruth = R_runNumber;
       recoEvent->mcChannelNumber = R_mcChannelNumber;
       
       recoEvent->weight_mc = R_weight_mc;
       recoEvent->weight_pileup = R_weight_pileup;	    
       recoEvent->KFactor_weight = R_KFactor_weight;

       recoEvent->met_met = R_met_met;
       recoEvent->met_phi = R_met_phi;
       
       if(R_mu_pt->size()==1){
	 recoEvent-> mu_pt  = R_mu_pt->at(0);
	 recoEvent-> mu_eta = R_mu_eta->at(0);
	 recoEvent-> mu_phi = R_mu_phi->at(0);
	 recoEvent-> mu_charge = R_mu_charge->at(0);
	 m_myRecoMap.insert(std::make_pair(R_eventNumber,recoEvent));//fill the map
       }
       
     }//truth loop
   }//if not Data and sampleFlag


   TFile* newfile;
   newfile = new TFile((config.OutputFileDir + "Files/" + dirYear +  "Truth/" + nameOfSample + "_"+ wzchannel + ".root").c_str(),"RECREATE");
   
   TH1::SetDefaultSumw2(kTRUE);
   Int_t nbins=7;
   Double_t x0 = TMath::Log10(50);
   Double_t x1 = TMath::Log10(100);
   Double_t x2 = 0.349485002;
   
   Double_t bins[]={x0,x1,x1+1*(0.25),x1+2*(0.25),x1+3*(0.25),3.,3.+x2,3.+2*x2};
   //Double_t bins[]={x0,2.5,2.5+x2,2.5+2*x2,3.};
      
   // Log binning in lepton pT
   printf("  calculating logarithmic pT binning...\n");
   float Ptmin = 25000.; float Ptmax = 3000000.;
   Ptmin = log10(Ptmin); Ptmax = log10(Ptmax);
   const int nPtEdge1 = 101;
   float PtBins1[nPtEdge1];
   float ptWidth = (Ptmax - Ptmin)/(nPtEdge1-1);
   for (int i = 0; i < nPtEdge1; i++){
     PtBins1[i  ] = Ptmin+(i  )*ptWidth;
     PtBins1[i  ] = pow(10,PtBins1[i  ]);
     //printf("bin  i=%3i  edge= %9.4f   width=%9.2f\n",i,PtBins1[i],ptWidth);
   }
   printf("  calculating logarithmic M binning...\n");
   float Massmin = 50000.; float Massmax = 7000000.;
   Massmin = log10(Massmin); Massmax = log10(Massmax);
   const int nMassEdge = 101;
   float MassBins[nMassEdge];
   float binWidth = (Massmax - Massmin)/(nMassEdge-1);
   for (int i = 0; i < nMassEdge; i++){
     MassBins[i  ] = Massmin+(i  )*binWidth;
     MassBins[i  ] = pow(10,MassBins[i  ]);
     //printf("bin  i=%3i  lo=%9.4f  hi=%9.4f \n",i,MassBins[i-1],MassBins[i]);
   }
   
  //W variables
  TH1D* h_WZ_m    = new  TH1D("h_WZ_m"  , "mass; m^{truth} [GeV]; Entries ",nMassEdge-1,MassBins );
  //TH1D* h_WZ_m    = new  TH1D("h_WZ_m"  , "mass; m^{truth} [GeV]; Entries ",200,50000,7000000);
  TH1D* h_WZ_mt   = new  TH1D("h_WZ_mt"  , "mass; m^{truth} [GeV]; Entries ",nMassEdge-1,MassBins );
  TH1D* h_WZ_pt   = new  TH1D("h_WZ_pt" , "pt; p_{T} [GeV]; Entries", 50, 0.,150);	   
  TH1D* h_WZ_eta  = new  TH1D("h_WZ_eta", "#eta; Entries", 50,-6.,6.);			   
  TH1D* h_WZ_phi  = new  TH1D("h_WZ_phi", "#phi [rad]; Entries", 50, -3.141592, 3.141592);   
	                                     
  //BORN LEVEL			     
  TH1D* h_WZ_dilep_m_born    = new  TH1D("h_WZ_dilep_m_born"  , "mass; m^{truth} [GeV]; Entries ",nMassEdge-1,MassBins );
  TH1D* h_WZ_dilep_pt_born   = new  TH1D("h_WZ_dilep_pt_born" , "pt; p_{T} [GeV]; Entries",nPtEdge1-1, PtBins1);
  TH1D* h_WZ_dilep_eta_born  = new  TH1D("h_WZ_dilep_eta_born", "#eta; Entries", 50,-6.,6.);			     
  TH1D* h_WZ_dilep_phi_born  = new  TH1D("h_WZ_dilep_phi_born", "#phi [rad]; Entries", 50, -3.141592, 3.141592);  
	                                     
  //BARE LEVEL			     
  TH1D* h_WZ_dilep_m_bare    = new  TH1D("h_WZ_dilep_m_bare"  , "mass; m^{truth} [GeV]; Entries ",nMassEdge-1,MassBins );
  TH1D* h_WZ_dilep_pt_bare   = new  TH1D("h_WZ_dilep_pt_bare" , "pt; p_{T} [GeV]; Entries",nPtEdge1-1, PtBins1);
  TH1D* h_WZ_dilep_eta_bare  = new  TH1D("h_WZ_dilep_eta_bare", "#eta; Entries", 50,-6.,6.);			  
  TH1D* h_WZ_dilep_phi_bare  = new  TH1D("h_WZ_dilep_phi_bare", "#phi [rad]; Entries", 50, -3.141592, 3.141592); 
  
  //DRESSED LEVEL			     
  TH1D* h_WZ_dilep_m_dres    = new  TH1D("h_WZ_dilep_m_dres"  , "mass; m^{truth} [GeV]; Entries",nMassEdge-1,MassBins );
  TH1D* h_WZ_dilep_pt_dres   = new  TH1D("h_WZ_dilep_pt_dres" , "pt; p_{T} [GeV]; Entries",nPtEdge1-1, PtBins1);
  TH1D* h_WZ_dilep_eta_dres  = new  TH1D("h_WZ_dilep_eta_dres", "#eta; Entries", 50,-6.,6.);			 
  TH1D* h_WZ_dilep_phi_dres  = new  TH1D("h_WZ_dilep_phi_dres", "#phi [rad]; Entries", 50, -3.141592, 3.141592); 

  TH1F* hTruMassW = new TH1F("hTruMassW"   , "Tru: all M    ; m^{tru} [GeV]; Entries"   ,nMassEdge-1,MassBins );

  // TH2D *PtmWT= new TH2D("PtmWT","",nbins,bins,nbins,bins);/*10,40000.,120000.);*/
  // TH2D *PtmW= new TH2D("PtmW","Pt Vs Invariant Mass",nbins,bins,nbins,bins);/*10,40000.,120000.);*/

  TH1D *h_mw  = new TH1D("h_mw","W Boson Truth InvariantMass ;M_{W}^{Tru}  [GeV]    ;Events ",nMassEdge-1,MassBins);//150000.
  TH1D *h_mwt = new TH1D("h_mwt","W Boson Truth transverse mass Mass; M_{W,T}  [GeV]; Events",nMassEdge-1,MassBins);
  TH1D *h_pwt = new TH1D("h_pwt","W Boson Truth transverse momentum; M_{W,T}  [GeV]; Events",nMassEdge-1,MassBins);
  
  // double eta_bins[]={0.0,0.21,0.42,0.64,0.86,1.08,1.30,1.52,1.74,1.96,2.18,2.5}; //andres
  double eta_bins[]={0.0,0.21,0.42,0.63,0.84,1.05,1.37,1.52,1.74,1.95,2.18,2.4};//8,7 TeV binning

  TH1D *h_eta = new TH1D("h_eta","Muon Pseudorapidity    ;#eta^{#mu}; Events",11,eta_bins);//22,-2.6,2.6);(// 13/52 bins default
  //TH1D *h_eta = new TH1D("h_eta","Muon Pseudorapidity    ;#eta^{#mu}; Events",22,-5.,5.);// 22,-2.6,2.6) 13/52 bins default
  
  TH2D *h_mw_pwt  = new TH2D("h_mw_pwt","Invariant Mass vs Ptw",nMassEdge-1,MassBins , nPtEdge1-1, PtBins1);
  TH2D *h_mwt_pwt = new TH2D("h_mwt_pwt","Transverse mass vs Ptw",nMassEdge-1,MassBins,nPtEdge1-1, PtBins1);
  TH2D *h_mw_pt   = new TH2D("h_mw_pt","Transverse mass vs Pt muon",nMassEdge-1,MassBins , nPtEdge1-1, PtBins1);
  TH2D *h_pwt_pt  = new TH2D("h_pwt_pt","Ptw Vs Pt#mu",nPtEdge1-1, PtBins1, nPtEdge1-1, PtBins1);

  TH2D *h_mw_mwt= new TH2D("h_mw_mwt","Invariant Mass vs transverse mass",nMassEdge-1,MassBins,nMassEdge-1,MassBins);
  TH2D *h_mw_deltaPhi= new TH2D("h_mw_deltaPhi","Invariant Mass vs #Delta#phi",nMassEdge-1,MassBins,nMassEdge-1,MassBins);/*10,40000.,120000.);*/
  TH2D *h_mwt_deltaPhi= new TH2D("h_mwt_deltaPhi","Delta Phi Vs Invariant Mass",nMassEdge-1,MassBins,nMassEdge-1,MassBins);/*10,40000.,120000.);*/

  TH1D *h_mLog=new TH1D("h_mLog","W Boson Truth Invariant Mass; M_{W}  [GeV]; Events",nbins,bins);
  TH1D *h_mLog_high = (TH1D*)h_mLog->Clone("h_mLog_high");

  TH1D *h_mwtLogCw=new TH1D("h_mwtLogCw","W Boson Truth Invariant Mass; M_{W,T}  [GeV]; Events",nbins,bins);
  TH1D *h_mwtLogCw_high = (TH1D*)h_mwtLogCw->Clone("h_mwtLogCw_high");

  TProfile *profmax_mwt = new TProfile("profmax_mwt","",nbins,bins,-6.6,-1.2);
  TProfile *profmin_mwt = new TProfile("profmin_mwt","",nbins,bins,-6.6,-1.2);
  TProfile *profmax_m = new TProfile("profmax_m","",nbins,bins,-6.6,-1.2);
  TProfile *profmin_m = new TProfile("profmin_m","",nbins,bins,-6.6,-1.2);
  
  profmax_mwt->SetErrorOption("s");
  profmin_mwt->SetErrorOption("s");
  profmax_m->SetErrorOption("s");
  profmin_m->SetErrorOption("s");

  //Cw factor calculation
  unsigned int nxbins = config.xBinsCw.size()-1;
  double xbins[nxbins+1];
  for(int i = 0; i<nxbins+1;i++) xbins[i]=config.xBinsCw.at(i);
  
  TH1D* hMassCrec    = new  TH1D("hMassCrec"   , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  TH1D* hMassCgen    = new  TH1D("hMassCgen"   , "mass; m_{T}^{reco} [GeV]; Entries ", nxbins,xbins);
  TH1D* hMassCstayG  = new  TH1D("hMassCstayG" , "mass; m_{T}^{reco} [GeV]; Entries ",nxbins,xbins);
  TH1D* hMassCleave  = new  TH1D("hMassCleave" , "mass; m_{T}^{reco} [GeV]; Entries ",nxbins,xbins);
  
  TH1D* hEtaCrec    = new  TH1D("hEtaCrec"   , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  TH1D* hEtaCgen    = new  TH1D("hEtaCgen"   , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  TH1D* hEtaCstayG  = new  TH1D("hEtaCstayG" , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  TH1D* hEtaCleave  = new  TH1D("hEtaCleave" , "eta; #eta_{#mu}^{reco}; Entries ", 11,eta_bins);
  
  Float_t W_m=0.0, lW_m=0.0, W_pt=0.0, mu_pt=0.0, lmu_pt = 0., nu_pt=0.0, mu_eta=0.0, lW_pt = 0.;
  Float_t nu_eta=0.0, mu_phi=0.0, nu_phi=0.0, mwt=0.0, lmwt=0.0, deltaphi = 0.;
  Float_t x1pdf=0., x2pdf=0., xmin=0., xmax=0.;
  
  //normalization
  double weightsSumWZ=0., weightsSumWZ_nocuts=0.;
  double weightsSum_5=0., weightsSum_6=0.;
  double weightsSum_7=0., weightsSum_8=0.;
  
  //Cw factor error calculation
  int iRecBin = 0, iGenBin=0;
  int iRecBinEta = 0, iGenBinEta=0;
  
  int counT=0;
  Long64_t nentries = fChain->GetEntries();
  if(config.NumberOfEvents!=-1) nentries=config.NumberOfEvents;//custom # events
  std::cout<<nentries<<"    "<<cut_pt1<<"   "<<cut_eta<<"    "<<cut_down<<"   "<<cut_up<<std::endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    if(jentry % 100000 == 0) std::cout<<"Processing event  "<<jentry<<"/"<<nentries<<std::endl;
    
    //truth-reco event matching
    bool recoSelection = false;
    if(/*sampleFlag &&*/ config.RecoMatching=="True"){
      Itr i = m_myRecoMap.find(eventNumber);
      if (i != m_myRecoMap.end()){ 
	std::shared_ptr<BosonReco> myRecoEvent = (*i).second;
	mapCounter=mapCounter+1;	    
	
	//Cw factor calculation
	double metReco        = myRecoEvent->met_met;//neutrino
	double metphiReco     = myRecoEvent->met_phi;	 
	double pt_mu_reco     = myRecoEvent->mu_pt;
	double eta_mu_reco    = myRecoEvent->mu_eta;
	double phi_mu_reco    = myRecoEvent->mu_phi;
	double charge_mu_reco = myRecoEvent->mu_charge;
	
	float dPhiReco=std::fabs(metphiReco - phi_mu_reco);
	if (dPhiReco>3.14159265) dPhiReco=2.*3.14159265-dPhiReco;
	double wtmReco = std::sqrt( 2.*( (pt_mu_reco) * metReco )*(1.-std::cos(dPhiReco)));
	
	iRecBin = hMassCgen->FindBin(wtmReco);
	iRecBinEta = hEtaCgen->FindBin(std::fabs(eta_mu_reco));
	
	if(pt_mu_reco > cut_pt1 && metReco > cut_met1){
	  if(std::fabs(eta_mu_reco) < cut_eta ){
	    if((charge_mu_reco) == chargeSelection){
	      if(wtmReco > cut_mwt1){
		recoSelection = true;
	      }
	    }
	  }
	  
	}
	
      } //check reco event has a truth event associated with it 
    }
    
    
    //avoid double counting for the inclusive sample
    bool massMatch = true;
    if(config.OnlyInclusive!="True" && (nameOfSample=="wplusmunu" || nameOfSample=="wminmunu") && MC_WZ_m>=120000.)
      massMatch=false;
    
    Int_t nx1 = PDFinfo_X1->size();
    Int_t nx2 = PDFinfo_X2->size();
    
    if(nx1 == 1 && nx2 == 1){	
      x1pdf=PDFinfo_X1->at(0);
      x2pdf=PDFinfo_X2->at(0);	
    }else if (nx1 != 1 || nx2 != 1){
      
      std::cout<<nx1<<"   "<<nx2<<"    :Warning"<<std::endl;
    }
    
    if(x1pdf >= x2pdf){
      xmax=x1pdf;
      xmin=x2pdf;
    }else if(x1pdf < x2pdf){
      xmax=x2pdf;
      xmin=x1pdf;
    }
    
    
    xmax=std::log10(xmax);
    xmin=std::log10(xmin);
    
    if (xmax < -7.0 || xmin < -7.0) std::cout<<xmax<<"   "<<xmin<<std::endl;
    
    // if ( weight_mc > 0) weight_mc =  1.0;
    // if ( weight_mc < 0) weight_mc = -1.0;
    if(massMatch){
      float ScaleGen = weight_pileup * KFactor_weight_truth * weight_mc;
      weightsSumWZ_nocuts+=weight_pileup*weight_mc;    

      //pileup weight normalisation
      weightsSum_5+= KFactor_weight_truth * weight_mc;
      weightsSum_6+= weight_pileup * KFactor_weight_truth * weight_mc;	      

      h_WZ_m    -> Fill(MC_WZ_m/1.,ScaleGen);
      h_WZ_pt   -> Fill(MC_WZ_pt/1.,ScaleGen);
      h_WZ_eta  -> Fill(MC_WZ_eta/1.,ScaleGen);
      h_WZ_phi  -> Fill(MC_WZ_phi/1.,ScaleGen);
      
      h_WZ_dilep_m_born    -> Fill(MC_WZ_dilep_m_born/1.,ScaleGen);
      h_WZ_dilep_pt_born   -> Fill(MC_WZ_dilep_pt_born/1.,ScaleGen);
      h_WZ_dilep_eta_born  -> Fill(MC_WZ_dilep_eta_born/1.,ScaleGen);
      h_WZ_dilep_phi_born  -> Fill(MC_WZ_dilep_phi_born/1.,ScaleGen);
      
      h_WZ_dilep_m_bare    -> Fill(MC_WZ_dilep_m_bare/1.,ScaleGen);
      h_WZ_dilep_pt_bare   -> Fill(MC_WZ_dilep_pt_bare/1.,ScaleGen);
      h_WZ_dilep_eta_bare  -> Fill(MC_WZ_dilep_eta_bare/1.,ScaleGen);
      h_WZ_dilep_phi_bare  -> Fill(MC_WZ_dilep_phi_bare/1.,ScaleGen);
      
      h_WZ_dilep_m_dres    -> Fill(MC_WZ_dilep_m_dres/1.,ScaleGen);
      h_WZ_dilep_pt_dres   -> Fill(MC_WZ_dilep_pt_dres/1.,ScaleGen);
      h_WZ_dilep_eta_dres  -> Fill(MC_WZ_dilep_eta_dres/1.,ScaleGen);
      h_WZ_dilep_phi_dres  -> Fill(MC_WZ_dilep_phi_dres/1.,ScaleGen);
	
      W_m=MC_WZ_m;
      lW_m=std::log10(W_m*0.001);
      W_pt=MC_WZ_pt;
      mu_pt=MC_WZmu_el_pt_born;
      nu_pt=MC_WZneutrino_pt_born;
      mu_eta=MC_WZmu_el_eta_born;
      nu_eta=MC_WZneutrino_eta_born;
      mu_phi=MC_WZmu_el_phi_born;
      nu_phi=MC_WZneutrino_phi_born;
      
      lW_pt=std::log10(W_pt*0.001);
      
      float dPhi=std::fabs(mu_phi-nu_phi);
      if (dPhi>3.141592) dPhi=2.*(3.141592)-dPhi;
      
      mwt = std::sqrt(2*(mu_pt*nu_pt)*(1-std::cos(dPhi)));
      lmwt = std::log10(mwt*0.001);
      lmu_pt=std::log10(mu_pt*0.001);
      
      h_WZ_mt -> Fill(mwt/1.,ScaleGen);

      //cw calculation
      iGenBin = hMassCgen->FindBin(mwt);
      iGenBinEta = hEtaCgen->FindBin(std::fabs(mu_eta));	
      
      if( W_m >= cut_down && W_m <= cut_up || true){
	hTruMassW->Fill(W_m/1.,ScaleGen);
	if(mu_pt > cut_pt1 && nu_pt > cut_met1){//kine cuts 1
	  if(fabs(mu_eta) < cut_eta){//kine cuts 2
	    if(mwt > cut_mwt1){

	      //pileup weight normalisation
	      weightsSum_7+= KFactor_weight_truth * weight_mc;
	      weightsSum_8+= weight_pileup * KFactor_weight_truth * weight_mc;
	      
	      h_mw->Fill(W_m/1.,ScaleGen);
	      h_mwt->Fill(mwt/1.,ScaleGen);	
	      h_eta->Fill(std::fabs(mu_eta),ScaleGen);
	      h_pwt->Fill(W_pt/1.,ScaleGen);
	      h_mw_pwt->Fill(W_m/1.,W_pt/1.,ScaleGen);
	      h_mwt_pwt->Fill(mwt/1.,W_pt/1.,ScaleGen);
	      h_mw_pt->Fill(W_m/1.,mu_pt/1.,ScaleGen);
	      h_pwt_pt->Fill(W_pt/1.,mu_pt/1.,ScaleGen);
	      h_mw_mwt->Fill(W_m/1.,mwt/1.,ScaleGen);
	      h_mw_deltaPhi->Fill(W_m/1.,dPhi,ScaleGen);
	      h_mwt_deltaPhi->Fill(mwt/1.,dPhi,ScaleGen);
	      
	      h_mLog->Fill(lW_m,ScaleGen);
	      h_mwtLogCw->Fill(lmwt,ScaleGen);
	      profmax_mwt->Fill(lmwt,xmax,ScaleGen);
	      profmin_mwt->Fill(lmwt,xmin,ScaleGen);
	      profmax_m->Fill(lW_m,xmax,ScaleGen);
	      profmin_m->Fill(lW_m,xmin,ScaleGen);
	      
	      //Cw factor calculation
	      hMassCgen->Fill(mwt,ScaleGen);
	      if( iGenBin == iRecBin && recoSelection  ) hMassCstayG->Fill(mwt,ScaleGen);
	      if( iGenBin == iRecBin && !recoSelection ) hMassCleave->Fill(mwt,ScaleGen);
	      if( iGenBin != iRecBin                   ) hMassCleave->Fill(mwt,ScaleGen);
	      hEtaCgen->Fill(std::fabs(mu_eta),ScaleGen);
	      if( iGenBinEta == iRecBinEta && recoSelection  ) hEtaCstayG->Fill(std::fabs(mu_eta),ScaleGen);
	      if( iGenBinEta == iRecBinEta && !recoSelection ) hEtaCleave->Fill(std::fabs(mu_eta),ScaleGen);
	      if( iGenBinEta != iRecBinEta                   ) hEtaCleave->Fill(std::fabs(mu_eta),ScaleGen);
	      
	      weightsSumWZ+=weight_pileup*weight_mc;
	      
	    }//transverse mass 
	  }//kinecuts1
	}//kinecuts2	
      }//mass range cut_pt1      
      
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
  }//event loop


  double pileup_norm_nocuts = weightsSum_5/weightsSum_6;
  double pileup_norm = weightsSum_7/weightsSum_8;
  
  double scale_nocuts=pileup_norm_nocuts*lumi_sample / tt2p;
  double scale=pileup_norm*lumi_sample / tt2p;
  
  std::cout<<scale<<"    "<<weightsSumWZ<<"    "<<weightsSumWZ_nocuts<<"    "<<lumi_sample<<"    "<<tt2p<<"        pileup ->:"<< pileup_norm  <<"       Normalisation"<<std::endl;
  

   h_WZ_m    -> Scale(scale_nocuts);
   h_WZ_mt   -> Scale(scale_nocuts);
   h_WZ_pt   -> Scale(scale_nocuts);
   h_WZ_eta  -> Scale(scale_nocuts);
   h_WZ_phi  -> Scale(scale_nocuts); 
   
   h_WZ_dilep_m_born    -> Scale(scale_nocuts);
   h_WZ_dilep_pt_born   -> Scale(scale_nocuts);
   h_WZ_dilep_eta_born  -> Scale(scale_nocuts);
   h_WZ_dilep_phi_born  -> Scale(scale_nocuts);
                                    					                                   
   h_WZ_dilep_m_bare    -> Scale(scale_nocuts);
   h_WZ_dilep_pt_bare   -> Scale(scale_nocuts);
   h_WZ_dilep_eta_bare  -> Scale(scale_nocuts);
   h_WZ_dilep_phi_bare  -> Scale(scale_nocuts);	  
   
   h_WZ_dilep_m_dres    -> Scale(scale_nocuts);
   h_WZ_dilep_pt_dres   -> Scale(scale_nocuts);
   h_WZ_dilep_eta_dres  -> Scale(scale_nocuts);
   h_WZ_dilep_phi_dres  -> Scale(scale_nocuts);

   hMassCgen  ->Scale(scale);
   hMassCstayG->Scale(scale);
   hMassCleave->Scale(scale);
   
   hEtaCgen  ->Scale(scale);
   hEtaCstayG->Scale(scale);
   hEtaCleave->Scale(scale);

   h_mw ->Scale(scale);
   h_mwt->Scale(scale);
   h_eta->Scale(scale);
   h_pwt->Scale(scale);
   
   newfile->cd();

   hTruMassW->Write("hTruMassW");
   h_mw ->Write("h_mw");
   h_mwt->Write("h_mwtw");
   h_eta->Write("h_etaw");
   h_pwt->Write("h_pwtw");
   h_mw_pwt->Write("h_mw_pwt");
   h_mwt_pwt->Write("h_mwt_pwt");
   h_mw_pt->Write("h_mw_pt");
   h_pwt_pt->Write("h_pwt_pt");
   h_mw_mwt->Write("h_mw_mwt");
   h_mw_deltaPhi->Write("h_mw_deltaPhi");
   h_mwt_deltaPhi->Write("h_mwt_deltaPhi");
   h_mLog->Write("h_mLog");
   h_mLog_high->Write("h_mLog_high");
   h_mwtLogCw->Write("h_mwtLogCw");
   h_mwtLogCw_high->Write("h_mwtLogCw_high");

   h_WZ_m    -> Write("h_WZ_m"  );
   h_WZ_mt   -> Write("h_WZ_mt" );
   h_WZ_pt   -> Write("h_WZ_pt" ); 
   h_WZ_eta  -> Write("h_WZ_eta"); 
   h_WZ_phi  -> Write("h_WZ_phi");   
   
   h_WZ_dilep_m_born    -> Write("h_WZ_dilep_m_born"  );
   h_WZ_dilep_pt_born   -> Write("h_WZ_dilep_pt_born" );   
   h_WZ_dilep_eta_born  -> Write("h_WZ_dilep_eta_born");   
   h_WZ_dilep_phi_born  -> Write("h_WZ_dilep_phi_born");

   h_WZ_dilep_m_bare    -> Write("h_WZ_dilep_m_bare"  );
   h_WZ_dilep_pt_bare   -> Write("h_WZ_dilep_pt_bare" ); 
   h_WZ_dilep_eta_bare  -> Write("h_WZ_dilep_eta_bare"); 
   h_WZ_dilep_phi_bare  -> Write("h_WZ_dilep_phi_bare");
 
   h_WZ_dilep_m_dres    -> Write("h_WZ_dilep_m_dres"  ); 
   h_WZ_dilep_pt_dres   -> Write("h_WZ_dilep_pt_dres" ); 
   h_WZ_dilep_eta_dres  -> Write("h_WZ_dilep_eta_dres"); 
   h_WZ_dilep_phi_dres  -> Write("h_WZ_dilep_phi_dres"); 

   hMassCrec    ->Write("hMassCrec");
   hMassCgen    ->Write("hMassCgen");
   hMassCstayG  ->Write("hMassCstayG");
   hMassCleave  ->Write("hMassCleave");

   hEtaCrec    ->Write("hEtaCrec");
   hEtaCgen    ->Write("hEtaCgen");
   hEtaCstayG  ->Write("hEtaCstayG");
   hEtaCleave  ->Write("hEtaCleave");

   profmax_mwt->Write("profmax_mwt");
   profmin_mwt->Write("profmin_mwt");
   profmax_m->Write("profmax_m");
   profmin_m->Write("profmin_m");

   delete hTruMassW;
   delete h_mw;
   delete h_mwt;
   delete h_eta;
   delete h_pwt;
   delete h_mw_pwt;
   delete h_mwt_pwt;
   delete h_mw_pt;
   delete h_pwt_pt;
   delete h_mw_mwt;
   delete h_mw_deltaPhi;
   delete h_mwt_deltaPhi;
   delete h_mLog;
   delete h_mLog_high;
   delete h_mwtLogCw;
   delete h_mwtLogCw_high;

   delete h_WZ_m;
   delete h_WZ_mt;
   delete h_WZ_pt;
   delete h_WZ_eta;
   delete h_WZ_phi;
   
   delete h_WZ_dilep_m_born;
   delete h_WZ_dilep_pt_born;
   delete h_WZ_dilep_eta_born;
   delete h_WZ_dilep_phi_born;

   delete h_WZ_dilep_m_bare;
   delete h_WZ_dilep_pt_bare;
   delete h_WZ_dilep_eta_bare;
   delete h_WZ_dilep_phi_bare;
 
   delete h_WZ_dilep_m_dres;
   delete h_WZ_dilep_pt_dres;
   delete h_WZ_dilep_eta_dres;
   delete h_WZ_dilep_phi_dres;

   delete hMassCrec;
   delete hMassCgen;
   delete hMassCstayG;
   delete hMassCleave;

   delete hEtaCrec;
   delete hEtaCgen;
   delete hEtaCstayG;
   delete hEtaCleave;

   delete profmax_mwt;
   delete profmin_mwt;
   delete profmax_m;
   delete profmin_m;

   newfile->Close("");
   delete newfile;
}
