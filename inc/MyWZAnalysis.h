#ifndef MyWZAnalysis_h
#define MyWZAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH2.h>
#include <TH1.h>
#include <TF1.h>

#include "vector"
#include "vector"
#include "vector"
#include <vector>

#include "ConfigSettings.h"
#include "TLorentzVector.h"

class MyWZAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TTree          *fChain_sumWeights;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain_truth;   //!pointer to the analyzed TTree or TChain

   std::string    nameOfSample;
   std::string    nameOfSystematic;
   // std::string    mass_range;
   Double_t       lumi_sample;
   Config config;

   Double_t       cut_eta;
   Double_t       cut_pt1;
   Double_t       cut_pt2;
   Double_t       cut_met1;
   Double_t       cut_met2;
   Double_t       cut_mwt1;
   Double_t       cut_mwt2;

   // Declaration of leaf types
   Float_t         weight_mc;
   Float_t         weight_pileup;
   Float_t         weight_leptonSF;
   Float_t         weight_globalLeptonTriggerSF;
   Float_t         weight_oldTriggerSF;
   //new variables
   Double_t        weight_KFactor;
   Float_t         bornMass_KFactor;
   Float_t         weight_jvt;
   //Float_t         weight_bTagSF_77;
   Float_t         weight_pileup_UP;
   Float_t         weight_pileup_DOWN;
   Float_t         weight_leptonSF_EL_SF_Trigger_UP;
   Float_t         weight_leptonSF_EL_SF_Trigger_DOWN;
   Float_t         weight_leptonSF_EL_SF_Reco_UP;
   Float_t         weight_leptonSF_EL_SF_Reco_DOWN;
   Float_t         weight_leptonSF_EL_SF_ID_UP;
   Float_t         weight_leptonSF_EL_SF_ID_DOWN;
   Float_t         weight_leptonSF_EL_SF_Isol_UP;
   Float_t         weight_leptonSF_EL_SF_Isol_DOWN;
   Float_t         weight_leptonSF_MU_SF_Trigger_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_Trigger_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Trigger_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_Trigger_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_DOWN;
   Float_t         weight_indiv_SF_EL_Trigger;
   Float_t         weight_indiv_SF_EL_Trigger_UP;
   Float_t         weight_indiv_SF_EL_Trigger_DOWN;
   Float_t         weight_indiv_SF_EL_Reco;
   Float_t         weight_indiv_SF_EL_Reco_UP;
   Float_t         weight_indiv_SF_EL_Reco_DOWN;
   Float_t         weight_indiv_SF_EL_ID;
   Float_t         weight_indiv_SF_EL_ID_UP;
   Float_t         weight_indiv_SF_EL_ID_DOWN;
   Float_t         weight_indiv_SF_EL_Isol;
   Float_t         weight_indiv_SF_EL_Isol_UP;
   Float_t         weight_indiv_SF_EL_Isol_DOWN;
   Float_t         weight_globalLeptonTriggerSF_EL_Trigger_UP;
   Float_t         weight_globalLeptonTriggerSF_EL_Trigger_DOWN;
   Float_t         weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP;
   Float_t         weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN;
   Float_t         weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP;
   Float_t         weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN;
   Float_t         weight_oldTriggerSF_EL_Trigger_UP;
   Float_t         weight_oldTriggerSF_EL_Trigger_DOWN;
   Float_t         weight_oldTriggerSF_MU_Trigger_STAT_UP;
   Float_t         weight_oldTriggerSF_MU_Trigger_STAT_DOWN;
   Float_t         weight_oldTriggerSF_MU_Trigger_SYST_UP;
   Float_t         weight_oldTriggerSF_MU_Trigger_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_ID;
   Float_t         weight_indiv_SF_MU_ID_STAT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_Isol;
   Float_t         weight_indiv_SF_MU_Isol_STAT_UP;
   Float_t         weight_indiv_SF_MU_Isol_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_Isol_SYST_UP;
   Float_t         weight_indiv_SF_MU_Isol_SYST_DOWN;
   //new variables
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_TTVA;
   Float_t         weight_indiv_SF_MU_TTVA_STAT_UP;
   Float_t         weight_indiv_SF_MU_TTVA_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_TTVA_SYST_UP;
   Float_t         weight_indiv_SF_MU_TTVA_SYST_DOWN;
   Float_t         weight_jvt_UP;
   Float_t         weight_jvt_DOWN;
   ULong64_t       eventNumber;
   UInt_t          runNumber;
   UInt_t          randomRunNumber;
   UInt_t          mcChannelNumber;
   Float_t         mu;
   UInt_t          backgroundFlags;
   UInt_t          hasBadMuon;
   unsigned long long weight_pileup_hash;

   Double_t        weight_kfactor_sys_ALPHAS__1down;
   Double_t        weight_kfactor_sys_ALPHAS__1up;
   Double_t        weight_kfactor_sys_BEAM_ENERGY__1down;
   Double_t        weight_kfactor_sys_BEAM_ENERGY__1up;
   Double_t        weight_kfactor_sys_CHOICE_HERAPDF20;
   Double_t        weight_kfactor_sys_CHOICE_NNPDF30;
   Double_t        weight_kfactor_sys_PDF_EV1;
   Double_t        weight_kfactor_sys_PDF_EV2;
   Double_t        weight_kfactor_sys_PDF_EV3;
   Double_t        weight_kfactor_sys_PDF_EV4;
   Double_t        weight_kfactor_sys_PDF_EV5;
   Double_t        weight_kfactor_sys_PDF_EV6;
   Double_t        weight_kfactor_sys_PDF_EV7;
   Double_t        weight_kfactor_sys_PDF_EW__1down;
   Double_t        weight_kfactor_sys_PDF_EW__1up;
   Double_t        weight_kfactor_sys_PDF__1down;
   Double_t        weight_kfactor_sys_PDF__1up;
   Double_t        weight_kfactor_sys_PI__1down;
   Double_t        weight_kfactor_sys_PI__1up;
   Double_t        weight_kfactor_sys_REDCHOICE_NNPDF30;
   Double_t        weight_kfactor_sys_SCALE_W__1down;
   Double_t        weight_kfactor_sys_SCALE_W__1up;
   Double_t        weight_kfactor_sys_SCALE_Z__1down;
   Double_t        weight_kfactor_sys_SCALE_Z__1up;


   std::vector<float>   *el_pt;
   std::vector<float>   *el_eta;
   std::vector<float>   *el_cl_eta; //new
   std::vector<float>   *el_phi;
   std::vector<float>   *el_e;
   std::vector<float>   *el_charge;
   std::vector<float>   *el_topoetcone20;
   std::vector<float>   *el_ptvarcone20;
   //new variables
   std::vector<float>   *el_d0sig;
   std::vector<float>   *el_delta_z0_sintheta;
   std::vector<int>     *el_true_type;
   std::vector<int>     *el_true_origin;
   std::vector<int>     *el_true_typebkg;
   std::vector<int>     *el_true_originbkg;
   Float_t         hadronic_pt;
   Float_t         hadronic_eta;
   Float_t         hadronic_phi;
   Float_t         hadronic_e;
   Double_t        sumET_PFO;
   std::vector<float>   *mu_pt;
   std::vector<float>   *mu_eta;
   std::vector<float>   *mu_phi;
   std::vector<float>   *mu_e;
   std::vector<float>   *mu_charge;
   //new variables
   std::vector<float>   *mu_topoetcone20;
   std::vector<float>   *mu_topoetcone30;
   std::vector<float>   *mu_topoetcone40;
   std::vector<float>   *mu_ptvarcone20;
   std::vector<float>   *mu_ptvarcone30;
   std::vector<float>   *mu_ptvarcone40;
   std::vector<float>   *mu_ptcone20;
   std::vector<float>   *mu_ptcone30;
   std::vector<float>   *mu_ptcone40;
   std::vector<int>     *muon_gradientIsolation;
   std::vector<int>     *muon_trigger_mu20_iloose;
   std::vector<int>     *muon_trigger_mu50;
   std::vector<float>   *mu_d0sig;
   std::vector<float>   *mu_delta_z0_sintheta;
   std::vector<int>     *mu_true_type;
   std::vector<int>     *mu_true_origin;
   std::vector<float>   *jet_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_e;
   std::vector<float>   *jet_mv2c00;
   std::vector<float>   *jet_mv2c10;
   std::vector<float>   *jet_mv2c20;
   std::vector<float>   *jet_ip3dsv1;
   std::vector<float>   *jet_jvt;
   std::vector<int>     *jet_truthflav; //new
   Float_t         met_met;
   Float_t         met_phi;

   Int_t           munu_2015;
   Int_t           mumu_2015;

   Int_t           munu_2016;
   Int_t           mumu_2016;

   Int_t           munu_2017;
   Int_t           mumu_2017;

   Int_t           munu_2018;
   Int_t           mumu_2018;

   Char_t          HLT_mu50;
   Char_t          HLT_mu20_iloose_L1MU15;
   std::vector<char>    *mu_trigMatch_HLT_mu50;
   std::vector<char>    *mu_trigMatch_HLT_mu20_iloose_L1MU15;
   Int_t           primaryVertices;

   // List of branches
   TBranch        *b_weight_mc;   //!
   TBranch        *b_weight_pileup;   //!
   TBranch        *b_weight_leptonSF;   //!
   TBranch        *b_weight_globalLeptonTriggerSF;   //!
   TBranch        *b_weight_oldTriggerSF;   //!


   TBranch        *b_weight_KFactor;   //!
   TBranch        *b_bornMass_KFactor;   //!
   TBranch        *b_weight_kfactor_sys_ALPHAS__1down;   //!
   TBranch        *b_weight_kfactor_sys_ALPHAS__1up;   //!
   TBranch        *b_weight_kfactor_sys_BEAM_ENERGY__1down;   //!
   TBranch        *b_weight_kfactor_sys_BEAM_ENERGY__1up;   //!
   TBranch        *b_weight_kfactor_sys_CHOICE_HERAPDF20;   //!
   TBranch        *b_weight_kfactor_sys_CHOICE_NNPDF30;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV1;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV2;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV3;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV4;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV5;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV6;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EV7;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EW__1down;   //!
   TBranch        *b_weight_kfactor_sys_PDF_EW__1up;   //!
   TBranch        *b_weight_kfactor_sys_PDF__1down;   //!
   TBranch        *b_weight_kfactor_sys_PDF__1up;   //!
   TBranch        *b_weight_kfactor_sys_PI__1down;   //!
   TBranch        *b_weight_kfactor_sys_PI__1up;   //!
   TBranch        *b_weight_kfactor_sys_REDCHOICE_NNPDF30;   //!
   TBranch        *b_weight_kfactor_sys_SCALE_W__1down;   //!
   TBranch        *b_weight_kfactor_sys_SCALE_W__1up;   //!
   TBranch        *b_weight_kfactor_sys_SCALE_Z__1down;   //!
   TBranch        *b_weight_kfactor_sys_SCALE_Z__1up;   //!

   TBranch        *b_weight_jvt;   //!
   TBranch        *b_weight_pileup_UP;   //!
   TBranch        *b_weight_pileup_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Trigger_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Trigger_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Reco_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Reco_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_ID_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_ID_DOWN;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Isol_UP;   //!
   TBranch        *b_weight_leptonSF_EL_SF_Isol_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Trigger_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_Isol_SYST_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_STAT_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_STAT_DOWN;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_SYST_UP;   //!
   TBranch        *b_weight_leptonSF_MU_SF_TTVA_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Trigger_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ID;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_DOWN;   //!

   TBranch        *b_weight_globalLeptonTriggerSF_EL_Trigger_UP;   //!
   TBranch        *b_weight_globalLeptonTriggerSF_EL_Trigger_DOWN;   //!
   TBranch        *b_weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP;   //!
   TBranch        *b_weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN;   //!
   TBranch        *b_weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP;   //!
   TBranch        *b_weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN;   //!
   TBranch        *b_weight_oldTriggerSF_EL_Trigger_UP;   //!
   TBranch        *b_weight_oldTriggerSF_EL_Trigger_DOWN;   //!
   TBranch        *b_weight_oldTriggerSF_MU_Trigger_STAT_UP;   //!
   TBranch        *b_weight_oldTriggerSF_MU_Trigger_STAT_DOWN;   //!
   TBranch        *b_weight_oldTriggerSF_MU_Trigger_SYST_UP;   //!
   TBranch        *b_weight_oldTriggerSF_MU_Trigger_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_LOWPT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_LOWPT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_Isol_SYST_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_MU_TTVA_SYST_DOWN;   //!
   TBranch        *b_weight_jvt_UP;   //!
   TBranch        *b_weight_jvt_DOWN;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_randomRunNumber;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_backgroundFlags;   //!
   TBranch        *b_hasBadMuon;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_cl_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_topoetcone20;   //!
   TBranch        *b_el_ptvarcone20;   //!
   TBranch        *b_el_d0sig;   //!
   TBranch        *b_el_delta_z0_sintheta;   //!
   TBranch        *b_el_true_type;   //!
   TBranch        *b_el_true_origin;   //!
   TBranch        *b_el_true_typebkg;   //!
   TBranch        *b_el_true_originbkg;   //!
   TBranch        *b_hadronic_pt;   //!
   TBranch        *b_hadronic_eta;   //!
   TBranch        *b_hadronic_phi;   //!
   TBranch        *b_hadronic_e;   //!
   TBranch        *b_sumET_PFO;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_topoetcone20;   //!
   TBranch        *b_mu_topoetcone30;   //!
   TBranch        *b_mu_topoetcone40;   //!
   TBranch        *b_mu_ptvarcone20;   //!
   TBranch        *b_mu_ptvarcone30;   //!
   TBranch        *b_mu_ptvarcone40;   //!
   TBranch        *b_mu_ptcone20;   //!
   TBranch        *b_mu_ptcone30;   //!
   TBranch        *b_mu_ptcone40;   //!
   TBranch        *b_muon_gradientIsolation;   //!
   TBranch        *b_muon_trigger_mu20_iloose;   //!
   TBranch        *b_muon_trigger_mu50;   //!
   TBranch        *b_mu_d0sig;   //!
   TBranch        *b_mu_delta_z0_sintheta;   //!
   TBranch        *b_mu_true_type;   //!
   TBranch        *b_mu_true_origin;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_mv2c00;   //!
   TBranch        *b_jet_mv2c10;   //!
   TBranch        *b_jet_mv2c20;   //!
   TBranch        *b_jet_ip3dsv1;   //!
   TBranch        *b_jet_jvt;   //!
   TBranch        *b_jet_truthflav;   //!
   TBranch        *b_met_met;   //!
   TBranch        *b_met_phi;   //!

   TBranch        *b_munu_2015;   //!
   TBranch        *b_mumu_2015;   //!

   TBranch        *b_munu_2016;   //!
   TBranch        *b_mumu_2016;   //!

   TBranch        *b_munu_2017;   //!
   TBranch        *b_mumu_2017;   //!

   TBranch        *b_munu_2018;   //!
   TBranch        *b_mumu_2018;   //!

   TBranch        *b_HLT_mu50;   //!
   TBranch        *b_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_mu_trigMatch_HLT_mu50;   //!
   TBranch        *b_mu_trigMatch_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_primaryVertices;   //!
   TBranch        *b_weight_pileup_hash; //!

   MyWZAnalysis(Config, std::string, std::string,
	   double, double, double, double,
	   double, double, double, double);
   virtual ~MyWZAnalysis();
   //virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, Config config, std::string systematic);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

 private:

   //configuration
   virtual void SetConfiguration(Config config);
   bool  DataYears2015p2016;
   bool  DataYears2017;
   bool  DataYears2018;
   bool  OnlyInclusiveAnalysis;
   bool  DoMultijetAnalysis;       
   bool  MJFitRegionSignal;	       
   bool  MJFitRegionFitR1;	       
   bool  MJFitRegionFitR2;	       
   bool  Dod0CutAnalysis;	       
   bool  Applyd0WeightAnalysis;    
   bool  BosonZSelection;	       
   bool  BosonWplusSelection;      
   bool  BosonWminusSelection;     
   bool  HasRecoilInfoAnalysis;    
   bool  OnTheFlyPileUpAnalysis;              
   bool  SETCalibrationAnalysis;  
   bool  InsituCorrectionAnalysis;
   bool  ResolResponseAnalysis;   
   bool  TruthMatchingAnalysis;               
   bool  SFVariationsAnalysis;    
   
   bool basicZEvent, basicWEvent;
   bool goodZEvent, goodWEvent;
   bool goodZEventHigh, goodWEventHigh;
   bool goodWEventRecoil, goodWEventRecoilHigh;
   TLorentzVector Zboson, Zboson1, Zboson2;
   TLorentzVector Wboson; double wtm=0.;
   TLorentzVector tboson, upfo, upfo_high;
   bool d0sigCut;
   
   double pt=0., lpt=0., calo20=0., track20=0., eta=0., phi=0., dPhi=0.,  z0=0.;
   double pt2=0.,  calo202=0., track202=0., eta2=0., phi2=0., dPhi2=0., z02=0.;
   double iso30=0., d0sig=0., d0sig2=0.;
   double met=0., met_recoil=0., phimet=0., lrmass=0., recoMass=0., recoMass_recoil=0.;
   double pileup = 1., skimming = 1., kFactor=1., mc_weight=1.;
   double SF=1., SF_isol=1., SF_trig=1., SF_ID=1., SF_ttva=1.;
   double d0sigWeight=1.;
   float mu_scaled=0.;
   double px=0., py=0.;
   double px_met=0., py_met=0.;
   double px_sum=0., py_sum=0.;
   double pT_tot=0.;
   double finalWeight=1.;

   //d0sig reweighting
   double shift = 0;
   virtual void Getd0sigFit(Config config);
   TF1 *fit_d0sig_gaussratio;

   //histograms and profiles
   virtual void setHist(Config config, std::string sf_name);
   virtual void setHistMultijet(Config config, std::string eta_bin);
   virtual void WriteHist(int k);
   virtual void WriteHistMultijet(int j);
   virtual void FillHistos(int iSF_variations, double var);
   virtual void FillHistosMultijet(int iEta_bins, double eta_low, double eta_high);
   virtual void LumiNormalisation(Config config, int it_sf_eta, bool isData,
				  double lumi_sample, double tt2p);

   std::vector<TFile*> newfile;

   std::vector<TH1I*> h_npvtx, h_npvtx_high;
   std::vector<TH2D*> PtRmass, PtRmass_high;
   std::vector<TH1D*> h_rmassLogCw, h_rmassLogCw_high, h_mu, h_mu_high, h_pt20Track, h_pt20Track_high,h_SumET;
   std::vector<TH1D*> h_SumET_high, h_SumET_eta1, h_SumET_eta2, h_SumET_eta3, h_SumET_eta4, h_SumET_eta5,h_pt20Calo;
   std::vector<TH1D*> h_pt20Calo_high, h_pt, h_pt_high;
		
   std::vector<TH1D*> h_eta, h_eta_high, h_eta_d, h_eta_high_d, h_phi, h_phi_high, h_phi_met, h_phi_met_high;	
   std::vector<TH1D*> h_delta_phi, h_delta_phi_high, h_met, h_met_high, h_met_recoil, h_met_recoil_high;
	
   std::vector<TH1D*> h_rmass, h_rmass_high, h_rmass_recoil, h_rmass_recoil_high, h_d0, h_d0_high, h_z0, h_z0_high;	
   std::vector<TH1D*> h_pt_boson_reco, h_pt_boson_reco_p, h_pt_boson, h_mass_boson, h_eta_boson, h_pt_boson_high, h_mass_boson_high, h_eta_boson_high;

   //no cut histos	
   TH1D *h_truth_mass_nocuts;
   TH1D *hPileUp;
   TH1I *h_nm;

   //Cw factor calculation
   TH1D* hMassCrec ;
   TH1D* hMassCgen ;
   TH1D* hMassCstay;
   TH1D* hMassCcome;
   
   TH1D* hEtaCrec ;
   TH1D* hEtaCgen ;
   TH1D* hEtaCstay;
   TH1D* hEtaCcome;

   //Cutflow
   TH1D *h_cutflow, *h_cutflow_high;
   
   //TProfiles
   TProfile *p_mwt_Tot ;
   TProfile *p_mwt_Rec ;
   TProfile *p_mwt_KF  ;
   TProfile *p_mwt_MC  ;
  
   TProfile *p_mwt_Iso           ;
   TProfile *p_mwt_Iso_sys_UP    ;
   TProfile *p_mwt_Iso_sys_DOWN  ;
   TProfile *p_mwt_Iso_stat_UP   ;
   TProfile *p_mwt_Iso_stat_DOWN ;
   
   TProfile *p_mwt_ID           ;
   TProfile *p_mwt_ID_sys_UP    ;
   TProfile *p_mwt_ID_sys_DOWN  ;
   TProfile *p_mwt_ID_stat_UP   ;
   TProfile *p_mwt_ID_stat_DOWN ;
   
   TProfile *p_mwt_Tri          ;
   TProfile *p_mwt_Tri_sys_UP   ;
   TProfile *p_mwt_Tri_sys_DOWN ;
   TProfile *p_mwt_Tri_stat_UP  ;
   TProfile *p_mwt_Tri_stat_DOWN;
   
   TProfile *p_mwt_ttv          ;
   TProfile *p_mwt_ttv_sys_UP   ;
   TProfile *p_mwt_ttv_sys_DOWN ;
   TProfile *p_mwt_ttv_stat_UP  ;
   TProfile *p_mwt_ttv_stat_DOWN;
   
   TProfile *p_mwt_PU     ;
   TProfile *p_mwt_PU_UP  ;
   TProfile *p_mwt_PU_DOWN;
   
   TProfile *p_mwt_jvt      ;
   TProfile *p_mwt_jvt_UP   ; 
   TProfile *p_mwt_jvt_DOWN ;   
   
   TProfile *p_m_Tot ;
   TProfile *p_m_Rec ;
   TProfile *p_m_KF  ;
   TProfile *p_m_MC  ;
   
   TProfile *p_m_Iso           ;
   TProfile *p_m_Iso_sys_UP    ;
   TProfile *p_m_Iso_sys_DOWN  ;
   TProfile *p_m_Iso_stat_UP   ;
   TProfile *p_m_Iso_stat_DOWN ;
   
   TProfile *p_m_ID          ;
   TProfile *p_m_ID_sys_UP   ;
   TProfile *p_m_ID_sys_DOWN ;
   TProfile *p_m_ID_stat_UP  ;
   TProfile *p_m_ID_stat_DOWN;
   
   TProfile *p_m_Tri          ;
   TProfile *p_m_Tri_sys_UP   ;
   TProfile *p_m_Tri_sys_DOWN ;
   TProfile *p_m_Tri_stat_UP  ;
   TProfile *p_m_Tri_stat_DOWN;
   
   TProfile *p_m_ttv          ;
   TProfile *p_m_ttv_sys_UP   ;
   TProfile *p_m_ttv_sys_DOWN ;
   TProfile *p_m_ttv_stat_UP  ;
   TProfile *p_m_ttv_stat_DOWN;
   
   TProfile *p_m_PU     ; 
   TProfile *p_m_PU_UP  ; 
   TProfile *p_m_PU_DOWN; 
   
   TProfile *p_m_jvt     ;
   TProfile *p_m_jvt_UP  ; 
   TProfile *p_m_jvt_DOWN;   
   
   TProfile *p_pt_Tot ;
   TProfile *p_pt_Rec ;
   TProfile *p_pt_KF  ;
   TProfile *p_pt_MC  ;
   
   TProfile *p_pt_Iso          ;
   TProfile *p_pt_Iso_sys_UP   ;
   TProfile *p_pt_Iso_sys_DOWN ;
   TProfile *p_pt_Iso_stat_UP  ;
   TProfile *p_pt_Iso_stat_DOWN;
   
   TProfile *p_pt_ID          ;
   TProfile *p_pt_ID_sys_UP   ;
   TProfile *p_pt_ID_sys_DOWN ;
   TProfile *p_pt_ID_stat_UP  ;
   TProfile *p_pt_ID_stat_DOWN;
   
   TProfile *p_pt_Tri          ;
   TProfile *p_pt_Tri_sys_UP   ;
   TProfile *p_pt_Tri_sys_DOWN ;
   TProfile *p_pt_Tri_stat_UP  ;
   TProfile *p_pt_Tri_stat_DOWN;
   
   TProfile *p_pt_ttv          ;
   TProfile *p_pt_ttv_sys_UP   ;
   TProfile *p_pt_ttv_sys_DOWN ;
   TProfile *p_pt_ttv_stat_UP  ;
   TProfile *p_pt_ttv_stat_DOWN;
   
   TProfile *p_pt_PU       ;
   TProfile *p_pt_PU_UP    ;
   TProfile *p_pt_PU_DOWN  ;
   
   TProfile *p_pt_jvt      ;
   TProfile *p_pt_jvt_UP   ; 
   TProfile *p_pt_jvt_DOWN ;   
   
   TProfile *p_eta_Tot ;
   TProfile *p_eta_Rec ;
   TProfile *p_eta_KF  ;
   TProfile *p_eta_MC  ;
   
   TProfile *p_eta_Iso          ;
   TProfile *p_eta_Iso_sys_UP   ;
   TProfile *p_eta_Iso_sys_DOWN ;
   TProfile *p_eta_Iso_stat_UP  ;
   TProfile *p_eta_Iso_stat_DOWN;
   
   TProfile *p_eta_ID           ;
   TProfile *p_eta_ID_sys_UP    ;
   TProfile *p_eta_ID_sys_DOWN  ;
   TProfile *p_eta_ID_stat_UP   ;
   TProfile *p_eta_ID_stat_DOWN ;
   
   TProfile *p_eta_Tri          ;
   TProfile *p_eta_Tri_sys_UP   ;
   TProfile *p_eta_Tri_sys_DOWN ;
   TProfile *p_eta_Tri_stat_UP  ;
   TProfile *p_eta_Tri_stat_DOWN;
   
   TProfile *p_eta_ttv          ;
   TProfile *p_eta_ttv_sys_UP   ;
   TProfile *p_eta_ttv_sys_DOWN ;
   TProfile *p_eta_ttv_stat_UP  ;
   TProfile *p_eta_ttv_stat_DOWN;
   
   TProfile *p_eta_PU       ;
   TProfile *p_eta_PU_UP    ;
   TProfile *p_eta_PU_DOWN  ;
   
   TProfile *p_eta_jvt      ;
   TProfile *p_eta_jvt_UP   ; 
   TProfile *p_eta_jvt_DOWN ;   
      
   //multijet background   
   std::vector<TH1D*> h_rmass_i; 
   std::vector<TH1D*> h_rmass_0; 
   std::vector<TH1D*> h_rmass_1; 
   std::vector<TH1D*> h_rmass_2; 
   std::vector<TH1D*> h_rmass_3; 
   std::vector<TH1D*> h_rmass_4; 
   std::vector<TH1D*> h_rmass_5;
   
   std::vector<TH1D*> h_pt_i;
   std::vector<TH1D*> h_pt_0;
   std::vector<TH1D*> h_pt_1;
   std::vector<TH1D*> h_pt_2;
   std::vector<TH1D*> h_pt_3;
   std::vector<TH1D*> h_pt_4;
   std::vector<TH1D*> h_pt_5;
   
   std::vector<TH1D*> h_eta_i; 
   std::vector<TH1D*> h_eta_0; 
   std::vector<TH1D*> h_eta_1; 
   std::vector<TH1D*> h_eta_2; 
   std::vector<TH1D*> h_eta_3; 
   std::vector<TH1D*> h_eta_4; 
   std::vector<TH1D*> h_eta_5; 
   
   std::vector<TH1D*> h_met_i;
   std::vector<TH1D*> h_met_0;
   std::vector<TH1D*> h_met_1;
   std::vector<TH1D*> h_met_2; 
   std::vector<TH1D*> h_met_3; 
   std::vector<TH1D*> h_met_4; 
   std::vector<TH1D*> h_met_5; 
   
   std::vector<TH1D*> h_iso30_i; 
   std::vector<TH1D*> h_iso30_0; 
   std::vector<TH1D*> h_iso30_1; 
   std::vector<TH1D*> h_iso30_2; 
   std::vector<TH1D*> h_iso30_3; 
   std::vector<TH1D*> h_iso30_4; 
   std::vector<TH1D*> h_iso30_5;

   std::vector<TH1D*> h_d0sig_i; 
   std::vector<TH1D*> h_d0sig_0; 
   std::vector<TH1D*> h_d0sig_1; 
   std::vector<TH1D*> h_d0sig_2; 
   std::vector<TH1D*> h_d0sig_3; 
   std::vector<TH1D*> h_d0sig_4; 
   std::vector<TH1D*> h_d0sig_5;

   std::vector<TH1D*> h_ptDmt_i; 
   std::vector<TH1D*> h_ptDmt_0; 
   std::vector<TH1D*> h_ptDmt_1; 
   std::vector<TH1D*> h_ptDmt_2; 
   std::vector<TH1D*> h_ptDmt_3; 
   std::vector<TH1D*> h_ptDmt_4; 
   std::vector<TH1D*> h_ptDmt_5;

   std::vector<TH1D*> h_eta_d_i;	    
   std::vector<TH1D*> h_eta_d_0;	    
   std::vector<TH1D*> h_eta_d_1;	    
   std::vector<TH1D*> h_eta_d_2;	    
   std::vector<TH1D*> h_eta_d_3;	    
   std::vector<TH1D*> h_eta_d_4;	    
   std::vector<TH1D*> h_eta_d_5;	    
   
   std::vector<TH1D*> h_phi_i;
   std::vector<TH1D*> h_phi_0;
   std::vector<TH1D*> h_phi_1;
   std::vector<TH1D*> h_phi_2;
   std::vector<TH1D*> h_phi_3;
   std::vector<TH1D*> h_phi_4;
   std::vector<TH1D*> h_phi_5;
   
   std::vector<TH1D*> h_phi_met_i;
   std::vector<TH1D*> h_phi_met_0;
   std::vector<TH1D*> h_phi_met_1;
   std::vector<TH1D*> h_phi_met_2;
   std::vector<TH1D*> h_phi_met_3;
   std::vector<TH1D*> h_phi_met_4;
   std::vector<TH1D*> h_phi_met_5;
   
   std::vector<TH1D*> h_delta_phi_i;
   std::vector<TH1D*> h_delta_phi_0;
   std::vector<TH1D*> h_delta_phi_1;
   std::vector<TH1D*> h_delta_phi_2;
   std::vector<TH1D*> h_delta_phi_3;
   std::vector<TH1D*> h_delta_phi_4;
   std::vector<TH1D*> h_delta_phi_5;
   
   std::vector<TH1D*> h_mu_i;
   std::vector<TH1D*> h_mu_0;
   std::vector<TH1D*> h_mu_1;
   std::vector<TH1D*> h_mu_2;
   std::vector<TH1D*> h_mu_3;
   std::vector<TH1D*> h_mu_4;
   std::vector<TH1D*> h_mu_5;
};

#endif

#ifdef MyWZAnalysis_cxx
MyWZAnalysis::MyWZAnalysis(Config myconfig, std::string samples, std::string systematics, 
			   double luminosity_samples,
			   double eta_cut, double pt_cut1, double pt_cut2, double met_cut1,
			   double met_cut2, double mwt_cut1, double mwt_cut2) : fChain(0),fChain_sumWeights(0),fChain_truth(0)
			   {
  //andres
  config = myconfig;
  nameOfSample = samples;
  nameOfSystematic = systematics;
  lumi_sample = luminosity_samples;
  cut_eta = eta_cut;
  cut_pt1 = pt_cut1;
  cut_pt2 = pt_cut2;
  cut_met1 = met_cut1;
  cut_met2 = met_cut2;
  cut_mwt1 = mwt_cut1;
  cut_mwt2 = mwt_cut2;

  TTree* tree;
  /* TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject((config.InputFileDir+nameOfSample+"/\*.root").c_str()); */
  /*     if (!f || !f->IsOpen()) { */
  /* 	f = new TFile((config.InputFileDir+nameOfSample+"/\*.root").c_str()); */
  /*     } */
  /*     f->GetObject(nameOfSystematic.c_str(),tree); */

  //Reco analysis
  TChain *chain = new TChain(nameOfSystematic.c_str(),"");
  chain->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
  tree=chain;
  Init(tree,config,nameOfSystematic);

  bool isData =false;
  if(nameOfSample=="data16" || nameOfSample=="data15" || nameOfSample=="data17" || nameOfSample=="data18") isData=true;

  if(!isData){
    //sum of weights
    TChain *chain_sumWeights = new TChain("sumWeights","");
    chain_sumWeights->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
    fChain_sumWeights=chain_sumWeights;
    
    //truth events, maaaps
    if(config.OnlyInclusive!="True" && config.TruthMatching=="True"){
      TChain *chain_truth = new TChain("truth","");
      chain_truth->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
      fChain_truth=chain_truth;
    }
  }//isData

}


MyWZAnalysis::~MyWZAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   if(!fChain_sumWeights) return;
   delete fChain_sumWeights->GetCurrentFile();
   if(!fChain_truth) return;
   delete fChain_truth->GetCurrentFile();
}

Int_t MyWZAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyWZAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyWZAnalysis::Init(TTree *tree, Config config, std::string nameOfSystematic)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   el_pt = 0;
   el_eta = 0;
   el_cl_eta = 0;
   el_phi = 0;
   el_e = 0;
   el_charge = 0;
   el_topoetcone20 = 0;
   el_ptvarcone20 = 0;
   el_d0sig = 0;
   el_delta_z0_sintheta = 0;
   el_true_type = 0;
   el_true_origin = 0;
   el_true_typebkg = 0;
   el_true_originbkg = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_e = 0;
   mu_charge = 0;
   mu_topoetcone20 = 0;
   mu_topoetcone30 = 0;
   mu_topoetcone40 = 0;
   mu_ptvarcone20 = 0;
   mu_ptvarcone30 = 0;
   mu_ptvarcone40 = 0;
   mu_ptcone20 = 0;
   mu_ptcone30 = 0;
   mu_ptcone40 = 0;
   muon_gradientIsolation = 0;
   muon_trigger_mu20_iloose = 0;
   muon_trigger_mu50 = 0;
   mu_d0sig = 0;
   mu_delta_z0_sintheta = 0;
   mu_true_type = 0;
   mu_true_origin = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_e = 0;
   jet_mv2c00 = 0;
   jet_mv2c10 = 0;
   jet_mv2c20 = 0;
   jet_ip3dsv1 = 0;
   jet_jvt = 0;
   jet_truthflav = 0;
   mu_trigMatch_HLT_mu50 = 0;
   mu_trigMatch_HLT_mu20_iloose_L1MU15 = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   bool isData =false;
   if(nameOfSample=="data16" || nameOfSample=="data15" || nameOfSample=="data17" || nameOfSample=="data18") isData=true;

   if(!isData){
     fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
     fChain->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
     fChain->SetBranchAddress("weight_leptonSF", &weight_leptonSF, &b_weight_leptonSF);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF", &weight_globalLeptonTriggerSF, &b_weight_globalLeptonTriggerSF);
     fChain->SetBranchAddress("weight_oldTriggerSF", &weight_oldTriggerSF, &b_weight_oldTriggerSF);
     fChain->SetBranchAddress("weight_KFactor", &weight_KFactor, &b_weight_KFactor);
     fChain->SetBranchAddress("bornMass_KFactor", &bornMass_KFactor, &b_bornMass_KFactor);
     fChain->SetBranchAddress("weight_jvt", &weight_jvt, &b_weight_jvt);
   }
     
   if(nameOfSystematic=="nominal" && !isData){
     fChain->SetBranchAddress("weight_kfactor_sys_ALPHAS__1down", &weight_kfactor_sys_ALPHAS__1down, &b_weight_kfactor_sys_ALPHAS__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_ALPHAS__1up", &weight_kfactor_sys_ALPHAS__1up, &b_weight_kfactor_sys_ALPHAS__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_BEAM_ENERGY__1down", &weight_kfactor_sys_BEAM_ENERGY__1down, &b_weight_kfactor_sys_BEAM_ENERGY__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_BEAM_ENERGY__1up", &weight_kfactor_sys_BEAM_ENERGY__1up, &b_weight_kfactor_sys_BEAM_ENERGY__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_CHOICE_HERAPDF20", &weight_kfactor_sys_CHOICE_HERAPDF20, &b_weight_kfactor_sys_CHOICE_HERAPDF20);
     fChain->SetBranchAddress("weight_kfactor_sys_CHOICE_NNPDF30", &weight_kfactor_sys_CHOICE_NNPDF30, &b_weight_kfactor_sys_CHOICE_NNPDF30);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV1", &weight_kfactor_sys_PDF_EV1, &b_weight_kfactor_sys_PDF_EV1);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV2", &weight_kfactor_sys_PDF_EV2, &b_weight_kfactor_sys_PDF_EV2);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV3", &weight_kfactor_sys_PDF_EV3, &b_weight_kfactor_sys_PDF_EV3);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV4", &weight_kfactor_sys_PDF_EV4, &b_weight_kfactor_sys_PDF_EV4);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV5", &weight_kfactor_sys_PDF_EV5, &b_weight_kfactor_sys_PDF_EV5);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV6", &weight_kfactor_sys_PDF_EV6, &b_weight_kfactor_sys_PDF_EV6);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EV7", &weight_kfactor_sys_PDF_EV7, &b_weight_kfactor_sys_PDF_EV7);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EW__1down", &weight_kfactor_sys_PDF_EW__1down, &b_weight_kfactor_sys_PDF_EW__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF_EW__1up", &weight_kfactor_sys_PDF_EW__1up, &b_weight_kfactor_sys_PDF_EW__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF__1down", &weight_kfactor_sys_PDF__1down, &b_weight_kfactor_sys_PDF__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_PDF__1up", &weight_kfactor_sys_PDF__1up, &b_weight_kfactor_sys_PDF__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_PI__1down", &weight_kfactor_sys_PI__1down, &b_weight_kfactor_sys_PI__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_PI__1up", &weight_kfactor_sys_PI__1up, &b_weight_kfactor_sys_PI__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_REDCHOICE_NNPDF30", &weight_kfactor_sys_REDCHOICE_NNPDF30, &b_weight_kfactor_sys_REDCHOICE_NNPDF30);
     fChain->SetBranchAddress("weight_kfactor_sys_SCALE_W__1down", &weight_kfactor_sys_SCALE_W__1down, &b_weight_kfactor_sys_SCALE_W__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_SCALE_W__1up", &weight_kfactor_sys_SCALE_W__1up, &b_weight_kfactor_sys_SCALE_W__1up);
     fChain->SetBranchAddress("weight_kfactor_sys_SCALE_Z__1down", &weight_kfactor_sys_SCALE_Z__1down, &b_weight_kfactor_sys_SCALE_Z__1down);
     fChain->SetBranchAddress("weight_kfactor_sys_SCALE_Z__1up", &weight_kfactor_sys_SCALE_Z__1up, &b_weight_kfactor_sys_SCALE_Z__1up);
     fChain->SetBranchAddress("weight_pileup_UP", &weight_pileup_UP, &b_weight_pileup_UP);
     fChain->SetBranchAddress("weight_pileup_DOWN", &weight_pileup_DOWN, &b_weight_pileup_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Trigger_UP", &weight_leptonSF_EL_SF_Trigger_UP, &b_weight_leptonSF_EL_SF_Trigger_UP);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Trigger_DOWN", &weight_leptonSF_EL_SF_Trigger_DOWN, &b_weight_leptonSF_EL_SF_Trigger_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Reco_UP", &weight_leptonSF_EL_SF_Reco_UP, &b_weight_leptonSF_EL_SF_Reco_UP);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Reco_DOWN", &weight_leptonSF_EL_SF_Reco_DOWN, &b_weight_leptonSF_EL_SF_Reco_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_ID_UP", &weight_leptonSF_EL_SF_ID_UP, &b_weight_leptonSF_EL_SF_ID_UP);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_ID_DOWN", &weight_leptonSF_EL_SF_ID_DOWN, &b_weight_leptonSF_EL_SF_ID_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Isol_UP", &weight_leptonSF_EL_SF_Isol_UP, &b_weight_leptonSF_EL_SF_Isol_UP);
     fChain->SetBranchAddress("weight_leptonSF_EL_SF_Isol_DOWN", &weight_leptonSF_EL_SF_Isol_DOWN, &b_weight_leptonSF_EL_SF_Isol_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_STAT_UP", &weight_leptonSF_MU_SF_Trigger_STAT_UP, &b_weight_leptonSF_MU_SF_Trigger_STAT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_STAT_DOWN", &weight_leptonSF_MU_SF_Trigger_STAT_DOWN, &b_weight_leptonSF_MU_SF_Trigger_STAT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_SYST_UP", &weight_leptonSF_MU_SF_Trigger_SYST_UP, &b_weight_leptonSF_MU_SF_Trigger_SYST_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Trigger_SYST_DOWN", &weight_leptonSF_MU_SF_Trigger_SYST_DOWN, &b_weight_leptonSF_MU_SF_Trigger_SYST_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_UP", &weight_leptonSF_MU_SF_ID_STAT_UP, &b_weight_leptonSF_MU_SF_ID_STAT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_DOWN", &weight_leptonSF_MU_SF_ID_STAT_DOWN, &b_weight_leptonSF_MU_SF_ID_STAT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_UP", &weight_leptonSF_MU_SF_ID_SYST_UP, &b_weight_leptonSF_MU_SF_ID_SYST_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_DOWN", &weight_leptonSF_MU_SF_ID_SYST_DOWN, &b_weight_leptonSF_MU_SF_ID_SYST_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP", &weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP, &b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN", &weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN, &b_weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP", &weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP, &b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN", &weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN, &b_weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_STAT_UP", &weight_leptonSF_MU_SF_Isol_STAT_UP, &b_weight_leptonSF_MU_SF_Isol_STAT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_STAT_DOWN", &weight_leptonSF_MU_SF_Isol_STAT_DOWN, &b_weight_leptonSF_MU_SF_Isol_STAT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_SYST_UP", &weight_leptonSF_MU_SF_Isol_SYST_UP, &b_weight_leptonSF_MU_SF_Isol_SYST_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_Isol_SYST_DOWN", &weight_leptonSF_MU_SF_Isol_SYST_DOWN, &b_weight_leptonSF_MU_SF_Isol_SYST_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_STAT_UP", &weight_leptonSF_MU_SF_TTVA_STAT_UP, &b_weight_leptonSF_MU_SF_TTVA_STAT_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_STAT_DOWN", &weight_leptonSF_MU_SF_TTVA_STAT_DOWN, &b_weight_leptonSF_MU_SF_TTVA_STAT_DOWN);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_SYST_UP", &weight_leptonSF_MU_SF_TTVA_SYST_UP, &b_weight_leptonSF_MU_SF_TTVA_SYST_UP);
     fChain->SetBranchAddress("weight_leptonSF_MU_SF_TTVA_SYST_DOWN", &weight_leptonSF_MU_SF_TTVA_SYST_DOWN, &b_weight_leptonSF_MU_SF_TTVA_SYST_DOWN);
     if(false){//remove/fix this?
       fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger", &weight_indiv_SF_EL_Trigger, &b_weight_indiv_SF_EL_Trigger);
       fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger_UP", &weight_indiv_SF_EL_Trigger_UP, &b_weight_indiv_SF_EL_Trigger_UP);
       fChain->SetBranchAddress("weight_indiv_SF_EL_Trigger_DOWN", &weight_indiv_SF_EL_Trigger_DOWN, &b_weight_indiv_SF_EL_Trigger_DOWN);
     }
     fChain->SetBranchAddress("weight_indiv_SF_EL_Reco", &weight_indiv_SF_EL_Reco, &b_weight_indiv_SF_EL_Reco);
     fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_UP", &weight_indiv_SF_EL_Reco_UP, &b_weight_indiv_SF_EL_Reco_UP);
     fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_DOWN", &weight_indiv_SF_EL_Reco_DOWN, &b_weight_indiv_SF_EL_Reco_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_EL_ID", &weight_indiv_SF_EL_ID, &b_weight_indiv_SF_EL_ID);
     fChain->SetBranchAddress("weight_indiv_SF_EL_ID_UP", &weight_indiv_SF_EL_ID_UP, &b_weight_indiv_SF_EL_ID_UP);
     fChain->SetBranchAddress("weight_indiv_SF_EL_ID_DOWN", &weight_indiv_SF_EL_ID_DOWN, &b_weight_indiv_SF_EL_ID_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_EL_Isol", &weight_indiv_SF_EL_Isol, &b_weight_indiv_SF_EL_Isol);
     fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_UP", &weight_indiv_SF_EL_Isol_UP, &b_weight_indiv_SF_EL_Isol_UP);
     fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_DOWN", &weight_indiv_SF_EL_Isol_DOWN, &b_weight_indiv_SF_EL_Isol_DOWN);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_EL_Trigger_UP", &weight_globalLeptonTriggerSF_EL_Trigger_UP, &b_weight_globalLeptonTriggerSF_EL_Trigger_UP);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_EL_Trigger_DOWN", &weight_globalLeptonTriggerSF_EL_Trigger_DOWN, &b_weight_globalLeptonTriggerSF_EL_Trigger_DOWN);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP", &weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP, &b_weight_globalLeptonTriggerSF_MU_Trigger_STAT_UP);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN", &weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN, &b_weight_globalLeptonTriggerSF_MU_Trigger_STAT_DOWN);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP", &weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP, &b_weight_globalLeptonTriggerSF_MU_Trigger_SYST_UP);
     fChain->SetBranchAddress("weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN", &weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN, &b_weight_globalLeptonTriggerSF_MU_Trigger_SYST_DOWN);
     fChain->SetBranchAddress("weight_oldTriggerSF_EL_Trigger_UP", &weight_oldTriggerSF_EL_Trigger_UP, &b_weight_oldTriggerSF_EL_Trigger_UP);
     fChain->SetBranchAddress("weight_oldTriggerSF_EL_Trigger_DOWN", &weight_oldTriggerSF_EL_Trigger_DOWN, &b_weight_oldTriggerSF_EL_Trigger_DOWN);
     fChain->SetBranchAddress("weight_oldTriggerSF_MU_Trigger_STAT_UP", &weight_oldTriggerSF_MU_Trigger_STAT_UP, &b_weight_oldTriggerSF_MU_Trigger_STAT_UP);
     fChain->SetBranchAddress("weight_oldTriggerSF_MU_Trigger_STAT_DOWN", &weight_oldTriggerSF_MU_Trigger_STAT_DOWN, &b_weight_oldTriggerSF_MU_Trigger_STAT_DOWN);
     fChain->SetBranchAddress("weight_oldTriggerSF_MU_Trigger_SYST_UP", &weight_oldTriggerSF_MU_Trigger_SYST_UP, &b_weight_oldTriggerSF_MU_Trigger_SYST_UP);
     fChain->SetBranchAddress("weight_oldTriggerSF_MU_Trigger_SYST_DOWN", &weight_oldTriggerSF_MU_Trigger_SYST_DOWN, &b_weight_oldTriggerSF_MU_Trigger_SYST_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID", &weight_indiv_SF_MU_ID, &b_weight_indiv_SF_MU_ID);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_UP", &weight_indiv_SF_MU_ID_STAT_UP, &b_weight_indiv_SF_MU_ID_STAT_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_DOWN", &weight_indiv_SF_MU_ID_STAT_DOWN, &b_weight_indiv_SF_MU_ID_STAT_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_UP", &weight_indiv_SF_MU_ID_SYST_UP, &b_weight_indiv_SF_MU_ID_SYST_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_DOWN", &weight_indiv_SF_MU_ID_SYST_DOWN, &b_weight_indiv_SF_MU_ID_SYST_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_LOWPT_UP", &weight_indiv_SF_MU_ID_STAT_LOWPT_UP, &b_weight_indiv_SF_MU_ID_STAT_LOWPT_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN", &weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN, &b_weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_LOWPT_UP", &weight_indiv_SF_MU_ID_SYST_LOWPT_UP, &b_weight_indiv_SF_MU_ID_SYST_LOWPT_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN", &weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN, &b_weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_Isol", &weight_indiv_SF_MU_Isol, &b_weight_indiv_SF_MU_Isol);
     fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_STAT_UP", &weight_indiv_SF_MU_Isol_STAT_UP, &b_weight_indiv_SF_MU_Isol_STAT_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_STAT_DOWN", &weight_indiv_SF_MU_Isol_STAT_DOWN, &b_weight_indiv_SF_MU_Isol_STAT_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_SYST_UP", &weight_indiv_SF_MU_Isol_SYST_UP, &b_weight_indiv_SF_MU_Isol_SYST_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_Isol_SYST_DOWN", &weight_indiv_SF_MU_Isol_SYST_DOWN, &b_weight_indiv_SF_MU_Isol_SYST_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA", &weight_indiv_SF_MU_TTVA, &b_weight_indiv_SF_MU_TTVA);
     fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_STAT_UP", &weight_indiv_SF_MU_TTVA_STAT_UP, &b_weight_indiv_SF_MU_TTVA_STAT_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_STAT_DOWN", &weight_indiv_SF_MU_TTVA_STAT_DOWN, &b_weight_indiv_SF_MU_TTVA_STAT_DOWN);
     fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_SYST_UP", &weight_indiv_SF_MU_TTVA_SYST_UP, &b_weight_indiv_SF_MU_TTVA_SYST_UP);
     fChain->SetBranchAddress("weight_indiv_SF_MU_TTVA_SYST_DOWN", &weight_indiv_SF_MU_TTVA_SYST_DOWN, &b_weight_indiv_SF_MU_TTVA_SYST_DOWN);
     fChain->SetBranchAddress("weight_jvt_UP", &weight_jvt_UP, &b_weight_jvt_UP);
     fChain->SetBranchAddress("weight_jvt_DOWN", &weight_jvt_DOWN, &b_weight_jvt_DOWN);
   }

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   if(!isData) fChain->SetBranchAddress("randomRunNumber", &randomRunNumber, &b_randomRunNumber);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("backgroundFlags", &backgroundFlags, &b_backgroundFlags);
   fChain->SetBranchAddress("hasBadMuon", &hasBadMuon, &b_hasBadMuon);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_cl_eta", &el_cl_eta, &b_el_cl_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_e", &el_e, &b_el_e);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_topoetcone20", &el_topoetcone20, &b_el_topoetcone20);
   fChain->SetBranchAddress("el_ptvarcone20", &el_ptvarcone20, &b_el_ptvarcone20);
   fChain->SetBranchAddress("el_d0sig", &el_d0sig, &b_el_d0sig);
   fChain->SetBranchAddress("el_delta_z0_sintheta", &el_delta_z0_sintheta, &b_el_delta_z0_sintheta);
   if(false){//remove/fix this?
     fChain->SetBranchAddress("el_true_type", &el_true_type, &b_el_true_type);
     fChain->SetBranchAddress("el_true_origin", &el_true_origin, &b_el_true_origin);
     fChain->SetBranchAddress("el_true_typebkg", &el_true_typebkg, &b_el_true_typebkg);
     fChain->SetBranchAddress("el_true_originbkg", &el_true_originbkg, &b_el_true_originbkg);
   }
   fChain->SetBranchAddress("hadronic_pt", &hadronic_pt, &b_hadronic_pt);
   fChain->SetBranchAddress("hadronic_eta", &hadronic_eta, &b_hadronic_eta);
   fChain->SetBranchAddress("hadronic_phi", &hadronic_phi, &b_hadronic_phi);
   fChain->SetBranchAddress("hadronic_e", &hadronic_e, &b_hadronic_e);
   fChain->SetBranchAddress("sumET_PFO", &sumET_PFO, &b_sumET_PFO);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_topoetcone20", &mu_topoetcone20, &b_mu_topoetcone20);
   fChain->SetBranchAddress("mu_topoetcone30", &mu_topoetcone30, &b_mu_topoetcone30);
   fChain->SetBranchAddress("mu_topoetcone40", &mu_topoetcone40, &b_mu_topoetcone40);
   fChain->SetBranchAddress("mu_ptvarcone20", &mu_ptvarcone20, &b_mu_ptvarcone20);
   fChain->SetBranchAddress("mu_ptvarcone30", &mu_ptvarcone30, &b_mu_ptvarcone30);
   fChain->SetBranchAddress("mu_ptvarcone40", &mu_ptvarcone40, &b_mu_ptvarcone40);
   fChain->SetBranchAddress("mu_ptcone20", &mu_ptcone20, &b_mu_ptcone20);
   fChain->SetBranchAddress("mu_ptcone30", &mu_ptcone30, &b_mu_ptcone30);
   fChain->SetBranchAddress("mu_ptcone40", &mu_ptcone40, &b_mu_ptcone40);
   fChain->SetBranchAddress("muon_gradientIsolation", &muon_gradientIsolation, &b_muon_gradientIsolation);
   if(false)//remove fix this
     fChain->SetBranchAddress("muon_trigger_mu20_iloose", &muon_trigger_mu20_iloose, &b_muon_trigger_mu20_iloose);
   fChain->SetBranchAddress("muon_trigger_mu50", &muon_trigger_mu50, &b_muon_trigger_mu50);
   fChain->SetBranchAddress("mu_d0sig", &mu_d0sig, &b_mu_d0sig);
   fChain->SetBranchAddress("mu_delta_z0_sintheta", &mu_delta_z0_sintheta, &b_mu_delta_z0_sintheta);
   if(!isData){
     fChain->SetBranchAddress("mu_true_type", &mu_true_type, &b_mu_true_type);
     fChain->SetBranchAddress("mu_true_origin", &mu_true_origin, &b_mu_true_origin);}
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_mv2c00", &jet_mv2c00, &b_jet_mv2c00);
   fChain->SetBranchAddress("jet_mv2c10", &jet_mv2c10, &b_jet_mv2c10);
   fChain->SetBranchAddress("jet_mv2c20", &jet_mv2c20, &b_jet_mv2c20);
   fChain->SetBranchAddress("jet_ip3dsv1", &jet_ip3dsv1, &b_jet_ip3dsv1);
   fChain->SetBranchAddress("jet_jvt", &jet_jvt, &b_jet_jvt);
   if(!isData) fChain->SetBranchAddress("jet_truthflav", &jet_truthflav, &b_jet_truthflav);
   fChain->SetBranchAddress("met_met", &met_met, &b_met_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);

   if(config.DataYears=="2015+2016"){
     if(config.WZSelection=="wminus" || config.WZSelection=="wplus"){
       fChain->SetBranchAddress("munu_2015", &munu_2015, &b_munu_2015);
       fChain->SetBranchAddress("munu_2016", &munu_2016, &b_munu_2016);
     }else if(config.WZSelection=="zmumu"){
       fChain->SetBranchAddress("mumu_2015", &mumu_2015, &b_mumu_2015);   
       fChain->SetBranchAddress("mumu_2016", &mumu_2016, &b_mumu_2016);
     }
   }

   if(config.DataYears=="2017"){
     if(config.WZSelection=="wminus" || config.WZSelection=="wplus")
       fChain->SetBranchAddress("munu_2017", &munu_2017, &b_munu_2017);
     else if(config.WZSelection=="zmumu")
       fChain->SetBranchAddress("mumu_2017", &mumu_2017, &b_mumu_2017);
   }

   if(config.DataYears=="2018"){
     if(config.WZSelection=="wminus" || config.WZSelection=="wplus")
       fChain->SetBranchAddress("munu_2018", &munu_2018, &b_munu_2018);
     else if(config.WZSelection=="zmumu")
       fChain->SetBranchAddress("mumu_2018", &mumu_2018, &b_mumu_2018);
   }

   fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu50", &mu_trigMatch_HLT_mu50, &b_mu_trigMatch_HLT_mu50);
   if(false){//remove/fix this?
     fChain->SetBranchAddress("HLT_mu20_iloose_L1MU15", &HLT_mu20_iloose_L1MU15, &b_HLT_mu20_iloose_L1MU15);
     fChain->SetBranchAddress("mu_trigMatch_HLT_mu20_iloose_L1MU15", &mu_trigMatch_HLT_mu20_iloose_L1MU15, &b_mu_trigMatch_HLT_mu20_iloose_L1MU15);
   }
   fChain->SetBranchAddress("primaryVertices", &primaryVertices, &b_primaryVertices);
   if(!isData) fChain->SetBranchAddress("weight_pileup_hash", &weight_pileup_hash, &b_weight_pileup_hash);

   Notify();
}

Bool_t MyWZAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyWZAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
/* Int_t MyWZAnalysis::Cut(Long64_t entry) */
/* { */
/* // This function may be called from Loop. */
/* // returns  1 if entry is accepted. */
/* // returns -1 otherwise. */
/*    return 1; */
/* } */
#endif // #ifdef MyWZAnalysis_cxx
