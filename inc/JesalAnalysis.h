//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  6 13:13:06 2019 by ROOT version 6.14/04
// from TTree nominal/tree
// found on file: /data/morales/atlas/Nov2018/mc16d/wminmunu/user.anramire.15915015._000001.mc_sel.root
//////////////////////////////////////////////////////////

#ifndef JesalAnalysis_h
#define JesalAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class JesalAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         weight_mc;
   Float_t         weight_pileup;
   Float_t         weight_leptonSF;
   Float_t         weight_globalLeptonTriggerSF;
   Float_t         weight_oldTriggerSF;
   Float_t         weight_jvt;
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
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP;
   Float_t         weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_Isol_SYST_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_STAT_DOWN;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_UP;
   Float_t         weight_leptonSF_MU_SF_TTVA_SYST_DOWN;
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
   Float_t         weight_indiv_SF_EL_Reco;
   Float_t         weight_indiv_SF_EL_Reco_UP;
   Float_t         weight_indiv_SF_EL_Reco_DOWN;
   Float_t         weight_indiv_SF_EL_ID;
   Float_t         weight_indiv_SF_EL_ID_UP;
   Float_t         weight_indiv_SF_EL_ID_DOWN;
   Float_t         weight_indiv_SF_EL_Isol;
   Float_t         weight_indiv_SF_EL_Isol_UP;
   Float_t         weight_indiv_SF_EL_Isol_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeID;
   Float_t         weight_indiv_SF_EL_ChargeID_UP;
   Float_t         weight_indiv_SF_EL_ChargeID_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeMisID;
   Float_t         weight_indiv_SF_EL_ChargeMisID_STAT_UP;
   Float_t         weight_indiv_SF_EL_ChargeMisID_STAT_DOWN;
   Float_t         weight_indiv_SF_EL_ChargeMisID_SYST_UP;
   Float_t         weight_indiv_SF_EL_ChargeMisID_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_ID;
   Float_t         weight_indiv_SF_MU_ID_STAT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_DOWN;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_STAT_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_UP;
   Float_t         weight_indiv_SF_MU_ID_SYST_LOWPT_DOWN;
   Float_t         weight_indiv_SF_MU_Isol;
   Float_t         weight_indiv_SF_MU_Isol_STAT_UP;
   Float_t         weight_indiv_SF_MU_Isol_STAT_DOWN;
   Float_t         weight_indiv_SF_MU_Isol_SYST_UP;
   Float_t         weight_indiv_SF_MU_Isol_SYST_DOWN;
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
   std::vector<float>   *el_pt;
   std::vector<float>   *el_eta;
   std::vector<float>   *el_cl_eta;
   std::vector<float>   *el_phi;
   std::vector<float>   *el_e;
   std::vector<float>   *el_charge;
   std::vector<float>   *el_topoetcone20;
   std::vector<float>   *el_ptvarcone20;
   std::vector<char>    *el_CF;
   std::vector<float>   *el_d0sig;
   std::vector<float>   *el_delta_z0_sintheta;
   std::vector<char>    *m_el_ECIDS;
   std::vector<double>  *m_el_ECIDSResult;
   std::vector<int>     *el_true_type;
   std::vector<int>     *el_true_origin;
   std::vector<int>     *el_true_firstEgMotherTruthType;
   std::vector<int>     *el_true_firstEgMotherTruthOrigin;
   std::vector<int>     *el_true_firstEgMotherPdgId;
   std::vector<char>    *el_true_isPrompt;
   std::vector<char>    *el_true_isChargeFl;
   std::vector<float>   *mu_pt;
   std::vector<float>   *mu_eta;
   std::vector<float>   *mu_phi;
   std::vector<float>   *mu_e;
   std::vector<float>   *mu_charge;
   std::vector<float>   *mu_topoetcone20;
   std::vector<float>   *mu_ptvarcone30;
   std::vector<float>   *mu_d0sig;
   std::vector<float>   *mu_delta_z0_sintheta;
   std::vector<int>     *mu_true_type;
   std::vector<int>     *mu_true_origin;
   std::vector<char>    *mu_true_isPrompt;
   std::vector<float>   *jet_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_e;
   std::vector<float>   *jet_mv2c00;
   std::vector<float>   *jet_mv2c10;
   std::vector<float>   *jet_mv2c20;
   std::vector<float>   *jet_ip3dsv1;
   std::vector<float>   *jet_jvt;
   std::vector<char>    *jet_passfjvt;
   std::vector<int>     *jet_truthflav;
   std::vector<int>     *jet_truthPartonLabel;
   std::vector<char>    *jet_isTrueHS;
   std::vector<int>     *jet_truthflavExtended;
   std::vector<float>   *jet_MV2c10mu;
   std::vector<float>   *jet_MV2c10rnn;
   std::vector<float>   *jet_DL1;
   std::vector<float>   *jet_DL1mu;
   std::vector<float>   *jet_DL1rnn;
   std::vector<float>   *jet_MV2cl100;
   std::vector<float>   *jet_MV2c100;
   std::vector<float>   *jet_DL1_pu;
   std::vector<float>   *jet_DL1_pc;
   std::vector<float>   *jet_DL1_pb;
   std::vector<float>   *jet_DL1mu_pu;
   std::vector<float>   *jet_DL1mu_pc;
   std::vector<float>   *jet_DL1mu_pb;
   std::vector<float>   *jet_DL1rnn_pu;
   std::vector<float>   *jet_DL1rnn_pc;
   std::vector<float>   *jet_DL1rnn_pb;
   Float_t         met_met;
   Float_t         met_phi;
   Int_t           munu_2017;
   Char_t          HLT_mu50;
   Char_t          HLT_mu26_ivarmedium;
   std::vector<char>    *mu_trigMatch_HLT_mu50;
   std::vector<char>    *mu_trigMatch_HLT_mu26_ivarmedium;
   ULong64_t       weight_pileup_hash;
   Float_t         bornMass_KFactor;
   Double_t        weight_KFactor;
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
   Int_t           primaryVertices;
   Float_t         hadronic_pt;
   Float_t         hadronic_eta;
   Float_t         hadronic_phi;
   Float_t         hadronic_e;
   Double_t        sumET_PFO;
   std::vector<float>   *mu_topoetcone30;
   std::vector<float>   *mu_topoetcone40;
   std::vector<float>   *mu_ptvarcone20;
   std::vector<float>   *mu_ptvarcone40;
   std::vector<float>   *mu_ptcone20;
   std::vector<float>   *mu_ptcone30;
   std::vector<float>   *mu_ptcone40;
   std::vector<int>     *muon_gradientIsolation;//it was twice declared
   std::vector<int>     *muon_trigger_mu20_iloose_L1MU15;
   std::vector<int>     *muon_trigger_mu50;

   // List of branches
   TBranch        *b_weight_mc;   //!
   TBranch        *b_weight_pileup;   //!
   TBranch        *b_weight_leptonSF;   //!
   TBranch        *b_weight_globalLeptonTriggerSF;   //!
   TBranch        *b_weight_oldTriggerSF;   //!
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
   TBranch        *b_weight_indiv_SF_EL_Reco;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Reco_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ID;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ID_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_Isol_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeID_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_STAT_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_STAT_DOWN;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_SYST_UP;   //!
   TBranch        *b_weight_indiv_SF_EL_ChargeMisID_SYST_DOWN;   //!
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
   TBranch        *b_el_CF;   //!
   TBranch        *b_el_d0sig;   //!
   TBranch        *b_el_delta_z0_sintheta;   //!
   TBranch        *b_m_el_ECIDS;   //!
   TBranch        *b_m_el_ECIDSResult;   //!
   TBranch        *b_el_true_type;   //!
   TBranch        *b_el_true_origin;   //!
   TBranch        *b_el_true_firstEgMotherTruthType;   //!
   TBranch        *b_el_true_firstEgMotherTruthOrigin;   //!
   TBranch        *b_el_true_firstEgMotherPdgId;   //!
   TBranch        *b_el_true_isPrompt;   //!
   TBranch        *b_el_true_isChargeFl;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_e;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_topoetcone20;   //!
   TBranch        *b_mu_ptvarcone30;   //!
   TBranch        *b_mu_d0sig;   //!
   TBranch        *b_mu_delta_z0_sintheta;   //!
   TBranch        *b_mu_true_type;   //!
   TBranch        *b_mu_true_origin;   //!
   TBranch        *b_mu_true_isPrompt;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_mv2c00;   //!
   TBranch        *b_jet_mv2c10;   //!
   TBranch        *b_jet_mv2c20;   //!
   TBranch        *b_jet_ip3dsv1;   //!
   TBranch        *b_jet_jvt;   //!
   TBranch        *b_jet_passfjvt;   //!
   TBranch        *b_jet_truthflav;   //!
   TBranch        *b_jet_truthPartonLabel;   //!
   TBranch        *b_jet_isTrueHS;   //!
   TBranch        *b_jet_truthflavExtended;   //!
   TBranch        *b_jet_MV2c10mu;   //!
   TBranch        *b_jet_MV2c10rnn;   //!
   TBranch        *b_jet_DL1;   //!
   TBranch        *b_jet_DL1mu;   //!
   TBranch        *b_jet_DL1rnn;   //!
   TBranch        *b_jet_MV2cl100;   //!
   TBranch        *b_jet_MV2c100;   //!
   TBranch        *b_jet_DL1_pu;   //!
   TBranch        *b_jet_DL1_pc;   //!
   TBranch        *b_jet_DL1_pb;   //!
   TBranch        *b_jet_DL1mu_pu;   //!
   TBranch        *b_jet_DL1mu_pc;   //!
   TBranch        *b_jet_DL1mu_pb;   //!
   TBranch        *b_jet_DL1rnn_pu;   //!
   TBranch        *b_jet_DL1rnn_pc;   //!
   TBranch        *b_jet_DL1rnn_pb;   //!
   TBranch        *b_met_met;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_munu_2017;   //!
   TBranch        *b_HLT_mu50;   //!
   TBranch        *b_HLT_mu26_ivarmedium;   //!
   TBranch        *b_mu_trigMatch_HLT_mu50;   //!
   TBranch        *b_mu_trigMatch_HLT_mu26_ivarmedium;   //!
   TBranch        *b_weight_pileup_hash;   //!
   TBranch        *b_bornMass_KFactor;   //!
   TBranch        *b_weight_KFactor;   //!
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
   TBranch        *b_primaryVertices;   //!
   TBranch        *b_hadronic_pt;   //!
   TBranch        *b_hadronic_eta;   //!
   TBranch        *b_hadronic_phi;   //!
   TBranch        *b_hadronic_e;   //!
   TBranch        *b_sumET_PFO;   //!
   TBranch        *b_mu_topoetcone30;   //!
   TBranch        *b_mu_topoetcone40;   //!
   TBranch        *b_mu_ptvarcone20;   //!
   TBranch        *b_mu_ptvarcone40;   //!
   TBranch        *b_mu_ptcone20;   //!
   TBranch        *b_mu_ptcone30;   //!
   TBranch        *b_mu_ptcone40;   //!
   TBranch        *b_muon_gradientIsolation;   //!
   //   TBranch        *b_muon_gradientIsolation;   //!
   TBranch        *b_muon_trigger_mu20_iloose_L1MU15;   //!
   TBranch        *b_muon_trigger_mu50;   //!

   JesalAnalysis(TTree *tree=0);
   virtual ~JesalAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef JesalAnalysis_cxx
JesalAnalysis::JesalAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/morales/atlas/Nov2018/mc16d/wminmunu/user.anramire.15915015._000001.mc_sel.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/morales/atlas/Nov2018/mc16d/wminmunu/user.anramire.15915015._000001.mc_sel.root");
      }
      f->GetObject("nominal",tree);

   }
   Init(tree);
}

JesalAnalysis::~JesalAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t JesalAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t JesalAnalysis::LoadTree(Long64_t entry)
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

void JesalAnalysis::Init(TTree *tree)
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
   el_CF = 0;
   el_d0sig = 0;
   el_delta_z0_sintheta = 0;
   m_el_ECIDS = 0;
   m_el_ECIDSResult = 0;
   el_true_type = 0;
   el_true_origin = 0;
   el_true_firstEgMotherTruthType = 0;
   el_true_firstEgMotherTruthOrigin = 0;
   el_true_firstEgMotherPdgId = 0;
   el_true_isPrompt = 0;
   el_true_isChargeFl = 0;
   mu_pt = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_e = 0;
   mu_charge = 0;
   mu_topoetcone20 = 0;
   mu_ptvarcone30 = 0;
   mu_d0sig = 0;
   mu_delta_z0_sintheta = 0;
   mu_true_type = 0;
   mu_true_origin = 0;
   mu_true_isPrompt = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_e = 0;
   jet_mv2c00 = 0;
   jet_mv2c10 = 0;
   jet_mv2c20 = 0;
   jet_ip3dsv1 = 0;
   jet_jvt = 0;
   jet_passfjvt = 0;
   jet_truthflav = 0;
   jet_truthPartonLabel = 0;
   jet_isTrueHS = 0;
   jet_truthflavExtended = 0;
   jet_MV2c10mu = 0;
   jet_MV2c10rnn = 0;
   jet_DL1 = 0;
   jet_DL1mu = 0;
   jet_DL1rnn = 0;
   jet_MV2cl100 = 0;
   jet_MV2c100 = 0;
   jet_DL1_pu = 0;
   jet_DL1_pc = 0;
   jet_DL1_pb = 0;
   jet_DL1mu_pu = 0;
   jet_DL1mu_pc = 0;
   jet_DL1mu_pb = 0;
   jet_DL1rnn_pu = 0;
   jet_DL1rnn_pc = 0;
   jet_DL1rnn_pb = 0;
   mu_trigMatch_HLT_mu50 = 0;
   mu_trigMatch_HLT_mu26_ivarmedium = 0;
   mu_topoetcone30 = 0;
   mu_topoetcone40 = 0;
   mu_ptvarcone20 = 0;
   mu_ptvarcone40 = 0;
   mu_ptcone20 = 0;
   mu_ptcone30 = 0;
   mu_ptcone40 = 0;
   muon_gradientIsolation = 0;
   //   muon_gradientIsolation = 0;
   muon_trigger_mu20_iloose_L1MU15 = 0;
   muon_trigger_mu50 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   fChain->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
   fChain->SetBranchAddress("weight_leptonSF", &weight_leptonSF, &b_weight_leptonSF);
   fChain->SetBranchAddress("weight_globalLeptonTriggerSF", &weight_globalLeptonTriggerSF, &b_weight_globalLeptonTriggerSF);
   fChain->SetBranchAddress("weight_oldTriggerSF", &weight_oldTriggerSF, &b_weight_oldTriggerSF);
   fChain->SetBranchAddress("weight_jvt", &weight_jvt, &b_weight_jvt);
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
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco", &weight_indiv_SF_EL_Reco, &b_weight_indiv_SF_EL_Reco);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_UP", &weight_indiv_SF_EL_Reco_UP, &b_weight_indiv_SF_EL_Reco_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Reco_DOWN", &weight_indiv_SF_EL_Reco_DOWN, &b_weight_indiv_SF_EL_Reco_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID", &weight_indiv_SF_EL_ID, &b_weight_indiv_SF_EL_ID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID_UP", &weight_indiv_SF_EL_ID_UP, &b_weight_indiv_SF_EL_ID_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ID_DOWN", &weight_indiv_SF_EL_ID_DOWN, &b_weight_indiv_SF_EL_ID_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol", &weight_indiv_SF_EL_Isol, &b_weight_indiv_SF_EL_Isol);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_UP", &weight_indiv_SF_EL_Isol_UP, &b_weight_indiv_SF_EL_Isol_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_Isol_DOWN", &weight_indiv_SF_EL_Isol_DOWN, &b_weight_indiv_SF_EL_Isol_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID", &weight_indiv_SF_EL_ChargeID, &b_weight_indiv_SF_EL_ChargeID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID_UP", &weight_indiv_SF_EL_ChargeID_UP, &b_weight_indiv_SF_EL_ChargeID_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeID_DOWN", &weight_indiv_SF_EL_ChargeID_DOWN, &b_weight_indiv_SF_EL_ChargeID_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID", &weight_indiv_SF_EL_ChargeMisID, &b_weight_indiv_SF_EL_ChargeMisID);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_STAT_UP", &weight_indiv_SF_EL_ChargeMisID_STAT_UP, &b_weight_indiv_SF_EL_ChargeMisID_STAT_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_STAT_DOWN", &weight_indiv_SF_EL_ChargeMisID_STAT_DOWN, &b_weight_indiv_SF_EL_ChargeMisID_STAT_DOWN);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_SYST_UP", &weight_indiv_SF_EL_ChargeMisID_SYST_UP, &b_weight_indiv_SF_EL_ChargeMisID_SYST_UP);
   fChain->SetBranchAddress("weight_indiv_SF_EL_ChargeMisID_SYST_DOWN", &weight_indiv_SF_EL_ChargeMisID_SYST_DOWN, &b_weight_indiv_SF_EL_ChargeMisID_SYST_DOWN);
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
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("randomRunNumber", &randomRunNumber, &b_randomRunNumber);
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
   fChain->SetBranchAddress("el_CF", &el_CF, &b_el_CF);
   fChain->SetBranchAddress("el_d0sig", &el_d0sig, &b_el_d0sig);
   fChain->SetBranchAddress("el_delta_z0_sintheta", &el_delta_z0_sintheta, &b_el_delta_z0_sintheta);
   fChain->SetBranchAddress("m_el_ECIDS", &m_el_ECIDS, &b_m_el_ECIDS);
   fChain->SetBranchAddress("m_el_ECIDSResult", &m_el_ECIDSResult, &b_m_el_ECIDSResult);
   fChain->SetBranchAddress("el_true_type", &el_true_type, &b_el_true_type);
   fChain->SetBranchAddress("el_true_origin", &el_true_origin, &b_el_true_origin);
   fChain->SetBranchAddress("el_true_firstEgMotherTruthType", &el_true_firstEgMotherTruthType, &b_el_true_firstEgMotherTruthType);
   fChain->SetBranchAddress("el_true_firstEgMotherTruthOrigin", &el_true_firstEgMotherTruthOrigin, &b_el_true_firstEgMotherTruthOrigin);
   fChain->SetBranchAddress("el_true_firstEgMotherPdgId", &el_true_firstEgMotherPdgId, &b_el_true_firstEgMotherPdgId);
   fChain->SetBranchAddress("el_true_isPrompt", &el_true_isPrompt, &b_el_true_isPrompt);
   fChain->SetBranchAddress("el_true_isChargeFl", &el_true_isChargeFl, &b_el_true_isChargeFl);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_e", &mu_e, &b_mu_e);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_topoetcone20", &mu_topoetcone20, &b_mu_topoetcone20);
   fChain->SetBranchAddress("mu_ptvarcone30", &mu_ptvarcone30, &b_mu_ptvarcone30);
   fChain->SetBranchAddress("mu_d0sig", &mu_d0sig, &b_mu_d0sig);
   fChain->SetBranchAddress("mu_delta_z0_sintheta", &mu_delta_z0_sintheta, &b_mu_delta_z0_sintheta);
   fChain->SetBranchAddress("mu_true_type", &mu_true_type, &b_mu_true_type);
   fChain->SetBranchAddress("mu_true_origin", &mu_true_origin, &b_mu_true_origin);
   fChain->SetBranchAddress("mu_true_isPrompt", &mu_true_isPrompt, &b_mu_true_isPrompt);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_mv2c00", &jet_mv2c00, &b_jet_mv2c00);
   fChain->SetBranchAddress("jet_mv2c10", &jet_mv2c10, &b_jet_mv2c10);
   fChain->SetBranchAddress("jet_mv2c20", &jet_mv2c20, &b_jet_mv2c20);
   fChain->SetBranchAddress("jet_ip3dsv1", &jet_ip3dsv1, &b_jet_ip3dsv1);
   fChain->SetBranchAddress("jet_jvt", &jet_jvt, &b_jet_jvt);
   fChain->SetBranchAddress("jet_passfjvt", &jet_passfjvt, &b_jet_passfjvt);
   fChain->SetBranchAddress("jet_truthflav", &jet_truthflav, &b_jet_truthflav);
   fChain->SetBranchAddress("jet_truthPartonLabel", &jet_truthPartonLabel, &b_jet_truthPartonLabel);
   fChain->SetBranchAddress("jet_isTrueHS", &jet_isTrueHS, &b_jet_isTrueHS);
   fChain->SetBranchAddress("jet_truthflavExtended", &jet_truthflavExtended, &b_jet_truthflavExtended);
   fChain->SetBranchAddress("jet_MV2c10mu", &jet_MV2c10mu, &b_jet_MV2c10mu);
   fChain->SetBranchAddress("jet_MV2c10rnn", &jet_MV2c10rnn, &b_jet_MV2c10rnn);
   fChain->SetBranchAddress("jet_DL1", &jet_DL1, &b_jet_DL1);
   fChain->SetBranchAddress("jet_DL1mu", &jet_DL1mu, &b_jet_DL1mu);
   fChain->SetBranchAddress("jet_DL1rnn", &jet_DL1rnn, &b_jet_DL1rnn);
   fChain->SetBranchAddress("jet_MV2cl100", &jet_MV2cl100, &b_jet_MV2cl100);
   fChain->SetBranchAddress("jet_MV2c100", &jet_MV2c100, &b_jet_MV2c100);
   fChain->SetBranchAddress("jet_DL1_pu", &jet_DL1_pu, &b_jet_DL1_pu);
   fChain->SetBranchAddress("jet_DL1_pc", &jet_DL1_pc, &b_jet_DL1_pc);
   fChain->SetBranchAddress("jet_DL1_pb", &jet_DL1_pb, &b_jet_DL1_pb);
   fChain->SetBranchAddress("jet_DL1mu_pu", &jet_DL1mu_pu, &b_jet_DL1mu_pu);
   fChain->SetBranchAddress("jet_DL1mu_pc", &jet_DL1mu_pc, &b_jet_DL1mu_pc);
   fChain->SetBranchAddress("jet_DL1mu_pb", &jet_DL1mu_pb, &b_jet_DL1mu_pb);
   fChain->SetBranchAddress("jet_DL1rnn_pu", &jet_DL1rnn_pu, &b_jet_DL1rnn_pu);
   fChain->SetBranchAddress("jet_DL1rnn_pc", &jet_DL1rnn_pc, &b_jet_DL1rnn_pc);
   fChain->SetBranchAddress("jet_DL1rnn_pb", &jet_DL1rnn_pb, &b_jet_DL1rnn_pb);
   fChain->SetBranchAddress("met_met", &met_met, &b_met_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("munu_2017", &munu_2017, &b_munu_2017);
   fChain->SetBranchAddress("HLT_mu50", &HLT_mu50, &b_HLT_mu50);
   fChain->SetBranchAddress("HLT_mu26_ivarmedium", &HLT_mu26_ivarmedium, &b_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu50", &mu_trigMatch_HLT_mu50, &b_mu_trigMatch_HLT_mu50);
   fChain->SetBranchAddress("mu_trigMatch_HLT_mu26_ivarmedium", &mu_trigMatch_HLT_mu26_ivarmedium, &b_mu_trigMatch_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("weight_pileup_hash", &weight_pileup_hash, &b_weight_pileup_hash);
   fChain->SetBranchAddress("bornMass_KFactor", &bornMass_KFactor, &b_bornMass_KFactor);
   fChain->SetBranchAddress("weight_KFactor", &weight_KFactor, &b_weight_KFactor);
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
   fChain->SetBranchAddress("primaryVertices", &primaryVertices, &b_primaryVertices);
   fChain->SetBranchAddress("hadronic_pt", &hadronic_pt, &b_hadronic_pt);
   fChain->SetBranchAddress("hadronic_eta", &hadronic_eta, &b_hadronic_eta);
   fChain->SetBranchAddress("hadronic_phi", &hadronic_phi, &b_hadronic_phi);
   fChain->SetBranchAddress("hadronic_e", &hadronic_e, &b_hadronic_e);
   fChain->SetBranchAddress("sumET_PFO", &sumET_PFO, &b_sumET_PFO);
   fChain->SetBranchAddress("mu_topoetcone30", &mu_topoetcone30, &b_mu_topoetcone30);
   fChain->SetBranchAddress("mu_topoetcone40", &mu_topoetcone40, &b_mu_topoetcone40);
   fChain->SetBranchAddress("mu_ptvarcone20", &mu_ptvarcone20, &b_mu_ptvarcone20);
   fChain->SetBranchAddress("mu_ptvarcone40", &mu_ptvarcone40, &b_mu_ptvarcone40);
   fChain->SetBranchAddress("mu_ptcone20", &mu_ptcone20, &b_mu_ptcone20);
   fChain->SetBranchAddress("mu_ptcone30", &mu_ptcone30, &b_mu_ptcone30);
   fChain->SetBranchAddress("mu_ptcone40", &mu_ptcone40, &b_mu_ptcone40);
   fChain->SetBranchAddress("muon_gradientIsolation", &muon_gradientIsolation, &b_muon_gradientIsolation);
//    fChain->SetBranchAddress("muon_gradientIsolation", &muon_gradientIsolation, &b_muon_gradientIsolation);
   fChain->SetBranchAddress("muon_trigger_mu20_iloose_L1MU15", &muon_trigger_mu20_iloose_L1MU15, &b_muon_trigger_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("muon_trigger_mu50", &muon_trigger_mu50, &b_muon_trigger_mu50);
   Notify();
}

Bool_t JesalAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JesalAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t JesalAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef JesalAnalysis_cxx
