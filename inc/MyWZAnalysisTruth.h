
#ifndef MyWZAnalysisTruth_h
#define MyWZAnalysisTruth_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
#include "vector"
#include <vector>

#include "ConfigSettings.h"

class MyWZAnalysisTruth {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   std::string     nameOfSample;
   Config config;

   TTree          *fChain_sumWeights;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain_reco;   //!pointer to the analyzed TTree or TChain

   Double_t        lumi_sample;
   Double_t        cut_up;
   Double_t        cut_down;
   Double_t        cut_eta;
   Double_t        cut_pt1;
   Double_t        cut_pt2;
   Double_t        cut_met1;
   Double_t        cut_met2;
   Double_t        cut_mwt1;
   Double_t        cut_mwt2;
   
   // Declaration of leaf types
   Float_t         weight_mc;
   ULong64_t       eventNumber;
   UInt_t          runNumber;
   Float_t         mu;
   Float_t         weight_pileup;
   Double_t        KFactor_weight_truth;
   Double_t        KFactor_bornMass_truth;
   UInt_t          randomRunNumber;
   UInt_t          mcChannelNumber;
   std::vector<float>   *PDFinfo_X1;
   std::vector<float>   *PDFinfo_X2;
   std::vector<int>     *PDFinfo_PDGID1;
   std::vector<int>     *PDFinfo_PDGID2;
   std::vector<float>   *PDFinfo_Q;
   std::vector<float>   *PDFinfo_XF1;
   std::vector<float>   *PDFinfo_XF2;
   Float_t         MC_WZneutrino_eta_born;
   Float_t         MC_WZneutrino_pt_born;
   Float_t         MC_WZneutrino_m_born;
   Float_t         MC_WZmu_el_phi_born;
   Float_t         MC_WZmu_el_pt_born;
   Float_t         MC_WZ_phi;
   Float_t         MC_WZ_pt;
   Float_t         MC_WZ_m;
   Float_t         MC_WZ_dilep_m_dres;
   Float_t         MC_WZ_dilep_eta_dres;
   Float_t         MC_WZneutrino_eta_dres;
   Float_t         MC_WZ_eta;
   Float_t         MC_WZmu_el_phi_bare;
   Float_t         MC_WZ_dilep_eta_bare;
   Float_t         MC_WZmu_el_m_bare;
   Float_t         MC_WZneutrino_phi_dres;
   Float_t         MC_WZ_dilep_eta_born;
   Float_t         MC_WZmu_el_m_dres;
   Float_t         MC_WZneutrino_pt_dres;
   Float_t         MC_WZneutrino_eta_bare;
   Float_t         MC_WZ_dilep_phi_born;
   Float_t         MC_WZ_dilep_pt_bare;
   Float_t         MC_WZmu_el_pt_bare;
   Float_t         MC_WZneutrino_phi_bare;
   Float_t         MC_WZneutrino_m_bare;
   Float_t         MC_WZ_dilep_pt_born;
   Float_t         MC_WZneutrino_phi_born;
   Float_t         MC_WZmu_el_m_born;
   Float_t         MC_WZneutrino_m_dres;
   Float_t         MC_WZneutrino_pt_bare;
   Float_t         MC_WZmu_el_eta_bare;
   Float_t         MC_WZmu_el_pt_dres;
   Float_t         MC_WZ_dilep_m_bare;
   Float_t         MC_WZ_dilep_pt_dres;
   Float_t         MC_WZ_dilep_m_born;
   Float_t         MC_WZmu_el_phi_dres;
   Float_t         MC_WZ_dilep_phi_bare;
   Float_t         MC_WZmu_el_eta_born;
   Float_t         MC_WZ_dilep_phi_dres;
   Float_t         MC_WZmu_el_eta_dres;

   // List of branches
   TBranch        *b_weight_mc;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_mu;   //!
   TBranch        *b_weight_pileup;   //!
   TBranch        *b_KFactor_weight_truth;   //!
   TBranch        *b_KFactor_bornMass_truth;   //!
   TBranch        *b_randomRunNumber;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_PDFinfo_X1;   //!
   TBranch        *b_PDFinfo_X2;   //!
   TBranch        *b_PDFinfo_PDGID1;   //!
   TBranch        *b_PDFinfo_PDGID2;   //!
   TBranch        *b_PDFinfo_Q;   //!
   TBranch        *b_PDFinfo_XF1;   //!
   TBranch        *b_PDFinfo_XF2;   //!
   TBranch        *b_MC_WZneutrino_eta_born;   //!
   TBranch        *b_MC_WZneutrino_pt_born;   //!
   TBranch        *b_MC_WZneutrino_m_born;   //!
   TBranch        *b_MC_WZmu_el_phi_born;   //!
   TBranch        *b_MC_WZmu_el_pt_born;   //!
   TBranch        *b_MC_WZ_phi;   //!
   TBranch        *b_MC_WZ_pt;   //!
   TBranch        *b_MC_WZ_m;   //!
   TBranch        *b_MC_WZ_dilep_m_dres;   //!
   TBranch        *b_MC_WZ_dilep_eta_dres;   //!
   TBranch        *b_MC_WZneutrino_eta_dres;   //!
   TBranch        *b_MC_WZ_eta;   //!
   TBranch        *b_MC_WZmu_el_phi_bare;   //!
   TBranch        *b_MC_WZ_dilep_eta_bare;   //!
   TBranch        *b_MC_WZmu_el_m_bare;   //!
   TBranch        *b_MC_WZneutrino_phi_dres;   //!
   TBranch        *b_MC_WZ_dilep_eta_born;   //!
   TBranch        *b_MC_WZmu_el_m_dres;   //!
   TBranch        *b_MC_WZneutrino_pt_dres;   //!
   TBranch        *b_MC_WZneutrino_eta_bare;   //!
   TBranch        *b_MC_WZ_dilep_phi_born;   //!
   TBranch        *b_MC_WZ_dilep_pt_bare;   //!
   TBranch        *b_MC_WZmu_el_pt_bare;   //!
   TBranch        *b_MC_WZneutrino_phi_bare;   //!
   TBranch        *b_MC_WZneutrino_m_bare;   //!
   TBranch        *b_MC_WZ_dilep_pt_born;   //!
   TBranch        *b_MC_WZneutrino_phi_born;   //!
   TBranch        *b_MC_WZmu_el_m_born;   //!
   TBranch        *b_MC_WZneutrino_m_dres;   //!
   TBranch        *b_MC_WZneutrino_pt_bare;   //!
   TBranch        *b_MC_WZmu_el_eta_bare;   //!
   TBranch        *b_MC_WZmu_el_pt_dres;   //!
   TBranch        *b_MC_WZ_dilep_m_bare;   //!
   TBranch        *b_MC_WZ_dilep_pt_dres;   //!
   TBranch        *b_MC_WZ_dilep_m_born;   //!
   TBranch        *b_MC_WZmu_el_phi_dres;   //!
   TBranch        *b_MC_WZ_dilep_phi_bare;   //!
   TBranch        *b_MC_WZmu_el_eta_born;   //!
   TBranch        *b_MC_WZ_dilep_phi_dres;   //!
   TBranch        *b_MC_WZmu_el_eta_dres;   //!

   MyWZAnalysisTruth(Config, std::string,
		double, double, double,
		double, double, double, 
		double, double, double, double);

   virtual ~MyWZAnalysisTruth();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyWZAnalysisTruth_cxx
MyWZAnalysisTruth::MyWZAnalysisTruth(Config myconfig, std::string samples,
			   double luminosity_samples, double up_cut, double down_cut,
			   double eta_cut, double pt_cut1, double pt_cut2, 
				     double met_cut1, double met_cut2, double mwt_cut1, double mwt_cut2) : fChain(0),fChain_sumWeights(0),fChain_reco(0)
{
  //andres
  config = myconfig;
  nameOfSample = samples;
  lumi_sample = luminosity_samples;
  cut_up = up_cut;
  cut_down = down_cut;
  cut_eta = eta_cut;
  cut_pt1 = pt_cut1;
  cut_pt2 = pt_cut2;
  cut_met1 = met_cut1;
  cut_met2 = met_cut2;
  cut_mwt1 = mwt_cut1;
  cut_mwt2 = mwt_cut2;

  TTree *tree;

  //Reco analysis
  TChain *chain = new TChain("truth","");
  chain->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
  tree=chain;
  Init(tree);

  //sum of weights
  TChain *chain_sumWeights = new TChain("sumWeights","");
  chain_sumWeights->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
  fChain_sumWeights=chain_sumWeights;
  
  //reco events, maaaps
  if(config.RecoMatching=="True"){
    TChain *chain_reco = new TChain("nominal","");
    chain_reco->Add((config.InputFileDir+nameOfSample+"/*.root").c_str());
    fChain_reco=chain_reco;
  }
  
}

MyWZAnalysisTruth::~MyWZAnalysisTruth()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyWZAnalysisTruth::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyWZAnalysisTruth::LoadTree(Long64_t entry)
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

void MyWZAnalysisTruth::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PDFinfo_X1 = 0;
   PDFinfo_X2 = 0;
   PDFinfo_PDGID1 = 0;
   PDFinfo_PDGID2 = 0;
   PDFinfo_Q = 0;
   PDFinfo_XF1 = 0;
   PDFinfo_XF2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("mu", &mu, &b_mu);
   fChain->SetBranchAddress("weight_pileup", &weight_pileup, &b_weight_pileup);
   fChain->SetBranchAddress("KFactor_weight_truth", &KFactor_weight_truth, &b_KFactor_weight_truth);
   fChain->SetBranchAddress("KFactor_bornMass_truth", &KFactor_bornMass_truth, &b_KFactor_bornMass_truth);
   fChain->SetBranchAddress("randomRunNumber", &randomRunNumber, &b_randomRunNumber);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("PDFinfo_X1", &PDFinfo_X1, &b_PDFinfo_X1);
   fChain->SetBranchAddress("PDFinfo_X2", &PDFinfo_X2, &b_PDFinfo_X2);
   fChain->SetBranchAddress("PDFinfo_PDGID1", &PDFinfo_PDGID1, &b_PDFinfo_PDGID1);
   fChain->SetBranchAddress("PDFinfo_PDGID2", &PDFinfo_PDGID2, &b_PDFinfo_PDGID2);
   fChain->SetBranchAddress("PDFinfo_Q", &PDFinfo_Q, &b_PDFinfo_Q);
   fChain->SetBranchAddress("PDFinfo_XF1", &PDFinfo_XF1, &b_PDFinfo_XF1);
   fChain->SetBranchAddress("PDFinfo_XF2", &PDFinfo_XF2, &b_PDFinfo_XF2);
   fChain->SetBranchAddress("MC_WZneutrino_eta_born", &MC_WZneutrino_eta_born, &b_MC_WZneutrino_eta_born);
   fChain->SetBranchAddress("MC_WZneutrino_pt_born", &MC_WZneutrino_pt_born, &b_MC_WZneutrino_pt_born);
   fChain->SetBranchAddress("MC_WZneutrino_m_born", &MC_WZneutrino_m_born, &b_MC_WZneutrino_m_born);
   fChain->SetBranchAddress("MC_WZmu_el_phi_born", &MC_WZmu_el_phi_born, &b_MC_WZmu_el_phi_born);
   fChain->SetBranchAddress("MC_WZmu_el_pt_born", &MC_WZmu_el_pt_born, &b_MC_WZmu_el_pt_born);
   fChain->SetBranchAddress("MC_WZ_phi", &MC_WZ_phi, &b_MC_WZ_phi);
   fChain->SetBranchAddress("MC_WZ_pt", &MC_WZ_pt, &b_MC_WZ_pt);
   fChain->SetBranchAddress("MC_WZ_m", &MC_WZ_m, &b_MC_WZ_m);
   fChain->SetBranchAddress("MC_WZ_dilep_m_dres", &MC_WZ_dilep_m_dres, &b_MC_WZ_dilep_m_dres);
   fChain->SetBranchAddress("MC_WZ_dilep_eta_dres", &MC_WZ_dilep_eta_dres, &b_MC_WZ_dilep_eta_dres);
   fChain->SetBranchAddress("MC_WZneutrino_eta_dres", &MC_WZneutrino_eta_dres, &b_MC_WZneutrino_eta_dres);
   fChain->SetBranchAddress("MC_WZ_eta", &MC_WZ_eta, &b_MC_WZ_eta);
   fChain->SetBranchAddress("MC_WZmu_el_phi_bare", &MC_WZmu_el_phi_bare, &b_MC_WZmu_el_phi_bare);
   fChain->SetBranchAddress("MC_WZ_dilep_eta_bare", &MC_WZ_dilep_eta_bare, &b_MC_WZ_dilep_eta_bare);
   fChain->SetBranchAddress("MC_WZmu_el_m_bare", &MC_WZmu_el_m_bare, &b_MC_WZmu_el_m_bare);
   fChain->SetBranchAddress("MC_WZneutrino_phi_dres", &MC_WZneutrino_phi_dres, &b_MC_WZneutrino_phi_dres);
   fChain->SetBranchAddress("MC_WZ_dilep_eta_born", &MC_WZ_dilep_eta_born, &b_MC_WZ_dilep_eta_born);
   fChain->SetBranchAddress("MC_WZmu_el_m_dres", &MC_WZmu_el_m_dres, &b_MC_WZmu_el_m_dres);
   fChain->SetBranchAddress("MC_WZneutrino_pt_dres", &MC_WZneutrino_pt_dres, &b_MC_WZneutrino_pt_dres);
   fChain->SetBranchAddress("MC_WZneutrino_eta_bare", &MC_WZneutrino_eta_bare, &b_MC_WZneutrino_eta_bare);
   fChain->SetBranchAddress("MC_WZ_dilep_phi_born", &MC_WZ_dilep_phi_born, &b_MC_WZ_dilep_phi_born);
   fChain->SetBranchAddress("MC_WZ_dilep_pt_bare", &MC_WZ_dilep_pt_bare, &b_MC_WZ_dilep_pt_bare);
   fChain->SetBranchAddress("MC_WZmu_el_pt_bare", &MC_WZmu_el_pt_bare, &b_MC_WZmu_el_pt_bare);
   fChain->SetBranchAddress("MC_WZneutrino_phi_bare", &MC_WZneutrino_phi_bare, &b_MC_WZneutrino_phi_bare);
   fChain->SetBranchAddress("MC_WZneutrino_m_bare", &MC_WZneutrino_m_bare, &b_MC_WZneutrino_m_bare);
   fChain->SetBranchAddress("MC_WZ_dilep_pt_born", &MC_WZ_dilep_pt_born, &b_MC_WZ_dilep_pt_born);
   fChain->SetBranchAddress("MC_WZneutrino_phi_born", &MC_WZneutrino_phi_born, &b_MC_WZneutrino_phi_born);
   fChain->SetBranchAddress("MC_WZmu_el_m_born", &MC_WZmu_el_m_born, &b_MC_WZmu_el_m_born);
   fChain->SetBranchAddress("MC_WZneutrino_m_dres", &MC_WZneutrino_m_dres, &b_MC_WZneutrino_m_dres);
   fChain->SetBranchAddress("MC_WZneutrino_pt_bare", &MC_WZneutrino_pt_bare, &b_MC_WZneutrino_pt_bare);
   fChain->SetBranchAddress("MC_WZmu_el_eta_bare", &MC_WZmu_el_eta_bare, &b_MC_WZmu_el_eta_bare);
   fChain->SetBranchAddress("MC_WZmu_el_pt_dres", &MC_WZmu_el_pt_dres, &b_MC_WZmu_el_pt_dres);
   fChain->SetBranchAddress("MC_WZ_dilep_m_bare", &MC_WZ_dilep_m_bare, &b_MC_WZ_dilep_m_bare);
   fChain->SetBranchAddress("MC_WZ_dilep_pt_dres", &MC_WZ_dilep_pt_dres, &b_MC_WZ_dilep_pt_dres);
   fChain->SetBranchAddress("MC_WZ_dilep_m_born", &MC_WZ_dilep_m_born, &b_MC_WZ_dilep_m_born);
   fChain->SetBranchAddress("MC_WZmu_el_phi_dres", &MC_WZmu_el_phi_dres, &b_MC_WZmu_el_phi_dres);
   fChain->SetBranchAddress("MC_WZ_dilep_phi_bare", &MC_WZ_dilep_phi_bare, &b_MC_WZ_dilep_phi_bare);
   fChain->SetBranchAddress("MC_WZmu_el_eta_born", &MC_WZmu_el_eta_born, &b_MC_WZmu_el_eta_born);
   fChain->SetBranchAddress("MC_WZ_dilep_phi_dres", &MC_WZ_dilep_phi_dres, &b_MC_WZ_dilep_phi_dres);
   fChain->SetBranchAddress("MC_WZmu_el_eta_dres", &MC_WZmu_el_eta_dres, &b_MC_WZmu_el_eta_dres);
   Notify();
}

Bool_t MyWZAnalysisTruth::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyWZAnalysisTruth::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyWZAnalysisTruth::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyWZAnalysisTruth_cxx
