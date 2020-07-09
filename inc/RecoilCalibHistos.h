#ifndef RECOILCALIBHISTOS_H
#define RECOILCALIBHISTOS_H

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TMath.h>
#include "TLorentzVector.h"

#include <string>
#include <vector>


class RecoilCalibHistos { 
 public: 
  RecoilCalibHistos();
  virtual ~RecoilCalibHistos();

  virtual void  initialize(TString dir=".",TString name="");
   virtual  void  execute(TLorentzVector tboson, //std::vector<int> truthCharge,
				std::vector<TLorentzVector> recoLeptons,  //std::vector<int> recoCharge,
				TLorentzVector upfo, double sumet, float mu, double recoWeight);
  virtual void  scalehists(double scale=1.);

  virtual void  finalize();
  
 private: 

  //histograms
  TFile* m_newfile;
  TH1D* h_mu, *h_Sumet, *h_ptvtrue, *h_ptv, *h_upar, *h_uperp, *h_upar_rec, *h_uperp_rec, *h_ut, *h_ux, *h_uy, *h_utphi, *h_met, *h_metphi, *h_upar_ptv, *h_upar_ptv_n;
  TH2D* h_uperp_pfo_vs_set, *h_upar_pfo_vs_set, *h_set_vs_ptvtrue, *h_set_vs_ptv, *h_bias_vs_ptvtrue, *h_bias_vs_ptv, *h_set_vs_ut, *h_bias_vs_set, *h_uperp_pfo_vs_pt;
  TProfile* p_ux_phi, *p_uy_phi, *p_ux_set, *p_uy_set, *p_bias_set, *p_bias_pt, *p_uperp_set, *p_uperp_pt;
  TProfile2D* p_uperp_vs_set_vs_ptvtrue, *p_upar_vs_set_vs_ptvtrue, *p_uperp_vs_set_vs_ptv, *p_upar_vs_set_vs_ptv;
  //TH3 if you do not know how to handle TProfile2D like me...
  TH3D* respcorr_vs_set_vs_ptv, *resolcorr_vs_set_vs_ptv;
 
  // Kinematic functions
  float deltaR(const float eta1,const float eta2,const float phi1,const float phi2) const;
  float deltaPhi(const float phi1,const float phi2) const;
  float deltaEta(const float eta1,const float eta2) const;


};

#endif //> !RECOILCALIBHISTOS_H
