#ifndef RECOILCALIBRATION_H
#define RECOILCALIBRATION_H

#include <TH2.h>
#include <TVector2.h>
#include <TF1.h>
#include <TProfile2D.h>
#include <TFile.h>

#include "TLorentzVector.h"

#include <string>
#include <vector>


class RecoilCalibration { 
 public: 
  RecoilCalibration();
  virtual ~RecoilCalibration();

  virtual void  initialize(TString calibFile, TString channel="z");
  virtual void  execute(TLorentzVector trueboson, double sumet, TLorentzVector& u, float& weight);
  virtual void  applySumetReweighting(double ptv, double sumet, float& weight);
  void applyUX(double sumet, TLorentzVector& u);
  void applyUY(double sumet, TLorentzVector& u);
  void applyUXbin(double sumet, TLorentzVector& u);
  void applyUYbin(double sumet, TLorentzVector& u);
  void applyResponseAndResolCorr(TLorentzVector trueboson, double sumet, TLorentzVector& u);
  virtual void  finalize();
  
 private: 
  TFile *calibfile;
  TProfile *hux, *huy;
  TH2F *h; //histo for SET reweighting
  TH2D *h2_resol; //histo for resolution correction
  TF1 *fux, *fuy;
  TH2D* h2_resp;
  TProfile2D *p2_resp, *p2_upar_data; //histo for response correction
  float m_pi;
  //histograms
  // Kinematic functions
  float deltaR(const float eta1,const float eta2,const float phi1,const float phi2) const;
  float deltaPhi(const float phi1,const float phi2) const;
  float deltaEta(const float eta1,const float eta2) const;


}; 

#endif //> !RECOILCALIBRATION_H
