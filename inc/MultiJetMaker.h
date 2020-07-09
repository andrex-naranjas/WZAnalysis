//MultiJetMaker includes
#ifndef MULTIJETMAKER_H
#define MULTIJETMAKER_H

#include "ConfigSettings.h"

#include <string>
#include <vector>
#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TMath.h>

#include "TFractionFitter.h"
#include "TObjArray.h"


class MultiJetMaker{

 public:
  MultiJetMaker();
  virtual ~MultiJetMaker();

  virtual void initialize(Config config);
  virtual void execute(Config config);
  virtual void finalize();

 private:

  TH1D* MultiJetShape(Config config, std::string kine, std::string eta_bin,std::string iso);

  void SingleNormalisation(TH1D *data, TH1D *mc0, TH1D *mc1,
			   double &fracMu,double &fracMC,
			   double &errMu, double &errMC,
			   std::string kine, std::string eta_bin);

  void MultiJetNorm(Config config, std::string kine, std::string eta_bin,std::string label);

  void MultijetFinal(Config config, int it, TH1D* h, std::string kine, std::string eta_bin);

  Int_t status;

  double scale=1.;

  std::vector<double> FinalMJFracMET; 
  std::vector<double> FinalMJFracMWT;

  std::vector<double> FinalMJFracErrorMET; 
  std::vector<double> FinalMJFracErrorMWT;

  std::vector<double> OriginalMJFracMET; 
  std::vector<double> OriginalMJFracMWT;

  std::vector<double> finalFractionAv;
  std::vector<double> finalFractionErrorAv;

  std::vector<int> status_met; 
  std::vector<int> status_mwt;

  std::string wzchannel;
  std::string dirYear;
  std::string fileDir;

  std::vector<double> multijet_eta_bin;
  std::vector<double> mc_eta_bin;
  std::vector<double> data_eta_bin;

  TH1D *hData_eta, *hMC_eta, *hMC_eta_multi;//final eta histos

  TFile *fdata0, *fmc0, *fdata1, *fmc1, *fdata2, *fmc2, *fout;
  TFile *fdata0_p, *fmc0_p;
  
  TH1D *hData, *hBG, *hData1, *hBG1, *hData2, *hBG2;
  TH1D *hData_nod0, *hBG_nod0;
  TH1D *hDataAntid0, *hBGAntid0, *hData1Antid0, *hBG1Antid0, *hData2Antid0, *hBG2Antid0;
  TH1D *hMultiSR, *hMultiFR1, *hMultiFR2;

  TH1D *hBG1_norm, *hMultiFR1_norm, *hData1_norm;
  TH1D *hBG2_norm, *hMultiFR2_norm, *hData2_norm;

  TH1D *hBG1_norm_scaled, *hMultiFR1_norm_scaled;
  TH1D *hBG2_norm_scaled, *hMultiFR2_norm_scaled;

  TH1D *h_sum_fr1, *h_sum_fr2;
  TH1D *h_sum_fr1_norm_scaled;
  TH1D *h_sum_fr1_scaled, *h_sum_fr2_scaled;

  TH1D *hBG1_scaled, *hMultiFR1_scaled;

  TH1D *result, *result_norm;

  double MultiFrac1=1., MCFrac1=1.;
  double MultiFrac2=1., MCFrac2=1.; 

  double MultiErr1=1., MCErr1=1.;
  double MultiErr2=1., MCErr2=1.; 
  
};

#endif
