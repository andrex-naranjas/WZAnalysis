//CwCalculator includes
#ifndef CWCALCULATOR_H
#define CWCALCULATOR_H

#include "ConfigSettings.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <iostream>

#include <TROOT.h>
#include "TLine.h"

#include "TFile.h"
#include "TH1D.h"

class CwCalculator{

 public:
  CwCalculator();
  virtual ~CwCalculator();

  virtual void initialize(Config config);
  virtual void execute(Config config);
  virtual void finalize();

 private:
  void LatexTable(std::string fileName, std::vector<std::string> eta_bin_name, std::string channel);
  std::vector<std::string> LatexBinName();


  TFile *fmcreco, *fmctruth, *fout;
  TH1D *hReco, *hTruth, *hMwtCstayTruth, *hMwtCleaveTruth, *hMwtCstayReco, *hMwtCcomeReco, *hCwz;
  TH1D *hRecoEta, *hTruthEta, *hEtaCstayTruth, *hEtaCleaveTruth, *hEtaCstayReco, *hEtaCcomeReco, *hCwzEta;
  TH1D *hPtywz, *hStabwz, *hPtywzEta, *hStabwzEta;
  TH1D *hCw_simple, *hCw_eta_simple;

  std::vector<double> Nrec_temp;
  std::vector<double> Ngen_temp;
  std::vector<double> NstayRec_temp;
  std::vector<double> NstayGen_temp;
  std::vector<double> NleaveGen_temp;
  std::vector<double> NComeRec_temp;
  std::vector<double> dCwOld_temp;
  std::vector<double> dCwNew_temp;
  std::vector<double> Cw_temp;

};

#endif //> !CWCALCULATOR
