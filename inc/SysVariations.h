//Systematic Uncertainties
#ifndef SYSVARIATIONS_H
#define SYSVARIATIONS_H

#include "ConfigSettings.h"

#include "TString.h"
#include "TFile.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPad.h"
#include "TPaveStats.h"

#include "TCanvas.h"
#include "TPad.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveText.h"


#include "TPad.h"
#include "TH1F.h"
#include "TH1.h"
#include <vector>

class SysVariations{

 public:
  SysVariations();
  virtual ~SysVariations();
  virtual void initialize(Config config);
  virtual void execute(Config config);

 private:
  void SysPlots(Config config,TH2D* h_up_eta, TH2D* h_dn_eta);
  void setstyle();
  bool isScalekFactor(std::string sys);
  bool isScaleEffFactor(std::string sys);
  TFile *newfile;
  TH1D *h_down, *h_downTEMP, *h_stat, *h_statTEMP, *h_nomTEMP, *h_up,*h_upTEMP;
  TH1D *h_downInd, *h_upInd;
  TH2D *h_downInd_eta, *h_upInd_eta;
  std::string wzchannel;
  std::string calib;
  std::string puname;
  std::string dirYear;
  std::vector<std::string> samples;
  std::vector<std::string> systematics_down;
  std::vector<std::string> systematics_up;
  std::vector<std::string> kine;
  std::string dirInclusive, tot;

  //virtual void getErrorUpDown(TH1D* dnT, TH1D* upT, TH1D* stat, TH1D* &finDn, TH1D* &finUp);

};

#endif
