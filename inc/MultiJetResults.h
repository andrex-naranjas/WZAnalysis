//MultiJetMaker includes
#ifndef MULTIJETRESULTS_H
#define MULTIJETRESULTS_H

#include "ConfigSettings.h"

#include <string>
#include <vector>
#include <iostream>
#include <TFile.h>
#include <TH1.h>

#include "TObjArray.h"

#include "TCanvas.h"
#include "TPad.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveText.h"
#include <TMath.h>
#include <TCanvas.h>
#include "TGaxis.h"

#include <string>
#include <vector>


class MultiJetResults{

 public:
  MultiJetResults();
  virtual ~MultiJetResults();

  virtual void initialize(Config config);
  virtual void execute(Config config);
  virtual void finalize();

 private:

  virtual void Plots(Config config);

  void ComparePlot(Config config, std::string kine, std::string option, std::string ratiolabel,
		   TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4,
		   std::string boson, std::string eta_bin, std::string label,
		   double ylow, double yhigh, bool ylog, std::string channel);

  void SinglePlot(Config config, std::string kine, std::string option,
		  TH1D *hist1, TH1D *hist2, TH1D *hist3,TH1D *hist4,
		  std::string cuts, std::string eta_bin, std::string label,
		  double ylow, double yhigh, bool ylog, std::string channel);

  void setstyle();

  void plotAxisLine(TH1D* hist, int lineColor, int markerColor,
		    int markerStyle, double markerSize,
		    TString title, TString xlabel, TString ylabel, bool xRange,
		    double xlow, double xhigh, bool yRange, double ylow, double yhigh);
  
  void ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
		     double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
		     int markerstyle, int color, double markersize, int divisions);

  void MultiTable(std::string fileName, std::vector<std::string> eta_bin_name, std::string channel);
  void TableCompare(std::string fileName, std::vector<std::string> eta_bin_name);

  void TableSys(std::string fileName, std::vector<std::string> eta_bin_name,std::string channel, std::string systematic);
  
  std::vector<std::string> LatexBinName();

  Int_t status;

  std::string wzchannel;
  std::string dirYear;
  std::string fileDir;
  std::string year;
  std::string dirInclusive, total;
  std::string systematic;

  std::vector<double> multijet_eta_bin_wminus;
  std::vector<double> mc_eta_bin_wminus;
  std::vector<double> data_eta_bin_wminus;
  std::vector<double> finalFractions_wminus;
  std::vector<double> finalErrorFractions_wminus;

  std::vector<double> multijet_eta_bin_wplus;
  std::vector<double> mc_eta_bin_wplus;
  std::vector<double> data_eta_bin_wplus;
  std::vector<double> finalFractions_wplus;
  std::vector<double> finalErrorFractions_wplus;

  std::vector<double> multijet_eta_bin_wminus_error;
  std::vector<double> mc_eta_bin_wplus_error;

  std::vector<double> multijet_eta_bin_wminus_d0down;
  std::vector<double> mc_eta_bin_wminus_d0down;

  std::vector<double> multijet_eta_bin_wplus_d0down;
  std::vector<double> mc_eta_bin_wplus_d0down;

  std::vector<double> multijet_eta_bin_wminus_d0up;
  std::vector<double> mc_eta_bin_wminus_d0up;

  std::vector<double> multijet_eta_bin_wplus_d0up;
  std::vector<double> mc_eta_bin_wplus_d0up;

  std::vector<double> multijet_eta_bin_wminus_xSecdown;
  std::vector<double> mc_eta_bin_wminus_xSecdown;

  std::vector<double> multijet_eta_bin_wplus_xSecdown;
  std::vector<double> mc_eta_bin_wplus_xSecdown;

  std::vector<double> multijet_eta_bin_wminus_xSecup;
  std::vector<double> mc_eta_bin_wminus_xSecup;

  std::vector<double> multijet_eta_bin_wplus_xSecup;
  std::vector<double> mc_eta_bin_wplus_xSecup;

  
  std::vector<std::vector<long int>> totalFlow;
  std::vector<std::string> names;

  TFile *fdata0, *fmc0, *fdata1, *fmc1, *fdata2, *fmc2, *fout;
  TFile *fIN;

  //signal
  TH1D *hData,*hBG,*hDataAntid0,*hBGAntid0,*hMultiSR;
  //FR1
  TH1D *hData1,*hBG1,*hData1Antid0,*hBG1Antid0,*hMultiFR1;
  TH1D *hData_nod0,*hBG_nod0;

  //fit results
  TH1D *hData1_norm,*hBG1_norm,*hMultiFR1_norm,*hBG1_norm_scaled;
  TH1D *hMultiFR1_norm_scaled,*h_sum_fr1_norm_scaled;

  TH1D *hBG1_scaled,*hMultiFR1_scaled,*h_sum_fr1_scaled,*h_sum_fr1;
  TH1D *result,*result_norm;

  //shape
  TH1D *hShapeP,*hShapeP_norm;

  //final results
  TH1D *hMultiSR_Final,*hData_Final,* hMC_Final,*hMultijetMC_Final;
  TH1D *hMC_eta,*hMC_eta_multi,*hData_eta;

  TH1D *data_numbers,*mc_numbers,*multijet_numbers,*fraction_numbers,*fractionerror_numbers;
  
};

#endif
