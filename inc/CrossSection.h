//CrossSection includes
#ifndef CROSSSECTION_H
#define CROSSSECTION_H

#include "ConfigSettings.h"

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include <string>
#include <vector>
#include <iostream>

#include <TROOT.h>
#include "TLine.h"

#include "TCanvas.h"
#include "TPad.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveText.h"

#include "TFile.h"
#include "TH1D.h"


class CrossSection{

 public:
  CrossSection();
  virtual ~CrossSection();

  virtual void initialize(Config config);
  virtual void execute(Config config);
  virtual void finalize();

 private:
  void InclusiveXsec(Config config);
  void EtaBinsXsection(Config config);
  double getDataLumi(Config config);
  double asymError(double xSec, double xSecP, double error, double errorP);

  void sysTable(std::vector<std::vector<double>> totalValues,
		std::vector<std::string> systematic, std::vector<std::string> header, std::string fileName);

  void xSecTableIncl(Config config, std::vector<std::vector<double>> totalValues, std::vector<std::string> column, std::vector<std::string> header, std::string fileName);

  void asyTable(Config config, std::vector<std::vector<double>> totalValues, 
		std::vector<std::string> column, std::vector<std::string> header,
		std::string fileName);

  void xSecTable(Config config, std::vector<std::vector<double>> totalValues,
		 std::vector<std::string> columnLabel, std::vector<std::string> header, std::string fileName);


  void xSecTableSum(std::vector<std::vector<double>> totalValues,std::vector<std::vector<double>> stat,
		    std::vector<std::vector<double>> sys_up, std::vector<std::vector<double>> sys_down, 
		    std::vector<std::string> column, std::vector<std::string> header, std::string fileName);
  
  void xSecPlot(Config config, TH1D *xSecHist, TH1D *xSecHistMC, std::string boson, double ylow, double yhigh);

  void plotAxisLine(TH1D* hist, int lineColor, int markerColor, int markerStyle, double markerSize,
		    TString title, TString xlabel, TString ylabel, bool xRange,
		    double xlow, double xhigh, bool yRange, double ylow, double yhigh);


  void ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
		     double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
		     int markerstyle, int color, double markersize, int divisions);


  std::vector<std::string> LatexBinName();
  std::vector<std::string> Names(std::string option);

  void setstyle();

  TFile *fdata, *fmc, *fbg, *fcw, *fsys, *fdatap, *fmcp, *fbgp, *fcwp, *fsysp, *fout;
  TFile *fmultijet, *fmultijetp;
  TH1D *hEtaBinXsec, *hEtaBinXsecP, *hEtaBinAsymm;
  TH1D *hEtaBinXsecMC, *hEtaBinXsecMCP, *hEtaBinAsymmMC;
  //debug
  double scale, scalep;


};

#endif
