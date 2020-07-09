//TruthPlots includes
#ifndef SFPLOTS_H
#define SFPLOTS_H

#include "ConfigSettings.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <THStack.h>

#include "TCanvas.h"
#include "TPad.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveText.h"

#include <string>
#include <vector>


class SFPlots{

 public:
  SFPlots();
  virtual ~SFPlots();

  virtual void  initialize(Config config);
  virtual void  execute(Config config);
  virtual void  setstyle();

 private:
  std::string dirYear, year;
  std::string dirInclusive, total;
  std::string puname;
  std::string systematic;
  std::string wzchannel, nameOfSample;

  TFile *fnom, *fup, *fdown;

  virtual void SinglePlot(Config config, std::string kine, std::string boson, std::string SF);

  virtual void CombinedPlot(Config config, std::string kine, std::string boson);
  
  TH1D* scaleXaxis(TProfile* hist, double scale, std::string kine);
  

  virtual void plotAxisLine(TH1D* hist, int lineColor, int markerColor,
			    int markerStyle, double markerSize,
			    TString title, TString xlabel, TString ylabel, bool xRange,
			    double xlow, double xhigh, bool yRange, double ylow, double yhigh);
};

#endif //> !PRINTPLOTS_H
