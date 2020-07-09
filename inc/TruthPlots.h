//TruthPlots includes
#ifndef TRUTHPLOTS_H
#define TRUTHPLOTS_H

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


class TruthPlots{

 public:
  TruthPlots();
  virtual ~TruthPlots();

  virtual void  initialize(Config config);
  virtual void  execute(Config config);
  virtual void  setstyle();

 private:

  std::string dirYear, year;
  std::string dirInclusive, total;
  std::string puname;
  std::string systematic;
  std::string wzchannel, nameOfSample;

  TFile *fmc, *fmcp, *fcw, *fcwp;

  virtual void SinglePlot(Config config, std::string kine, std::string boson);
  
  virtual void LevelsPlot(Config config, std::string kine, std::string boson);

  virtual void MatchPlot(Config config, std::string kine, std::string boson);

  virtual void CwPlots(Config config, std::string factor);

  TH1D* scaleXaxis(TH1D* hist, double scale, std::string kine);
  
  virtual void plotAxisLine(TH1D* hist, int lineColor, int markerColor,
			    int markerStyle, double markerSize,
			    TString title, TString xlabel, TString ylabel, bool xRange,
			    double xlow, double xhigh, bool yRange, double ylow, double yhigh);

  virtual void profAxisLine(TProfile* hist, int lineColor, int markerColor,
			    int markerStyle, double markerSize,
			    TString title, TString xlabel, TString ylabel, bool xRange,
			    double xlow, double xhigh, bool yRange, double ylow, double yhigh);

  virtual void ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
			     double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
			     int markerstyle, int color, double markersize, int divisions);

  virtual void profRatioSettings(TProfile* prof, double min, double max, TString xlabel, TString ylabel,
			     double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
			     int markerstyle, int color, double markersize, int divisions);

};

#endif //> !PRINTPLOTS_H
