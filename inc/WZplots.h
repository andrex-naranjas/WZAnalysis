#ifndef WZPLOTS_H
#define WZPLOTS_H

#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h"//kColors
#include "THStack.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TAttText.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TString.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include <string>
#include <vector>
#include <algorithm>//std::find

#include "AtlasStyle.h"
#include "AtlasUtils.h"

#include "ConfigSettings.h"


class WZplots{

 public:
  WZplots();
  virtual ~WZplots();

  UInt_t   textsize, fontNumber;
  Double_t g_legend_width;
  UInt_t   legendTextSize;
  Double_t g_legend_height_per_entry;
  Double_t ratioPadFraction;
  double ratio_leftbottommargin_rightmargin;
  std::string outFileName;
  bool   useAsymmetricErrors;
  int    errorBandFillStyle;
  double errorBandLineWeight;
  double errorBandLineSpacing;
  int    errorBandColor;
  double GeV, xScale, xmin, xmax;
  bool logy, logx, yRangeUserDefined;
  double ymin, ymax;
  //bool ratio = false;
  bool doRatio, fakeData, p0Fit;
  double yminRatio, ymaxRatio;
  unsigned int nDivisionsRatioYaxis; 
  double legend_x1, legend_y1, titleX, titleY, atlasX, atlasY;
  bool labelATLAS, labelStatus, labelStatus2, labelLumi;
  std::string titleLaTeX;
  TString yAxisTitle, xAxisUnitsString;
  std::string lumiString;
  std::string statusString;
  std::string status2String;
  std::string xAxisTitle;
  std::string yAxisTitleRatioPlot;
  std::string errorBandLegendEntry;
  bool errorBandLegendEntryTwoLines;
  std::string errorBandLegendEntryLineTwo;
  std::string binWidthPrecisionFormatString;
  std::vector<TH1F*> legendExcludes; // list of sample histograms to exclude form the legend : addHistToListNoLegend
  std::vector<TH1F*> histList;
  TH1F* h_data;
  TH1F* h_sys;
  TH1F* h_err_up;
  TH1F* h_err_dn;
  //FILE *file;
  Config configPlotsWZ;

  //Functions!
  void SetPlotParameters(Config config, std::string nameOfFile, double scale, double minx, double maxx,
			 bool Ylog, bool Xlog, double ratioYmin, double ratioYmax,
			 bool userDefinedY, double miny, double maxy, double x1_legend,
			 double y1_legend, double Xatlas, double Yatlas, TString latexTitle,
			 TString units, std::string titleXaxis, std::string YlabelPrecision,
			 TH1F *sys, TH1F *data, bool fitP0);

  void printPlot(FILE *fp, std::string names);

  // forward-declarations
  TGraphAsymmErrors* buildAsymmetricErrors(TH1F* stack, TH1F* errup, TH1F* errdn, bool doratio);
  
  TH1F* buildStackSum();
  TH1F* buildData();
  TH1F* buildRatio(TH1F* data, TH1F* sum); 
  TCanvas* buildCanvas();
  TH1F* buildMainCanvas();
  TPad* buildMainPad(TCanvas* c);
  TPad* buildRatioPad(TCanvas* c);
  THStack* buildStack(TPad* mainpad);
  TH1F* buildSystematicBand();
  TH1F* buildSystematicBandRatio(TH1F* stacksum);
  void  configureRatioFrameAxes(TH1F* frame); 
  TString buildBinSizeString(TH1F* h);
  void removeXaxis(THStack* s);
  TLegend* buildLegend(Config config, TH1F* data, TH1F* sys, TGraphAsymmErrors* = NULL);
  void ConfigureLegend(TLegend* legend, std::string header = "", int nheaderlines=1);
  void drawLumiAndStatus();
  //not forward declarations?

  
  virtual TH1F* scaleXaxis(TH1F* hist, double scale);
  virtual void scaleXaxis(std::vector<TH1F*>* stackList, double scale);
  void SetControlPlotStyle(Double_t ratio_frac=0.3, Double_t margin=0.05, Int_t font=43, Double_t size=25);
  void ATLAS_LABEL_unscaled(Double_t x,Double_t y,unsigned int scaling = 3);
  int getProcessColor(std::string proc);
  std::string getProcessLegendEntry(std::string proc);
  virtual void addHistToList(TH1F* hist, std::string processName, bool legend);
  virtual void addHistToList(TH1F* hist, std::string processName);
  virtual void addHistToList(std::vector<TH1F*>* s, TH1F* h, std::string name, double scale);
  void addHistToListNoLegend(TH1F* hist, std::string processName);
  virtual void addSignalToStack(TH1F* hist, std::string processName);
  virtual void addSignalToStack(std::vector<TH1F*>* s, TH1F* h, std::string name, double scale);

};

#endif
