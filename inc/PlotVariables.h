#ifndef PLOTVARIABLES_H
#define PLOTVARIABLES_H

#include <iostream>
#include "TROOT.h"
#include "TStyle.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TColor.h" // kColors
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
#include <algorithm> //std::find


UInt_t   textsize;
UInt_t   fontNumber;
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
double GeV;
double xScale;
double xmin;
double xmax;
bool logy;
bool logx;
bool   yRangeUserDefined;
double ymin;
double ymax;
//bool ratio = false;
bool doRatio;
bool fakeData;
double yminRatio;
double ymaxRatio;
unsigned int nDivisionsRatioYaxis;

double legend_x1;
double legend_y1;
double titleX;
double titleY;
double atlasX;
double atlasY;

bool labelATLAS = true;
bool labelStatus = true;
bool labelStatus2 = false;
bool labelLumi = true;

std::string titleLaTeX;
TString yAxisTitle;
TString xAxisUnitsString;
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

#endif
