#include "WZplots.h"

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
#include "TPave.h"
#include "TLegendEntry.h"
#include "TAttText.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include <string>
#include <vector>
#include <algorithm> //std::find


WZplots::WZplots()
{
}

WZplots::~WZplots(){}

void WZplots::SetPlotParameters(Config config, std::string nameOfFile, double scale, double minx, double maxx,
			   bool Ylog, bool Xlog, double ratioYmin, double ratioYmax,
			   bool userDefinedY, double miny, double maxy, double x1_legend,
			   double y1_legend, double Xatlas, double Yatlas, TString latexTitle,
			   TString units, std::string titleXaxis, std::string YlabelPrecision,
			   TH1F *sys, TH1F *data, bool fitP0){

  textsize = 25; // in pixel
  fontNumber = 43; // fixed-pixel-size Helvetica with no italics or bold 
  g_legend_width = 0.2;
  legendTextSize = 20;
  g_legend_height_per_entry = 0.06; // fraction of the pad?
  ratioPadFraction = 0.3;
  ratio_leftbottommargin_rightmargin = 2.3;
  outFileName = "DEFAULTPLOTOUTFILENAME.pdf";
  
  useAsymmetricErrors = false;
  errorBandFillStyle   = 3004;//3002; // 3345 black hatching w/ transparent background - displays well in all common image formats.
  errorBandLineWeight  = 0.75;
  errorBandLineSpacing = 1;
  errorBandColor       = kBlack;
  GeV = 1000;
  xScale = 1;  // scaling for x-axis, mostly for scaling MeV histograms to GeV for display
  xmin = 0;    // x-axis range minimum (after rescaling)
  xmax = 1000;  // x-axis range maximum (after rescaling)
  logy = false;
  logx = false;//andres
  yRangeUserDefined    = false;
  ymin                 = 0.001; // not zero in case of log plot
  ymax                 = 1;
  //bool ratio = false;
  doRatio = true;
  fakeData = false; // 
  p0Fit=false;
  yminRatio = 0.69;
  ymaxRatio = 1.31;
  nDivisionsRatioYaxis = 6;//andres

  legend_x1 = 0.60;//0.60 default
  legend_y1 = 0.45;//.45 default
  titleX = 0.15; // unused, preserved for back-compatibility
  titleY = 0.80; // unused, preserved for back-compatibility
  atlasX = 0.3;
  atlasY = 0.78;
  
  labelATLAS = true;
  labelStatus = true;
  labelStatus2 = false;
  labelLumi = true;

  titleLaTeX = "DEFAULTTITLE";
  yAxisTitle = "Entries";
  xAxisUnitsString = "XUNITS";
  lumiString = "37.7 fb^{-1}";       // This will be parsed as LaTeX.
  //statusString = "Analysis";       // "Internal", and "Preliminary" are the official choices
  statusString = "Internal";       // "Internal", and "Preliminary" are the official choices
  status2String = "SIMSTATUS";    // "Simulation," if the data histogram shows MC. 
  xAxisTitle = "XAXIS";     // this is parsed as TLaTeX (ROOT's busted LaTeX imitations) _{xx} is subscript, ^{xx} is superscript
  yAxisTitleRatioPlot = "Data/Pred."; // 
  errorBandLegendEntry = "MC Stat. #oplus Syst. Unc."; 
  errorBandLegendEntryTwoLines = false;
  errorBandLegendEntryLineTwo = "#oplus Syst. Unc.";
  binWidthPrecisionFormatString = "%.0f"; // printf format, " %.{N}f", e.g., "%.0f", where N is the number of decimal points - usually zero, 1 or 2 for eta plots

  // legendExcludes.clear(); // list of sample histograms to exclude form the legend : addHistToListNoLegend
  // histList.clear();
 
  h_data = NULL; //data histogram
  h_sys = NULL;  //systematics histogram
  h_err_up = NULL; // upper asymm. errors
  h_err_dn = NULL; // lower asymm. errors

  //change of variables, custom for each plot
  configPlotsWZ=config;
  if(configPlotsWZ.DataYears=="2017") lumiString = "44.3 fb^{-1}";
  else if(configPlotsWZ.DataYears=="2018") lumiString = "59.9 fb^{-1}";
  else if(configPlotsWZ.DataYears=="2015+2016") lumiString = "36.1 fb^{-1}";
  else if(configPlotsWZ.DataYears=="Full") lumiString = "79.6 fb^{-1}";
  outFileName = nameOfFile;
  xScale=scale;
  xmin=minx;
  xmax=maxx;
  logy=Ylog;
  logx=Xlog;
  p0Fit=fitP0;
  yminRatio=ratioYmin;
  ymaxRatio=ratioYmax;
  yRangeUserDefined=userDefinedY;
  ymin=miny;
  ymax=maxy;
  legend_x1=x1_legend;
  legend_y1=y1_legend;
  atlasX=Xatlas;
  atlasY=Yatlas;
  titleLaTeX=latexTitle;
  xAxisUnitsString=units;
  xAxisTitle=titleXaxis;
  binWidthPrecisionFormatString=YlabelPrecision;
  h_data = data;
  h_sys = sys;

}


TGraphAsymmErrors* WZplots::buildAsymmetricErrors(TH1F* st, TH1F* eu, TH1F* ed, bool doratio)
{
  if(!st) return NULL;
  if(!eu) return NULL;
  if(!ed) return NULL;
  std::cout<<"     ********************************* Andres *****************************      "<<std::endl;
  //TH1F* stack = scaleXaxis(st,xScale); // stack comes pre-scaled
  TH1F* stack = st;
  TH1F* errup = scaleXaxis(eu,xScale);
  TH1F* errdn = scaleXaxis(ed,xScale);
  
  // for now we trust that the three input histograms have identical binning. 
  /*unsigned*/ int nbins = stack->GetNbinsX();
  std::vector<double> xval;     xval.clear();   xval.resize(nbins);
  std::vector<double> yval;     yval.clear();   yval.resize(nbins);
  std::vector<double> xerrdn; xerrdn.clear(); xerrdn.resize(nbins);
  std::vector<double> xerrup; xerrup.clear(); xerrup.resize(nbins);
  std::vector<double> yerrdn; yerrdn.clear(); yerrdn.resize(nbins);
  std::vector<double> yerrup; yerrup.clear(); yerrup.resize(nbins);
  for(int b=0; b < nbins; b++)
    {
      int r = b+1; // root bins start from 1, not 0
      double scale = doratio ? 1.0/stack->GetBinContent(r)  : 1;
      double binwidth = stack->GetBinWidth(r); // scaled stack	
      xval.at(b) = stack->GetBinLowEdge(r) + 0.5*binwidth;
      yval.at(b) = doratio ? 1 : stack->GetBinContent(r);
      xerrdn.at(b) = 0.5*binwidth;
      xerrup.at(b) = 0.5*binwidth;
      yerrdn.at(b) = scale * fabs( errdn->GetBinContent(r) );
      yerrup.at(b) = scale * fabs( errup->GetBinContent(r) ); 
    }
  TGraphAsymmErrors* err = new TGraphAsymmErrors(nbins,&(xval[0]),&(yval[0]),&(xerrdn[0]),&(xerrup[0]),&(yerrdn[0]),&(yerrup[0]));
  err->SetFillStyle(errorBandFillStyle);
  return err;
}

void WZplots::printPlot(FILE* fp, std::string names)
{
  // format histsograms ---------------------------------------------------------
  // if(fakeData && !h_data) h_data = makeFakeData(&histList);
  gStyle->SetErrorX(1);
  std::cout << "building histograms and stacks" << std::endl;
  TH1F* h_dataScaled = buildData();
  //h_dataScaled->SetMarkerSize(1.);//andres
  TH1F* h_sum = buildStackSum();   // sum up the stack for the data/MC ratio - not scaled in x yet.
  TH1F* h_ratio = NULL;

  double chindf = h_dataScaled->Chi2Test(h_sum,"CHI2/NDF");
  double chi = h_dataScaled->Chi2Test(h_sum,"CHI2");
  double ndf = chi/chindf;
  fprintf(fp,"*********************************************\n");	    
  fprintf(fp,(names+"\n").c_str());
  fprintf(fp,"Chi2: %g    Chi2/NDF: %g    NDF: %g \n",chi,chindf,ndf);
  std::cout<<chi<<"    "<<chindf<<"   "<<ndf<<std::endl;

  if(doRatio) h_ratio = buildRatio(h_dataScaled,h_sum); // data/MC ratio

  std::cout << "building canvas and pads" << std::endl;
  TCanvas* canvas = buildCanvas(); 
  TPad* mainPad = buildMainPad(canvas);
  TPad* ratioPad = NULL;
  if(doRatio) ratioPad = buildRatioPad(canvas); 
  
  std::cout << "building stack & error bands" << std::endl;
  THStack* stack = buildStack(mainPad);
  TH1F* sysScaled = NULL;
  TH1F* sysRatio = NULL;
  TGraphAsymmErrors* asymmErr = NULL; 
  TGraphAsymmErrors* asymmErrRatio = NULL;
  if(         !useAsymmetricErrors) sysScaled     = buildSystematicBand();
  if(doRatio && !useAsymmetricErrors) sysRatio      = buildSystematicBandRatio(h_sum);
  if(          useAsymmetricErrors) asymmErr      = buildAsymmetricErrors(h_sum,h_err_up,h_err_dn,false);
  if(doRatio &&  useAsymmetricErrors) asymmErrRatio = buildAsymmetricErrors(h_sum,h_err_up,h_err_dn,true);
  // main pad: draw stack, then error band above, then data on top
  std::cout << "building drawing main pad" << std::endl;
  mainPad->cd();
  stack->Draw();
  if(!useAsymmetricErrors && sysScaled) {sysScaled->Draw("E2 same");}
  if( useAsymmetricErrors && asymmErr ) {asymmErr->Draw("2 same");}
  if(            h_dataScaled         ) {h_dataScaled->Draw("pe X0 same");}

  // ratio pad : check for pad, draw error band, then draw data on top
  TH1F* ratioFrame = NULL;
  if(doRatio && ratioPad)
    {
      ratioPad->cd();
      ratioFrame = ratioPad->DrawFrame(xmin,yminRatio,xmax,ymaxRatio); // DrawFrame returns a TH1F*, so you can set axis properties on that. Insane.
      configureRatioFrameAxes(ratioFrame); 
      ratioFrame->Draw();
      ratioPad->Update();
      canvas->Update();
      // red line on ratio = 1
      TLine* ratioOneLine = new TLine(/*ratioPad->GetUxmin()*/xmin,1,/*ratioPad->GetUxmax()*/xmax,1);//andres
      ratioOneLine->SetLineColor(kRed + 1);
      ratioOneLine->SetLineStyle(1);
      ratioOneLine->SetLineWidth(3);
      ratioOneLine->Draw();
      // draw data on top of systematic band on top of ratio-one line
      if(     sysRatio && !useAsymmetricErrors) {sysRatio->Draw("E2 SAME");}
      if(asymmErrRatio &&  useAsymmetricErrors) {asymmErrRatio->Draw("2 same");}
      
      if(h_ratio && !p0Fit) {h_ratio->Draw("pe X0 SAME");}
      if(h_ratio && p0Fit) {
	h_ratio->Draw("pe X0 SAME");
	TF1 *p0 = new TF1("p0","pol0"); p0->SetLineColor(kBlue-3); 
	h_ratio->Fit("p0");

	if (fp!=NULL) {
	  for (int i=0;i<p0->GetNpar();i++) {
	    Float_t value = p0->GetParameter(i);
	    fprintf(fp,"p0 fit:  %f \n",value);
	  }
	}


      }
    }

  // set titles and ranges ---------------------------------------------------
  yAxisTitle += buildBinSizeString(h_dataScaled);
  // I shoudl probably redo these titles with a TFrame
  // clean this up.
  stack->GetYaxis()->SetTitle(yAxisTitle); 
  stack->GetYaxis()->SetTitleOffset(1.4);
  stack->GetXaxis()->SetTitle(xAxisTitle.c_str()); 
  stack->GetXaxis()->SetTitleOffset(1.25);
  if(h_dataScaled) 
    {
	h_dataScaled->GetYaxis()->SetTitle(yAxisTitle);
	h_dataScaled->GetYaxis()->SetTitleOffset();
	h_dataScaled->GetXaxis()->SetTitle(xAxisTitle.c_str());
	h_dataScaled->GetXaxis()->SetTitleOffset(1.25);
    }
  if(sysScaled ) 
    {
      sysScaled->GetYaxis()->SetTitle(yAxisTitle);
      sysScaled->GetYaxis()->SetTitleOffset();
      sysScaled->GetXaxis()->SetTitle(xAxisTitle.c_str());
      sysScaled->GetXaxis()->SetTitleOffset(1.25);
    }
  
  if(!yRangeUserDefined) 
    {
      double maxHeight = stack->GetMaximum();
      if(sysScaled && !useAsymmetricErrors)
	{
	  for(/*unsigned*/ int b=1; b<(1+sysScaled->GetXaxis()->GetNbins()); b++) // loop over error histogram bins (bin 0 is underflow)
	    {
	      double errHeight = sysScaled->GetBinContent(b) + sysScaled->GetBinError(b); // this presumes symmetric errors
	      if(errHeight > maxHeight) maxHeight = errHeight;
	    }
	}
      if(asymmErr && useAsymmetricErrors)
	{
	  int nbins = asymmErr->GetN();
	  Double_t* val   = asymmErr->GetY();
	  Double_t* errup = asymmErr->GetEYhigh(); // get arrays of graph values & upper errors
	  for(/*unsigned*/ int b=0; b<(nbins); b++) // loop over graph points 
	    {
	      double errHeight = val[b] + errup[b]; 
			if(errHeight > maxHeight) maxHeight = errHeight;
	    }
	}
      stack->SetMaximum(logy ? maxHeight*10 : maxHeight*(1.05)); // add some vertical clearance to y range 
      stack->SetMinimum(ymin);
    }
  else 
    {
      stack->SetMaximum(ymax);
      stack->SetMinimum(ymin);
    }
  if(doRatio) removeXaxis(stack);

// build legend ------------------------------------------------------------
  TLegend* legend = buildLegend(configPlotsWZ, h_dataScaled,sysScaled,asymmErr);
  ConfigureLegend(legend);
  
  mainPad->cd();
  mainPad->Draw();
  legend->Draw();
  // myText(titleX,titleY,kBlack,titleLaTeX.c_str()); // this is added in drawLumiAndStatus, now.
  drawLumiAndStatus();
  
  if(doRatio) 
    {
      ratioPad->cd();
      ratioPad->Draw();
    }

  // print canvas ------------------------------------------------------------
  canvas->Update();
  canvas->Print(outFileName.c_str()/*,outFilePrintOpt.c_str()*/); // uncomment if printing to PDF? 

} // closing printPlot()


void WZplots::configureRatioFrameAxes(TH1F* frame)
{
  if(!frame)return;
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetXaxis()->SetTitleOffset(3.4);
  frame->GetXaxis()->SetLimits(xmin,xmax);
  frame->GetXaxis()->SetNdivisions(6,0,0,true);
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle(yAxisTitleRatioPlot.c_str());
  frame->GetYaxis()->SetLimits(yminRatio,ymaxRatio);
  frame->GetYaxis()->SetNdivisions(nDivisionsRatioYaxis,0,0,true); // less ticks on the ratio axis
  frame->GetYaxis()->SetLabelSize(17.5);//andres
  return;
}

TString WZplots::buildBinSizeString(TH1F* data)
{
  if(!data) return "";
  TString bss = "";
  double bin_width = data->GetXaxis()->GetBinWidth(1);
  std::string format = " / " + binWidthPrecisionFormatString;
  bss = TString::Format(format.c_str(),bin_width);
  bss += " " + xAxisUnitsString;
  return bss;
}
	/*
// this convolutes the errors of the data and MC, by default.
TH1F* buildRatio(TH1F* data, TH1F* sum)
{
if(!data) return NULL;
if(!sum) return NULL;
TH1F* hratio = (TH1F*)data->Clone("ratio");
hratio->Divide(sum);
return hratio;
}
	*/ 

TH1F* WZplots::buildRatio(TH1F* data, TH1F* sum)
{
  if(!data) return NULL;
  if(!sum) return NULL;
  TH1F* hratio = (TH1F*)data->Clone("ratio");
  for(/*unsigned*/ int b=1; b< (2+(data->GetNbinsX())); b++) // loop over bins from underflow to overflow
    {
      double dat = data->GetBinContent(b);
      if(dat == 0) continue;
      double err = data->GetBinError(b);
      double pred = sum->GetBinContent(b);
      if(pred == 0) continue;
      double ratioval = dat/pred;
      double ratioerr = err/pred;
      //std::cout<<ratioval<<"    "<<pred<<"     "<<data<<std::endl;
      // hratio->SetBinContent(b,1); // testing
      hratio->SetBinContent(b,ratioval);
      hratio->SetBinError(b,ratioerr);
	}
  return hratio;
}

TH1F* WZplots::buildStackSum()
{
  // if(!histList) return NULL; // we've always got the histList, now.
  if(histList.size() < 1) 
    {
      std::cout << "(!!) WARNING - NO HISTOGRAMS ADDED TO STACK (!!)" << std::endl;
      return NULL;
    }
  TH1F* sum = (TH1F*)( (histList.at(0))->Clone("sum") );
  for(unsigned int h=1;h<histList.size();h++)
    {
      TH1F* hist = (TH1F*)(histList.at(h));
      sum->Add(hist);
    }
  sum->SetMarkerSize(0); // remove the markers fo when we clone & plot on the ratio pad
  return sum;
}

TH1F* WZplots::buildSystematicBandRatio(TH1F* sum)
{
  if(!h_sys) return NULL;
  if(!sum) return NULL;
  TH1F* ratiosys = (TH1F*)h_sys->Clone("ratioSys");
  ratiosys->Divide(sum);
  TH1F* ratiosysScaled = scaleXaxis(ratiosys,xScale);
  ratiosysScaled->SetFillColor(errorBandColor);
  ratiosysScaled->SetFillStyle(errorBandFillStyle);
  ratiosysScaled->SetMarkerSize(0);
  return ratiosysScaled;
}

TH1F* WZplots::buildSystematicBand()
{
  if (!h_sys) return NULL;
  TH1F* scaled = scaleXaxis(h_sys,xScale);
  scaled->SetFillColor(errorBandColor);
  scaled->SetFillStyle(errorBandFillStyle);
  scaled->SetMarkerSize(0);
  return scaled;
}

TH1F* WZplots::buildData()
{
  if(!h_data) return NULL;
  h_data->SetFillColorAlpha(kWhite,0); // transparent white - only neeed for the legend 
  //h_data->SetMarkerSize(0.5);//andres
  TH1F* d = scaleXaxis(h_data,xScale);
  return d;
}

TPad* WZplots::buildRatioPad(TCanvas* canvas)
{
  if(!doRatio) return NULL;
  double padFraction = 0.3;
  double margin = 0.05;
  double left_bottom_margin = ratio_leftbottommargin_rightmargin*margin;
  TPad* pad = new TPad("ratioPad","ratioPad",0,0,1,padFraction);
  pad->SetTopMargin(0); // clearance between main and ratio pad is handled in the top pad
  pad->SetRightMargin(margin);
  pad->SetLeftMargin(left_bottom_margin);
  pad->SetBottomMargin(left_bottom_margin/padFraction);
  pad->SetGrid(0,2);
  pad->SetLogx(logx);//andres
  canvas->cd();
  pad->Draw();
  return pad;
}

TPad* WZplots::buildMainPad(TCanvas* canvas)
{
  double ratioPadFraction = 0; // default for no ratio plot
  double mainPadBottomMargin = 0.1; // default for no ratio plot
  if(doRatio) 
	{
	ratioPadFraction = 0.3;
	mainPadBottomMargin = 0.02;
	}
  double mainPadFraction = 1-ratioPadFraction;
  double rightmargin = 0.05;
  double bottommargin = 0.02;
  TPad* pad = new TPad("mainPad","mainPad",0, 1-mainPadFraction, 1,1);
  const Double_t left_bottom_margin = ratio_leftbottommargin_rightmargin*rightmargin;
  if(!doRatio) bottommargin = left_bottom_margin;
  pad->SetBottomMargin(bottommargin);
  pad->SetRightMargin(rightmargin);
  pad->SetLeftMargin(left_bottom_margin);
  pad->SetTopMargin(rightmargin/mainPadFraction);
  canvas->cd();
  pad->SetLogy(logy);
  pad->SetLogx(logx);//andres
  pad->Draw();
  return pad;
}

TCanvas* WZplots::buildCanvas()
{
  unsigned int canvasHeightPixels = 600; // default for all plots
  unsigned int canvasWidthPixels = 800; // default for non-ratio plots
  if(doRatio) {canvasWidthPixels = 600;} // default for ratio plots
  TCanvas* canvas = new TCanvas("CGcanvas","CGcanvas",canvasWidthPixels,canvasHeightPixels);
  return canvas;
}

TH1F* WZplots::scaleXaxis(TH1F* hist, double scale)
{
  if(scale == 1) return hist;
  double binMin = hist->GetXaxis()->GetXmin();
  double binMax = hist->GetXaxis()->GetXmax();
  double nBins  = hist->GetXaxis()->GetNbins();
  double binMinScaled = binMin*scale;
  double binMaxScaled = binMax*scale;
  TH1F* scaled = NULL;
//andres 
  if(logx){ 
    printf("  calculating logarithmic binning...\n");
    binMinScaled = log10(binMinScaled); binMaxScaled = log10(binMaxScaled);
    double scaledBins[(int)nBins+1];
    double binWidth = (binMaxScaled - binMinScaled)/(nBins);
    for (int i = 0; i < nBins+1; i++){
      scaledBins[i  ] = binMinScaled+(i  )*binWidth;
      scaledBins[i  ] = pow(10,scaledBins[i  ]);
    }
    scaled = new TH1F(hist->GetName(),hist->GetTitle(),nBins,scaledBins);
    std::cout<<"Scaled histogram!!!!!!!!!"<<std::endl;
  }else {
    scaled = new TH1F(hist->GetName(),hist->GetTitle(),nBins,binMinScaled,binMaxScaled);
  }
  
  //scaled->Sumw2();
  for(unsigned int b = 1; b< nBins+1; b++) // look at bin numbering scheme documented in TH1::GetBin()
    {
      double val = hist->GetBinContent(b);
      double err = hist->GetBinError(b);
      scaled->SetBinContent(b,val);
      scaled->SetBinError(b,err);
    }
  return scaled;
}

void WZplots::scaleXaxis(std::vector<TH1F*>* stackList, double scale)
{
  for(unsigned int h=0;h<stackList->size();h++) // loop over stack histograms (histList?)
    {
      stackList->at(h) = scaleXaxis(stackList->at(h), scale); // reset pointers 
    }
  return;
} // closing scaleXaxis()

void WZplots::ConfigureLegend(TLegend* legend, std::string header, int nheaderlines) 
{


  // const Double_t legend_x2 =legend_x1 + (g_legend_width+(legend->GetNColumns()-1)*legend->GetColumnSeparation())*legend->GetNColumns();
  // const Double_t legend_y2 = legend_y1 + g_legend_height_per_entry*(legend->GetNRows()+(header != "" ? nheaderlines : 0));


  // std::cout<<"         ******************************************************************** PPPPUUUUUUUUUUUUUTO*********************************************************************"<<std::endl;
  // std::cout<<legend_x1<<"      "<<legend_y1<<"     "<<legend_x2<<"     "<<legend_y2<<std::endl;

  // legend->SetX1NDC(legend_x1);
  // legend->SetX2NDC(legend_x2);
  // legend->SetY1NDC(legend_y1);
  // legend->SetY2NDC(legend_y2);

  legend->SetTextAlign(12);
  
  legend->SetTextFont(fontNumber);
  legend->SetTextSize(legendTextSize);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  if(header != "")  {legend->SetHeader(header.c_str());}
} // closing ConfigureLegend


void WZplots::removeXaxis(THStack* stack)
{
  if(!stack) return;
  stack->GetXaxis()->SetLabelSize(0);
  return;
}

void WZplots::SetControlPlotStyle(Double_t ratio_frac, Double_t margin, Int_t font, Double_t size)
{
std::cout << "\nApplying ControlPlot style settings...\n" << std::endl ;
 
 TGaxis* tga = new TGaxis(); 
 tga->SetMaxDigits(4); // this sets a global variable only accessible throufh this class
 
 // AtlasStyle.C default font 42 (Helvetica regular) scales atutomatically, so 
 // it shrinks in the bottom pad. 43 is Helvetica in fixed-pixel-width. 
 
 gStyle->SetTextFont(font);
 
 gStyle->SetTextSize(size);
 gStyle->SetLabelFont(font,"x");
 gStyle->SetTitleFont(font,"x");
 gStyle->SetLabelFont(font,"y");
 gStyle->SetTitleFont(font,"y");
 gStyle->SetLabelFont(font,"z");
 gStyle->SetTitleFont(font,"z");
 
 gStyle->SetLabelSize(size,"x");
 gStyle->SetTitleSize(size,"x");
 gStyle->SetLabelSize(size,"y");
 gStyle->SetTitleSize(size,"y");
 gStyle->SetLabelSize(size,"z");
 gStyle->SetTitleSize(size,"z");
 
 gStyle->SetHatchesSpacing(errorBandLineSpacing);
 gStyle->SetHatchesLineWidth(errorBandLineWeight);
 gStyle->SetLineScalePS(2);
 return;	
}

void WZplots::ATLAS_LABEL_unscaled(Double_t x,Double_t y,unsigned int scaling) 
{
  TLatex l; 
  l.SetNDC();
  unsigned int font = 40 + scaling; // 70 is Helvetica, bold, italics. +3 is fixed-pixel-width (no scaling to pads)
  l.SetTextFont(font);
  l.SetTextColor(kBlack);
  std::string labeltext = "#bf{#it{ATLAS}}";
  //std::string labeltext = "CCDY";
  if(labelStatus) labeltext += (" " +  statusString); // e.g., "Internal"
  l.DrawLatex(x,y,labeltext.c_str());
}

int WZplots::getProcessColor(std::string proc)
{
  int color = kGreen;
  if      (proc == "w_enu"        ){color = kCyan-10; }
  else if (proc == "w_taunu"      ){color = kMagenta+1; }
  else if (proc == "wplus_taunu"  ){color = kMagenta+1;;}
  else if (proc == "wminus_taunu" ){color = kMagenta+1;}
  else if (proc == "w_munu"       ){color = kGreen+2; }
  else if (proc == "wplus_munu"   ){color = kGreen+2; }
  else if (proc == "wminus_munu"  ){color = kGreen+2; }
  else if (proc == "z_ee"         ){color = kBlue-3; }
  else if (proc == "z_mumu"       ){color = kGreen-7; }
  else if (proc == "z_tautau"     ){color = kYellow-3; }
  else if (proc == "top"          ){color = kOrange+1; }
  else if (proc == "multijet"     ){color = kRed+1; }
  else if (proc == "diboson"      ){color = kAzure+7; }
  // less used, but no reason to remove them.
  else if (proc == "wplus_enu"    ){color = kBlue+1;}
  else if (proc == "wminus_enu"   ){color = kRed+1;}
  // else if (proc == "wplus_taunu"  ){color = kYellow;}
  // else if (proc == "wminus_taunu" ){color = kMagenta+1;}
  else if (proc == "ttbar"        ){color = kOrange+1;}
  else if (proc == "JF17"         ){color = kRed-1; }
  else if (proc == "JF23"         ){color = kRed-2; }
  else if (proc == "JF35"         ){color = kRed-3; }
  else if (proc == "JF50"         ){color = kRed-4; }
  else                             {color = kGreen; std::cout << "(!!) UNRECOGNIZED PROCESS NAME (!!)" << std::endl;}
  return color;
}

//andres

std::string WZplots::getProcessLegendEntry(std::string proc)
{
  std::string legendEntry = "what?";
  if      (proc == "w_enu"       ){legendEntry = "W^{+}#rightarrowe^{+}#nu";             }
  else if (proc == "w_taunu"     ){legendEntry = "W#rightarrow#tau#nu";          }
  else if (proc == "w_munu"      ){legendEntry = "W#rightarrow#mu#nu";           }
  else if (proc == "z_ee"        ){legendEntry = "Z#rightarrowe^{+}e^{-}";       }
  else if (proc == "z_mumu"      ){legendEntry = "Z#rightarrow#mu^{+}#mu^{-}";   }
  else if (proc == "z_tautau"    ){legendEntry = "Z#rightarrow#tau^{+}#tau^{-}"; }
  else if (proc == "top"         ){legendEntry = "Top";           }
  else if (proc == "multijet"    ){legendEntry = "Multijet";                      }
  else if (proc == "diboson"     ){legendEntry = "Diboson";                       }
  // less used, but no reason to remove them.
  else if (proc == "ttbar"       ){legendEntry = "t#bar{t}";   }
  else if (proc == "wplus_enu"   ){legendEntry = "W^{+}#rightarrowe^{+}#nu_{e}";  }
  else if (proc == "wminus_enu"  ){legendEntry = "W^{-}#rightarrowe^{-}#nu_{e}";  }
  else if (proc == "wplus_munu"   ){legendEntry = "W^{+}#rightarrow#mu^{+}#nu_{#mu}";  }
  else if (proc == "wminus_munu"  ){legendEntry = "W^{-}#rightarrow#mu^{-}#nu_{#mu}";  }
  else if (proc == "wplus_taunu" ){legendEntry = "W^{+}#rightarrow#tau^{+}#nu_{#tau}"; }
  else if (proc == "wminus_taunu"){legendEntry = "W^{-}#rightarrow#tau^{-}#nu_{#tau}";}
  else if (proc == "JF17"        ){legendEntry = "JF17";}
  else if (proc == "JF23"        ){legendEntry = "JF23";}
  else if (proc == "JF35"        ){legendEntry = "JF35";}
  else if (proc == "JF50"        ){legendEntry = "JF50";}
  else                            {legendEntry = "( ?? )"; std::cout << "(!!) UNRECOGNIZED PROCESS NAME (!!)" <<proc<< std::endl;}
  return legendEntry;
}

void WZplots::drawLumiAndStatus()
{
  double titleY = atlasY + 0.06;
  myText(atlasX,titleY,kBlack,titleLaTeX.c_str());
  // const TString prefix="#lower[-0.2]{#scale[0.6]{#int}}Ldt = "; // old style
  const TString prefix = "13 TeV, ";
  const TString lumi = lumiString; 
  const TString lumiLaTeX = prefix+lumi;
  if(labelATLAS) ATLAS_LABEL_unscaled(atlasX,atlasY);
  double status2Y = atlasY;
  if(labelStatus2) status2Y -= 0.05;
  if(labelStatus2) myText(atlasX,status2Y,kBlack,status2String.c_str());
  double lumiY = status2Y - 0.08;
  if(labelLumi) myText(atlasX,lumiY,kBlack,lumiLaTeX);
}

void WZplots::addHistToList(TH1F* hist, std::string processName, bool legend)
{
  hist = scaleXaxis(hist,xScale);
  int color = getProcessColor(processName);
  std::string title = getProcessLegendEntry(processName);
  hist->SetTitle(title.c_str());
  hist->SetFillColor(color);
  histList.push_back(hist);
  if(!legend) legendExcludes.push_back(hist);
  return;
}

// trivial overload to invoke with default legend choice
void WZplots::addHistToList(TH1F* hist, std::string processName)
{
  addHistToList(hist, processName, true); // do add it to the legend, by default
  return;
}

void WZplots::addHistToListNoLegend(TH1F* hist, std::string processName)
{
  addHistToList(hist, processName, false); // false tells it to add this one to the list of samples to be excluded form the legend
  return;
}

// this is a trivial overload to preserve back-compatibility with preexisting templates.
void WZplots::addHistToList(std::vector<TH1F*>* s, TH1F* h, std::string name, double scale)
{
  addHistToList(h,name);
  return;
}

void WZplots::addSignalToStack(TH1F* hist, std::string processName)
{
  hist = scaleXaxis(hist,xScale);
  int color = kWhite;
  std::string title = getProcessLegendEntry(processName);
  hist->SetTitle(title.c_str());
  hist->SetFillColor(color);
  histList.push_back(hist);
  return;
}

// this is a trivial overload to preserve back-compatibility with preesixsting templates.
void WZplots::addSignalToStack(std::vector<TH1F*>* s, TH1F* h, std::string name, double scale)
{
addSignalToStack(h,name);
return;
}

THStack* WZplots::buildStack(TPad* pad)
{
  THStack* stack = new THStack("stack","stack");
  for (unsigned int h=0;h<histList.size();h++)
    {
      TH1F* hist = histList.at(h);
      std::string addOption = "hist";
      stack->Add(hist,addOption.c_str());
    }
  pad->cd();
  stack->Draw(); // have to draw the stack in order to initialize it for further modification.
  stack->GetXaxis()->SetLimits(xmin, xmax);//andres

  if(yRangeUserDefined) stack->GetYaxis()->SetLimits(ymin,ymax);
  
  return stack;
}

TLegend* WZplots::buildLegend(Config config, TH1F* data, TH1F* sys, TGraphAsymmErrors* asym)
{

  const Double_t legend_x2 =legend_x1 + g_legend_width; //(g_legend_width+(legend->GetNColumns()-1)*legend->GetColumnSeparation())*legend->GetNColumns();
  const Double_t legend_y2 = legend_y1 + g_legend_height_per_entry +0.35;//g_legend_height_per_entry*(legend->GetNRows()+(header != "" ? nheaderlines : 0));

  TLegend * legend = new TLegend(legend_x1,legend_y1,legend_x2,legend_y2);
  // if(data) legend->AddEntry(data,"data 2015 (#sqrt{s}=13TeV)"); // old style
  const char *DataString="";
  if(configPlotsWZ.DataYears=="2017") DataString = "Data 2017";
  else if(configPlotsWZ.DataYears=="2015+2016") DataString = "Data 2015+16";
  else if(configPlotsWZ.DataYears=="2018") DataString = "Data 2018";
  else if(configPlotsWZ.DataYears=="Full") DataString = "Data 2015+16+17";

  if(data) legend->AddEntry(data,DataString); // new style, per Andreas Hoecker 
  if(sys) legend->AddEntry(sys,errorBandLegendEntry.c_str(),"f");
  if(asym) legend->AddEntry(asym,errorBandLegendEntry.c_str(),"f");
  if(errorBandLegendEntryTwoLines) legend->AddEntry((TObject*)0,errorBandLegendEntryLineTwo.c_str(),"");
  for (unsigned int h=histList.size();h>0;h--)
    {
      TH1F* hist = (TH1F*)(histList.at(h-1));
      std::vector<TH1F*>::iterator e = std::find(legendExcludes.begin(),legendExcludes.end(),hist); // search for this one in the list of hist.'s to exclude form the legend ... 
      if(e != legendExcludes.end()) continue;       // this one was on the list to exclude from the legend - move on.
      legend->AddEntry(hist, hist->GetTitle(),"f"); // "f" draws the fill color in a small box with no marker
    }
  return legend;
}
