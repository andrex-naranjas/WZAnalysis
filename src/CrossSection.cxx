//Class to calculate cross sections
#ifndef CROSSSECTION_CXX
#define CROSSSECTION_CXX

#include "CrossSection.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>


CrossSection::CrossSection()
{
}

CrossSection::~CrossSection(){}

void CrossSection::initialize(Config config){

  std::string dirYear="", year="";
  if(config.DataYears=="2015+2016"){ dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017"){ dirYear="2017/"; year="17";}
  if(config.DataYears=="2018"){ dirYear="2018/"; year="18";}
  if(config.DataYears=="Full"){ dirYear="Full/"; year="Full";}

  std::string dirInclusive, total;
  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}

  std::string puname="";
  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
  std::string systematic;
  systematic="_nominal";
 
  std::string wzchannel="", nameOfSample="";
  //zmumu
  if(config.WZSelection=="zmumu"){
    wzchannel="z"; nameOfSample="zmumu";
    fdata  = new TFile((config.OutputFileDir+"Files/"+ dirYear + "data" + year +"_z.root").c_str());
    fmc    = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fbg    = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + "backgrounds_" + nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fcw    = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
    fsys   = new TFile((config.OutputFileDir+"Sys/"+ dirYear + dirInclusive +"sys_"+wzchannel+total+".root").c_str());
    fout   = new TFile((config.OutputFileDir+"xSection/"+ dirYear + dirInclusive +"xSection_"+nameOfSample+total+".root").c_str(),"RECREATE");
  }

  //wplus
  if(config.WZSelection=="wplus" || config.WZSelection=="combined" || true){
    wzchannel="wplus"; nameOfSample="wplusmunu";
    fdata  = new TFile((config.OutputFileDir+"Files/"+ dirYear + "data" + year + "_wplus.root").c_str());
    fmc    = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fbg    = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + "backgrounds_"+ nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fcw    = new TFile((config.OutputFileDir+"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
    fsys   = new TFile((config.OutputFileDir+"Sys/"+ dirYear + dirInclusive +"sys_"+wzchannel + total + ".root").c_str());
    //multijet
    fmultijet = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/multijet_final.root").c_str());

    //fout   = new TFile((config.OutputFileDir+"xSection/"+ dirYear + dirInclusive +"xSection_"+nameOfSample+total+".root").c_str(),"RECREATE");
  }

  //wminus
  if(config.WZSelection=="wminus" || config.WZSelection=="combined" || true){
    wzchannel="wminus"; nameOfSample="wminmunu";
    //fdatap = new TFile((config.OutputFileDir+"Files/"+ dirYear + "data" + year +"_wminus.root").c_str());
    fdatap = new TFile((config.OutputFileDir+"Files/"+ dirYear + "data" + year +"_wminus.root").c_str());
    fmcp   = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fbgp   = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + "backgrounds_" + nameOfSample + "_" + wzchannel + systematic + puname +total+".root").c_str());
    fcwp   = new TFile((config.OutputFileDir+"Files/"+ dirYear +"Truth/"+ dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str());
    fsysp  = new TFile((config.OutputFileDir+"Sys/"+ dirYear + dirInclusive +"sys_"+wzchannel + total +".root").c_str());
    fout   = new TFile((config.OutputFileDir+"xSection/"+ dirYear + dirInclusive +"xSection"+total+".root").c_str(),"RECREATE");
  }                


  //check if files are open/exist
  if(!fdata->IsOpen() || !fmc->IsOpen() || !fbg->IsOpen() || !fcw->IsOpen() || !fsys->IsOpen()){
    std::cout<<"Cross section option on, but files not present, bye!"<<std::endl;
    exit(1);}

  if(!fdatap->IsOpen() || !fmcp->IsOpen() || !fbgp->IsOpen() || !fcwp->IsOpen() || !fsysp->IsOpen()){
    std::cout<<"Cross section option on, but files not present, bye!"<<std::endl;
    exit(1);}


  //debug
  scale  = 1.065;//w+
  scalep = 1.045;//w-
  std::cout<<"******    Cross Section Calculation    ******"<<std::endl;
  return;
}


void CrossSection::execute(Config config){

  InclusiveXsec(config);
  EtaBinsXsection(config);
  
  setstyle();
  xSecPlot(config,hEtaBinXsec,hEtaBinXsecMC,"wplus",200.,3000.);
  xSecPlot(config,hEtaBinXsecP,hEtaBinXsecMCP,"wminus",100.,2500.);
  xSecPlot(config,hEtaBinAsymm,hEtaBinAsymmMC,"asymmetry",0.,.3);

  return;
}

void CrossSection::InclusiveXsec(Config config){

  TH1D *hData, *hMC, *hBG, *hReco, *hTruth, *hSysUp, *hSysDown;
  TH1D *hDataP, *hMCP, *hBGP, *hRecoP, *hTruthP, *hSysUpP, *hSysDownP;
  TH1D *hMultijet, *hMultijetP;
  double nData=0., nMC=0., nBG=0., nReco=0., nTruth=0., Cwz=0.;
  double nDataP=0., nMCP=0., nBGP=0., nRecoP=0., nTruthP=0., CwzP=0.;
  double nMultijet=0., nMultijetP=0.;
  TH1D *hStatError, *hStatErrorP;
  double statErrorCw=0., statErrorCwP;

  if(config.WZSelection=="wplus" || config.WZSelection=="combined" || true){
    hData    = (TH1D*)fdata->GetObjectUnchecked("h_rmass");
    hMC      = (TH1D*)fmc  ->GetObjectUnchecked("h_rmass");
    hBG      = (TH1D*)fbg  ->GetObjectUnchecked("h_rmass");
    hReco    = (TH1D*)fcw  ->GetObjectUnchecked("hMassCrec");
    hTruth   = (TH1D*)fcw  ->GetObjectUnchecked("hMassCgen");
    hSysUp   = (TH1D*)fsys ->GetObjectUnchecked("h_upInd_Total");
    hSysDown = (TH1D*)fsys ->GetObjectUnchecked("h_downInd_Total");
    hMultijet= (TH1D*)fmultijet ->GetObjectUnchecked("eta_wplus");
    hStatError=(TH1D*)fcw ->GetObjectUnchecked("hCwzErrorsNew");

    nData=hData->Integral(0,hData->GetNbinsX()+1);
    nMC=hMC->Integral(0,hMC->GetNbinsX()+1);
    nBG=hBG->Integral(0,hBG->GetNbinsX()+1);
    nReco=hReco->Integral(0,hReco->GetNbinsX()+1);
    nTruth=hTruth->Integral(0,hTruth->GetNbinsX()+1);
    nMultijet=hMultijet->Integral(0,hMultijet->GetNbinsX()+1);
    nBG=nBG+nMultijet;//include the multijet bg
  }
  
  if(config.WZSelection=="wminus" || true){
    hDataP     = (TH1D*)fdatap->GetObjectUnchecked("h_rmass");
    hMCP       = (TH1D*)fmcp  ->GetObjectUnchecked("h_rmass");
    hBGP       = (TH1D*)fbgp  ->GetObjectUnchecked("h_rmass");
    hRecoP     = (TH1D*)fcwp  ->GetObjectUnchecked("hMassCrec");
    hTruthP    = (TH1D*)fcwp  ->GetObjectUnchecked("hMassCgen");
    hSysUpP    = (TH1D*)fsysp ->GetObjectUnchecked("h_upInd_Total");
    hSysDownP  = (TH1D*)fsysp ->GetObjectUnchecked("h_downInd_Total");
    hMultijetP = (TH1D*)fmultijet->GetObjectUnchecked("eta_wminus");
    hStatErrorP= (TH1D*)fcwp  ->GetObjectUnchecked("hCwzErrorsNew");

    nDataP=hDataP->Integral(0,hDataP->GetNbinsX()+1);
    nMCP=hMCP->Integral(0,hMCP->GetNbinsX()+1);
    nBGP=hBGP->Integral(0,hBGP->GetNbinsX()+1);
    nRecoP=hRecoP->Integral(0,hRecoP->GetNbinsX()+1);
    nTruthP=hTruthP->Integral(0,hTruthP->GetNbinsX()+1);
    nMultijetP=hMultijetP->Integral(0,hMultijetP->GetNbinsX()+1);
    nBGP=nBGP+nMultijetP;//include the multijet bg
  }

  if(config.WZSelection=="combined"){
    hData->Add(hDataP); hMC->Add(hMCP); hBG->Add(hBGP);
    hReco->Add(hRecoP); hTruth->Add(hTruthP);
    hSysUp->Add(hSysUpP); hSysDown->Add(hSysDownP);
    nData=nData+nDataP; nMC=nMC+nMCP; nBG=nBG+nBGP;
    nReco=nReco+nRecoP; nTruth=nTruth+nTruthP;
  }

  std::vector<std::string> systematic, header;
  std::vector<double> value;
  std::vector<std::vector<double>> totalValues; totalValues.clear();
  
  double sysSumUp=0., sysSumDn=0., sysSumUpP=0., sysSumDnP=0.;

  systematic=Names("systematics");

  //statistical systematics (Cw factor)
  for(int k = 1; k<hStatError->GetXaxis()->GetNbins()+1;k++)
    statErrorCw= statErrorCw + std::pow(hStatError->GetBinContent(k),2);  
  statErrorCw=std::sqrt(statErrorCw);

  for(int k = 1; k<hSysUp->GetXaxis()->GetNbins()+1;k++){
    value.push_back(hSysUp->GetBinContent(k)*100.);
    sysSumUp=sysSumUp+hSysUp->GetBinContent(k)*hSysUp->GetBinContent(k);}
  sysSumUp+=std::pow(statErrorCw,2);//add in quadrature cw stat error
    //systematic.push_back(hSysUp->GetXaxis()->GetBinLabel(k));}
  value.push_back(statErrorCw*100.);//add stat error
  totalValues.push_back(value);
  value.clear();
  header.push_back("$W^{+} (up)$");

  for(int k = 1; k<hSysDown->GetXaxis()->GetNbins()+1;k++){
    value.push_back(hSysDown->GetBinContent(k)*100.);
    sysSumDn=sysSumDn+hSysDown->GetBinContent(k)*hSysDown->GetBinContent(k);}
  sysSumDn+=std::pow(statErrorCw,2);//add in quadrature cw stat error

  value.push_back(statErrorCw*100.);//add stat error
  totalValues.push_back(value);
  value.clear();
  header.push_back("$W^{+} (down)$");

  //negative channel
  //statistical systematics (Cw factor)
  for(int k = 1; k<hStatErrorP->GetXaxis()->GetNbins()+1;k++)
    statErrorCwP= statErrorCwP + std::pow(hStatErrorP->GetBinContent(k),2);  
  statErrorCwP=std::sqrt(statErrorCwP);

  for(int k = 1; k<hSysUpP->GetXaxis()->GetNbins()+1;k++){
    value.push_back(hSysUpP->GetBinContent(k)*100.);
    sysSumUpP=sysSumUpP+hSysUpP->GetBinContent(k)*hSysUpP->GetBinContent(k);}
  sysSumUpP+=std::pow(statErrorCwP,2);//add in quadrature cw stat error

  value.push_back(statErrorCwP*100.);//add stat error
  totalValues.push_back(value);
  value.clear();
  header.push_back("$W^{-} (up)$");

  for(int k = 1; k<hSysDownP->GetXaxis()->GetNbins()+1;k++){
    value.push_back(hSysDownP->GetBinContent(k)*100.);
    sysSumDnP=sysSumDnP+hSysDownP->GetBinContent(k)*hSysDownP->GetBinContent(k);}
  sysSumDnP+=std::pow(statErrorCwP,2);//add in quadrature cw stat error

  value.push_back(statErrorCwP*100.);//add stat error
  totalValues.push_back(value);
  value.clear();
  header.push_back("$W^{-} (down)$");
      
  sysTable(totalValues, systematic, header,"inclusive_sys");

  std::cout<<"nData:    "<<nData <<"    nMC:   "<<nMC<<"     nBG:   "<<nBG<<"      nReco:   "<<nReco<<"     nTruth:    "<<nTruth<<std::endl;

  nReco=nMC*scale;
  nRecoP=nMCP*scalep;

  Cwz=nReco/nTruth;   CwzP=nRecoP/nTruthP;
  double lumiData=getDataLumi(config);
  double xSection=0., xSectionMC=0., xSectionDen=0., xSectionNum=0., xSectionNumMC=0.;
  double xSectionP=0., xSectionMCP=0., xSectionDenP=0., xSectionNumP=0., xSectionNumMCP=0.;

  xSectionNum   = nData-nBG*scale;
  xSectionNumMC = nReco;
  xSectionDen   = lumiData * Cwz;
  xSection      = xSectionNum/xSectionDen;
  xSectionMC    = xSectionNumMC/xSectionDen;

  xSectionNumP   = nDataP-nBGP*scalep;
  xSectionNumMCP = nRecoP;
  xSectionDenP   = lumiData * CwzP;
  xSectionP      = xSectionNumP/xSectionDenP;
  xSectionMCP    = xSectionNumMCP/xSectionDenP;

  std::cout<<"xSection:  "<<xSection<<"     "<<"Cwz    "<<Cwz<<"     "<<"nData:    "<<nData<<"      nReco:    "<<nReco<<"     nTruth:    "<<nTruth<<"      nBG:    "<<nBG<<"      lumiData:    "<<lumiData<<std::endl;

  double statError =std::sqrt(nData)/( Cwz * 44307.4);
  double statErrorP=std::sqrt(nDataP)/( CwzP * 44307.4);



  std::vector<std::string> boson; boson.clear(); boson.push_back("$\\sigma_{W^{+}}$[nb]"); boson.push_back("$\\sigma_{W^{-}}$[nb]"); boson.push_back("$As_{\\mu}$");
  header.clear(); totalValues.clear();//recycle vectors, be green!
  header=Names("xSection_inc");
  std::vector<double> wplus_val, wminus_val, asymmetry;
  wplus_val.clear(); wminus_val.clear();
  //push_back values in the correct order (order matters)
  wplus_val.push_back(xSection*.001); wplus_val.push_back(xSectionMC*.001); wplus_val.push_back(Cwz); wplus_val.push_back(std::sqrt(sysSumUp)*100.); wplus_val.push_back(std::sqrt(sysSumDn)*100.); wplus_val.push_back((statError*100)/xSection); wplus_val.push_back(nData); wplus_val.push_back(nReco); wplus_val.push_back(nTruth); wplus_val.push_back(nBG);

  wminus_val.push_back(xSectionP*.001); wminus_val.push_back(xSectionMCP*.001); wminus_val.push_back(CwzP); wminus_val.push_back(std::sqrt(sysSumUpP)*100.); wminus_val.push_back(std::sqrt(sysSumDnP)*100.); wminus_val.push_back((statErrorP*100)/xSection); wminus_val.push_back(nDataP); wminus_val.push_back(nRecoP); wminus_val.push_back(nTruthP); wminus_val.push_back(nBGP);

  //asymmetry
  double val_asy   =(xSection-xSectionP)/(xSection+xSectionP);
  double val_asyMC =(xSectionMC-xSectionMCP)/(xSectionMC+xSectionMCP);
  double dat_asy   =(nData-nDataP)/(nData+nDataP);
  double rec_asy   =(nReco-nRecoP)/(nReco+nRecoP);
  double tru_asy   =(nTruth-nTruthP)/(nTruth+nTruthP);
  double bac_asy   =(nBG-nBGP)/(nBG+nBGP);
  double up_asy    =asymError(xSection,xSectionP,std::sqrt(sysSumUp)*xSection,std::sqrt(sysSumUpP)*xSectionP);
  double dn_asy    =asymError(xSection,xSectionP,std::sqrt(sysSumDn)*xSection,std::sqrt(sysSumDnP)*xSectionP);
  double stat_asy  =asymError(xSection,xSectionP,statError,statErrorP);

  asymmetry.push_back(val_asy); asymmetry.push_back(val_asyMC); asymmetry.push_back(0); asymmetry.push_back(up_asy*100.); asymmetry.push_back(dn_asy*100.); asymmetry.push_back(stat_asy*100.);
  asymmetry.push_back(0); asymmetry.push_back(0); asymmetry.push_back(0); asymmetry.push_back(0);

  totalValues.push_back(wplus_val); totalValues.push_back(wminus_val); totalValues.push_back(asymmetry);

  std::vector<std::string> eta_bin_names=LatexBinName();
  xSecTableIncl(config,totalValues,boson,header,"xSection");

  return;
}


void CrossSection::EtaBinsXsection(Config config){
  
  TH1D *hData, *hMC, *hBG, *hReco, *hTruth, *hEtaCsimple;
  TH1D *hDataP, *hMCP, *hBGP, *hRecoP, *hTruthP, *hEtaCsimpleP;
  TH1D *hMultijet, *hMultijetP;
  double nData=0., nMC=0., nBG=0., nReco=0., nTruth=0., Cwz=0.;
  double nDataP=0., nMCP=0., nBGP=0., nRecoP=0., nTruthP=0., CwzP=0.;
  TH2D *hEtaSysUp, *hEtaSysUpP, *hEtaSysDown, *hEtaSysDownP;
  double nMultijet=0., nMultijetP=0.;
  TH1D *hStatError, *hStatErrorP;
  double statError=0., statErrorP;

  if(config.WZSelection=="wplus" || config.WZSelection=="combined"  || true){
    hData  = (TH1D*)fdata->GetObjectUnchecked("h_eta");
    hMC    = (TH1D*)fmc  ->GetObjectUnchecked("h_eta");
    hBG    = (TH1D*)fbg  ->GetObjectUnchecked("h_eta");
    hReco  = (TH1D*)fcw  ->GetObjectUnchecked("hEtaCrec");
    hTruth = (TH1D*)fcw  ->GetObjectUnchecked("hEtaCgen");
    //multijet
    hMultijet= (TH1D*)fmultijet ->GetObjectUnchecked("eta_wplus");
    //stat error
    hStatError=(TH1D*)fcw ->GetObjectUnchecked("hCwzErrorsNew");
    
    hEtaSysUp   = (TH2D*)fsys ->GetObjectUnchecked("h_upInd_eta");
    hEtaSysDown = (TH2D*)fsys ->GetObjectUnchecked("h_downInd_eta");
    hEtaCsimple = (TH1D*)fcw  ->GetObjectUnchecked("hCw_eta_simple");
    
    nData=hData->Integral(0,hData->GetNbinsX()+1);
    nMC=hMC->Integral(0,hMC->GetNbinsX()+1);
    nBG=hBG->Integral(0,hBG->GetNbinsX()+1);
    nReco=hReco->Integral(0,hReco->GetNbinsX()+1);
    nTruth=hTruth->Integral(0,hTruth->GetNbinsX()+1);
    nMultijet=hMultijet->Integral(0,hMultijet->GetNbinsX()+1);
    nBG=nBG+nMultijet;//include the multijet bg
  }
  
  if(config.WZSelection=="wminus" || true){
    hDataP  = (TH1D*)fdatap->GetObjectUnchecked("h_eta");
    hMCP    = (TH1D*)fmcp  ->GetObjectUnchecked("h_eta");
    hBGP    = (TH1D*)fbgp  ->GetObjectUnchecked("h_eta");
    hRecoP  = (TH1D*)fcwp  ->GetObjectUnchecked("hEtaCrec");
    hTruthP = (TH1D*)fcwp  ->GetObjectUnchecked("hEtaCgen");
    //multijet
    hMultijetP= (TH1D*)fmultijet->GetObjectUnchecked("eta_wminus");
    //stat error
    hStatErrorP=(TH1D*)fcwp ->GetObjectUnchecked("hCwzErrorsNew");

    hEtaSysUpP   = (TH2D*)fsysp ->GetObjectUnchecked("h_upInd_eta");
    hEtaSysDownP = (TH2D*)fsysp ->GetObjectUnchecked("h_downInd_eta");
    hEtaCsimpleP = (TH1D*)fcwp  ->GetObjectUnchecked("hCw_eta_simple");

    nDataP=hDataP->Integral(0,hDataP->GetNbinsX()+1);
    nMCP=hMCP->Integral(0,hMCP->GetNbinsX()+1);
    nBGP=hBGP->Integral(0,hBGP->GetNbinsX()+1);
    nRecoP=hRecoP->Integral(0,hRecoP->GetNbinsX()+1);
    nTruthP=hTruthP->Integral(0,hTruthP->GetNbinsX()+1);
    nMultijetP=hMultijetP->Integral(0,hMultijetP->GetNbinsX()+1);
    nBGP=nBGP+nMultijetP;//include the multijet bg
  }

  if(config.WZSelection=="combined"){
    hData->Add(hDataP); hMC->Add(hMCP); hBG->Add(hBGP);
    hReco->Add(hRecoP); hTruth->Add(hTruthP);
    hEtaCsimple ->Add(hEtaCsimpleP);
    nData=nData+nDataP; nMC=nMC+nMCP; nBG=nBG+nBGP;
    nReco=nReco+nRecoP; nTruth=nTruth+nTruthP;
  }
  
  double lumiData=getDataLumi(config);

  unsigned int nx = hData->GetXaxis()->GetNbins();
  double* xbins = new double[hData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hData->GetXaxis()->GetNbins()]=hData->GetXaxis()->GetBinUpEdge(hData->GetXaxis()->GetNbins());

  hEtaBinXsec     = new TH1D("hEtaBinXsec"   ,"xSection",nx,xbins);
  hEtaBinXsecP    = new TH1D("hEtaBinXsecP"  ,"xSection",nx,xbins);
  hEtaBinAsymm    = new TH1D("hEtaBinAsymm"  ,"xSection",nx,xbins);

  hEtaBinXsecMC   = new TH1D("hEtaBinXsecMC"   ,"xSection",nx,xbins);
  hEtaBinXsecMCP  = new TH1D("hEtaBinXsecMCP"  ,"xSection",nx,xbins);
  hEtaBinAsymmMC  = new TH1D("hEtaBinAsymmMC"  ,"xSection",nx,xbins);

  //systable
  std::vector<std::string> systematic, header;  systematic=Names("systematics");
  std::vector<std::vector<double>> totalValues, totalValuesP;
  totalValues.clear(); totalValuesP.clear();

  //xSection table
  std::vector<std::string> eta_bin; eta_bin.clear();
  std::vector<std::string> header_eta; header_eta.clear();
  std::vector<std::vector<double>> xSec_wplusEta, xSec_wminEta, asymmetry;
  xSec_wplusEta.clear(); xSec_wminEta.clear(); asymmetry.clear();
  header_eta=Names("xSection");


  //summary
  std::vector<std::vector<double>> xSec_summary;      xSec_summary.clear();
  std::vector<std::vector<double>> xSec_summary_up;   xSec_summary_up.clear();
  std::vector<std::vector<double>> xSec_summary_down; xSec_summary_down.clear();
  std::vector<std::vector<double>> xSec_summary_stat; xSec_summary_stat.clear();

  std::vector<double> summary_wplus, summary_wminus, summary_asy; summary_wplus.clear(); summary_wminus.clear(); summary_asy.clear(); 
  std::vector<double> summary_wplus_up, summary_wminus_up, summary_asy_up; summary_wplus_up.clear(); summary_wminus_up.clear(); summary_asy_up.clear();
  std::vector<double> summary_wplus_down, summary_wminus_down, summary_asy_down;summary_wplus_down.clear(); summary_wminus_down.clear(); summary_asy_down.clear();
  std::vector<double> summary_wplus_stat, summary_wminus_stat, summary_asy_stat; summary_wplus_stat.clear(); summary_wminus_stat.clear(); summary_asy_stat.clear();

  
  for(int k=1;k<hData->GetXaxis()->GetNbins()+1;k++){//every eta bin
    
    //systable
    std::vector<double> valueUp; valueUp.clear();
    std::vector<double> valueDn; valueDn.clear();
    header.push_back("Up"); header.push_back("Down");

    double sysSumUp=0., sysSumDn=0., sysSumUpP=0., sysSumDnP=0.;

    std::vector<std::vector<double>> single_bin_sys; single_bin_sys.clear();

    double sumSys=0.;
    for(int l=1; l<hEtaSysDown->GetYaxis()->GetNbins()+1;l++){
      double sysTemp = std::pow((std::fabs(hEtaSysDown->GetBinContent(k,l))+std::fabs(hEtaSysUp->GetBinContent(k,l)))*0.5,2);
      valueUp.push_back((hEtaSysDown->GetBinContent(k,l))*100.); valueDn.push_back((hEtaSysUp->GetBinContent(k,l))*100.);
      sysSumUp=sysSumUp+hEtaSysUp->GetBinContent(k,l)*hEtaSysUp->GetBinContent(k,l);
      sysSumDn=sysSumDn+hEtaSysDown->GetBinContent(k,l)*hEtaSysDown->GetBinContent(k,l);
      sumSys=sumSys+sysTemp;      
    }
    totalValues.push_back(valueUp);   totalValues.push_back(valueDn);
    valueUp.clear(); valueDn.clear();

    double sumSysP=0.;
    for(int l=1; l<hEtaSysDownP->GetYaxis()->GetNbins()+1;l++){
      double sysTempP = std::pow((std::fabs(hEtaSysDownP->GetBinContent(k,l))+std::fabs(hEtaSysUpP->GetBinContent(k,l)))*0.5,2);
      valueUp.push_back((hEtaSysDownP->GetBinContent(k,l))*100.); valueDn.push_back((hEtaSysUpP->GetBinContent(k,l)*100.));
      sysSumUpP=sysSumUp+hEtaSysUpP->GetBinContent(k,l)*hEtaSysUpP->GetBinContent(k,l);
      sysSumDnP=sysSumDn+hEtaSysDownP->GetBinContent(k,l)*hEtaSysDownP->GetBinContent(k,l);
      sumSysP=sumSysP+sysTempP;
    }
    totalValuesP.push_back(valueUp);  totalValuesP.push_back(valueDn);
    valueUp.clear(); valueDn.clear();

    //Cw factor MC statistics
    double statErrorCw  = hStatError->GetBinContent(k);
    double statErrorCwP = hStatErrorP->GetBinContent(k);

    sysSumUp+=std::pow(statErrorCw,2);	
    sysSumDn+=std::pow(statErrorCwP,2);
    
    sysSumUpP+=std::pow(statErrorCw,2);
    sysSumDnP+=std::pow(statErrorCwP,2);

    //single eta bin systematic
    for(int l=1; l<hEtaSysDown->GetYaxis()->GetNbins()+1;l++){
      valueUp.push_back((hEtaSysDown->GetBinContent(k,l))*100.);
      valueDn.push_back((hEtaSysUp->GetBinContent(k,l)*100.));
    }
    single_bin_sys.push_back(valueUp);  single_bin_sys.push_back(valueDn);
    valueUp.clear(); valueDn.clear();

    for(int l=1; l<hEtaSysDownP->GetYaxis()->GetNbins()+1;l++){
      valueUp.push_back((hEtaSysDownP->GetBinContent(k,l))*100.);
      valueDn.push_back((hEtaSysUpP->GetBinContent(k,l)*100.));
    }
    single_bin_sys.push_back(valueUp);  single_bin_sys.push_back(valueDn);
    valueUp.clear(); valueDn.clear();

    std::vector<std::string> systematic, header; systematic.clear(); header.clear();
    systematic = Names("systematics");
    header.push_back("$W^{+} (up)$");    header.push_back("$W^{+} (down)$");
    header.push_back("$W^{-} (up)$");    header.push_back("$W^{-} (down)$");
    //systematic table for every bin
    std::string bin_eta_string; bin_eta_string = "eta_"+std::to_string(k)+"_sys";
    sysTable(single_bin_sys, systematic, header,bin_eta_string);


    //xSection histos
    double width = hData->GetBinWidth(k); if(k==hData->GetXaxis()->GetNbins()) width = 2.4-2.18;
    double data  = hData->GetBinContent(k);
    double back  = hBG->GetBinContent(k)*scale;
    double multi = hMultijet->GetBinContent(k)*scale;
    back = back + multi;//include multijet
    double reco  = hMC->GetBinContent(k)*scale;
    double truth = hTruth->GetBinContent(k);
    //double unfo  = hEtaCsimple->GetBinContent(k)*scale;
    double unfo=reco/truth;
    double xSection = (data-back)/(width * lumiData * unfo);
    double xSectionMC = reco/(width * lumiData * unfo);
    sumSys=std::sqrt(sumSys);//*0.01*xSection;
    hEtaBinXsec->SetBinError(k,sumSys*xSection);
    hEtaBinXsec->SetBinContent(k,xSection);
    hEtaBinXsecMC->SetBinError(k,sumSys*xSectionMC);
    hEtaBinXsecMC->SetBinContent(k,xSectionMC);

    std::cout<<"test   "<<data<<"   "<<back<<"   "<<lumiData<<"    "<<unfo<<"    "<<xSection<<std::endl;

    double dataP  = hDataP->GetBinContent(k);
    double backP  = hBGP->GetBinContent(k)*scalep;
    double multiP = hMultijetP->GetBinContent(k)*scalep;
    backP=backP+multiP;//include multijet
    double recoP  = hMCP->GetBinContent(k)*scalep;
    double truthP = hTruthP->GetBinContent(k);
    //double unfoP  = hEtaCsimpleP->GetBinContent(k)*scale;
    double unfoP=recoP/truthP;
    double xSectionP = (dataP-backP)/(width * lumiData * unfoP);
    double xSectionMCP = recoP/(width * lumiData * unfoP);

    std::cout<<recoP/truthP<<"  Parrito   "<<unfoP<<std::endl;

    sumSysP=std::sqrt(sumSysP);//*0.01*xSectionP;//what's this?
    hEtaBinXsecP->SetBinError(k,sumSysP*xSectionP);
    hEtaBinXsecP->SetBinContent(k,xSectionP);

    hEtaBinXsecMCP->SetBinError(k,sumSys*xSectionMCP);
    hEtaBinXsecMCP->SetBinContent(k,xSectionMCP);


    //data statistical error
    double statError  =std::sqrt(nData)/(unfo * 44307.4);
    double statErrorP =std::sqrt(nDataP)/(unfoP * 44307.4);

    //xSection table
    eta_bin.push_back(hData->GetXaxis()->GetBinLabel(k));
    std::vector<double> wplus_val, wminus_val, asymm_val;
    wplus_val.clear(); wminus_val.clear(); asymm_val.clear();
    //push_back values in the correct order (order really matters)
    wplus_val.push_back(xSection*.001); wplus_val.push_back(xSectionMC*.001); wplus_val.push_back(unfo); wplus_val.push_back(std::sqrt(sysSumUp)*100.); wplus_val.push_back(std::sqrt(sysSumDn)*100.); wplus_val.push_back((statError*100)/xSection); wplus_val.push_back(data); wplus_val.push_back(reco); wplus_val.push_back(truth); wplus_val.push_back(back);
    
    wminus_val.push_back(xSectionP*.001); wminus_val.push_back(xSectionMCP*.001); wminus_val.push_back(unfoP); wminus_val.push_back(std::sqrt(sysSumUpP)*100.); wminus_val.push_back(std::sqrt(sysSumDnP)*100.); wminus_val.push_back((statErrorP*100)/xSectionP); wminus_val.push_back(dataP); wminus_val.push_back(recoP); wminus_val.push_back(truthP); wminus_val.push_back(backP);

    //*******asymmetry*******
    double val_asy=(xSection-xSectionP)/(xSection+xSectionP);
    double val_asyMC=(xSectionMC-xSectionMCP)/(xSectionMC+xSectionMCP);
    double dat_asy=(data-dataP)/(data+dataP);
    double rec_asy=(reco-recoP)/(reco+recoP);
    double tru_asy=(truth-truthP)/(truth+truthP);
    double bac_asy=(back-backP)/(back+backP);
    double up_asy = asymError(xSectionMC,xSectionMCP,std::sqrt(sysSumUp)*xSectionMC,std::sqrt(sysSumUpP)*xSectionMCP); //up_asy=up_asy/val_asyMC;
    double dn_asy = asymError(xSectionMC,xSectionMCP,std::sqrt(sysSumDn)*xSectionMC,std::sqrt(sysSumDnP)*xSectionMCP); //dn_asy=dn_asy/val_asyMC;
    double stat_as= asymError(xSectionMC,xSectionMCP,statError,statErrorP);

    hEtaBinAsymm->SetBinError(k,dn_asy*val_asy);//(up_asy+dn_asy)*0.5);
    hEtaBinAsymm->SetBinContent(k,val_asy);

    hEtaBinAsymmMC->SetBinError(k,dn_asy*val_asyMC);//(up_asy+dn_asy)*0.5);
    hEtaBinAsymmMC->SetBinContent(k,val_asyMC);
 
    asymm_val.push_back(val_asy); asymm_val.push_back(val_asyMC); asymm_val.push_back(up_asy*100.); asymm_val.push_back(dn_asy*100.);
    asymm_val.push_back((stat_as*100));

    xSec_wplusEta.push_back(wplus_val); xSec_wminEta.push_back(wminus_val); asymmetry.push_back(asymm_val);

    //summary
    //push_back values in the correct order (order really matters) summary
    summary_wplus.push_back(xSection*.001);        summary_wplus_up.push_back(std::sqrt(sysSumUp));  
    summary_wminus.push_back(xSectionP*.001);      summary_wminus_up.push_back(std::sqrt(sysSumUpP));
    summary_asy.push_back(val_asy);                summary_asy_up.push_back(up_asy*100);          
    
    summary_wplus_down.push_back(std::sqrt(sysSumDn));   summary_wplus_stat.push_back(statError);  
    summary_wminus_down.push_back(std::sqrt(sysSumDnP)); summary_wminus_stat.push_back(statErrorP);
    summary_asy_down.push_back(dn_asy*100);                  summary_asy_stat.push_back(up_asy*100);
    
    
  }//eta bins loop!

    xSec_summary.push_back(summary_wminus);        xSec_summary_up.push_back(summary_wminus_up); 
    xSec_summary.push_back(summary_wplus);         xSec_summary_up.push_back(summary_wplus_up);  
    xSec_summary.push_back(summary_asy);           xSec_summary_up.push_back(summary_asy_up);    

    xSec_summary_down.push_back(summary_wminus_down);   xSec_summary_stat.push_back(summary_wminus_stat); 
    xSec_summary_down.push_back(summary_wplus_down);    xSec_summary_stat.push_back(summary_wplus_stat);  
    xSec_summary_down.push_back(summary_asy_down);      xSec_summary_stat.push_back(summary_asy_stat);    
 
     
    sysTable(totalValues,  systematic, header,"eta_wplus");
    sysTable(totalValuesP, systematic, header,"eta_wminus");
  
  std::vector<std::string> eta_bin_names=LatexBinName();

  xSecTable(config,xSec_wplusEta,eta_bin_names,header_eta,"xSection_wplusEta");
  xSecTable(config,xSec_wminEta,eta_bin_names,header_eta,"xSection_wminEta");

  header_eta.clear(); header_eta = Names("Asymmetry");
  asyTable(config,asymmetry,eta_bin_names,header_eta,"xSection_Asymmetry");

  //summary results
  std::vector<std::string> header_eta_sum; header_eta_sum.clear(); header_eta_sum =Names("summary");
  xSecTableSum(xSec_summary,xSec_summary_stat,xSec_summary_up,xSec_summary_down,eta_bin_names,header_eta_sum,"xSection_Eta_summary");

  nReco=nMC;
  Cwz=nReco/nTruth;

  double xSection=0., xSectionDen=0., xSectionNum=0. ;

  xSectionNum = nData-nBG;
  xSectionDen = lumiData * Cwz;
  xSection = xSectionNum/xSectionDen;

  return;
}


double CrossSection::getDataLumi(Config config){

  //lumi normalisation
  double lumiData=0.;
  if(config.DataYears=="2015"){
    lumiData=3219.56;
  }else if(config.DataYears=="2017"){
    lumiData=44307.4;
  }else if(config.DataYears=="2018"){
    lumiData=59937.2;
  }else if(config.DataYears=="2015+2016"){
    lumiData=  32988.1+3219.56;}

  return lumiData;
}


double CrossSection::asymError(double xSec, double xSecP, double error, double errorP){

  double firstFac  = 2. / (std::pow((xSec+xSecP),2));
  double secondFac = std::pow(xSecP,2)*std::pow(error,2) + std::pow(xSec,2)*std::pow(errorP,2);
  secondFac = std::sqrt(secondFac);

  std::cout<<"    hermos parrito"<<firstFac*secondFac<<std::endl;

  return firstFac*secondFac;
}


void CrossSection::sysTable(std::vector<std::vector<double>> totalValues, std::vector<std::string> systematic, std::vector<std::string> header, std::string fileName){

  if(header.size()!=totalValues.size()) return;

  std::string bin_label="";
  if(fileName=="inclusively_sys") bin_label=" inclusive";
  else if(fileName=="eta_1_sys")  bin_label=" $0.00\\leq\\vert\\eta^{\\mu}\\vert<0.21$";
  else if(fileName=="eta_2_sys")  bin_label=" $0.21\\leq\\vert\\eta^{\\mu}\\vert<0.42$";
  else if(fileName=="eta_3_sys")  bin_label=" $0.42\\leq\\vert\\eta^{\\mu}\\vert<0.63$";
  else if(fileName=="eta_4_sys")  bin_label=" $0.63\\leq\\vert\\eta^{\\mu}\\vert<0.84$";
  else if(fileName=="eta_5_sys")  bin_label=" $0.84\\leq\\vert\\eta^{\\mu}\\vert<1.05$";
  else if(fileName=="eta_6_sys")  bin_label=" $1.05\\leq\\vert\\eta^{\\mu}\\vert<1.37$";
  else if(fileName=="eta_7_sys")  bin_label=" $1.37\\leq\\vert\\eta^{\\mu}\\vert<1.52$";
  else if(fileName=="eta_8_sys")  bin_label=" $1.52\\leq\\vert\\eta^{\\mu}\\vert<1.74$";
  else if(fileName=="eta_9_sys")  bin_label=" $1.74\\leq\\vert\\eta^{\\mu}\\vert<1.95$";
  else if(fileName=="eta_10_sys") bin_label=" $1.95\\leq\\vert\\eta^{\\mu}\\vert<2.18$";
  else if(fileName=="eta_11_sys") bin_label=" $2.18\\leq\\vert\\eta^{\\mu}\\vert<2.40$";

  int nHeader = (int)header.size(); int nValues = (int)totalValues.at(1).size();
  std::ofstream sysLatex(("./"+fileName+".tex").c_str());

  sysLatex << "\\begin{table}[h!]" << std::endl;
  sysLatex << "\\begin{center}" << std::endl;
  sysLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) sysLatex<<" c "; 
  sysLatex<<" }\\hline \\hline " << std::endl;

  sysLatex<< "Systematic & ";
  for(int k=0; k < nHeader; k++){
    sysLatex<<std::setw(6)<<std::right;
    if(k< nHeader-1) sysLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) sysLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nValues; j++){
    sysLatex<<systematic.at(j)<<" &";
    for(int l=0;l < nHeader; l++){
      sysLatex<<std::setw(6)<<std::right;
      //if(l< nHeader-1) sysLatex<<printf("%.2f",(totalValues.at(l)).at(j))<<" & ";
      if(std::fabs((totalValues.at(l)).at(j)) < 0.000001) (totalValues.at(l)).at(j) =0.;
      if(l< nHeader-1) sysLatex<<"\\num{"<<std::fabs((totalValues.at(l)).at(j))<<"} & ";
      if(l==nHeader-1) sysLatex<<"\\num{"<<std::fabs((totalValues.at(l)).at(j))<<"} \\\\"<<std::endl;
    }
  }
  
  sysLatex << "\\hline \\hline" << std::endl;
  sysLatex << "\\end{tabular}" << std::endl;
  sysLatex << "\\end{center}" << std::endl;
sysLatex << "\\caption{This table contains the contribution of each individual uncertainty in the "<<bin_label<<" bin, following the definition given in \\cref{eq:indtotsys}. The first column indicates the name of the systematic; the second and third columns, corresponding to the $W^{+}\\rightarrow\\mu^{+}\\nu$, are the up and down variation respectively; the fourth and fifth columns, corresponding to the $W^{-}\\rightarrow\\mu^{-}\\nu$, are the up and down variation respectively. In this table the symmetric uncertainty definition is not considered. All the values are presented in percentage.}" << std::endl;  
  sysLatex << "\\label{tab:inclxSec}" << std::endl;
  sysLatex << "\\end{table}" << std::endl;
  sysLatex.close();

  return;
}

void CrossSection::xSecTableIncl(Config config, std::vector<std::vector<double>> totalValues, std::vector<std::string> column, std::vector<std::string> header, std::string fileName){

  //if(header.size()!=totalValues.at(1).size()) return;

  int nHeader = (int)header.size(); int nValues = (int)totalValues.size();
  std::ofstream xSecLatex(("./"+fileName+".tex").c_str());

  std::string bosonlabel, tablabel;
  if(fileName=="xSection_wplusEta"){
    bosonlabel="$W^{+}\\rightarrow\\mu^{+}\\nu$";
    tablabel="wplus";
  }else if(fileName=="xSection_wminEta"){
    bosonlabel="$W^{-}\\rightarrow\\mu^{-}\\nu$";
    tablabel="wminus";
  }

  xSecLatex << "\\begin{table}[tb!]" << std::endl;
  xSecLatex << "\\begin{center}" << std::endl;
  xSecLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) xSecLatex<<" c "; 
  xSecLatex<<" }\\hline \\hline " << std::endl;

  xSecLatex<< " & ";
  for(int k=0; k < nHeader; k++){
    xSecLatex<<std::setw(6)<<std::right;
    if(k< nHeader-1) xSecLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) xSecLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nValues; j++){
    xSecLatex<<column.at(j)<<" &";
    for(int l=0;l < nHeader; l++){
      if(std::isnan((totalValues.at(j)).at(l)) || std::isinf((totalValues.at(j)).at(l)))(totalValues.at(j)).at(l)=0.;
      xSecLatex<<std::setw(6)<<std::right;
      if(l< nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} & ";
      if(l==nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} \\\\"<<std::endl;
    }
  }
  
  xSecLatex << "\\hline \\hline"<< std::endl;
  xSecLatex << "\\end{tabular}" << std::endl;
  xSecLatex << "\\end{center}"  << std::endl;
  xSecLatex << "\\caption{This table contains the inclusive cross section and muon charge asymmetry results. The first column displays the measured quantity; the second and third columns are the data and predicted cross sections and charge asymmetry computed with Equations (\\ref{eq:master}) and (\\ref{eq:as}), respectively; the fourth column contains the unfolding factor; the fifth, sixth and seventh columns show the up, down and statistical uncertainties in percentage, respectively; last right four columns contain the data, reconstructed, truth and background number of events in that order.}" << std::endl;  
  xSecLatex << "\\label{tab:xsec"<<tablabel<<"}" << std::endl;  
  xSecLatex << "\\end{table}" << std::endl;
  xSecLatex.close();

  return;
}

void CrossSection::xSecTable(Config config, std::vector<std::vector<double>> totalValues, std::vector<std::string> column, std::vector<std::string> header, std::string fileName){

  if(header.size()!=totalValues.at(1).size()) return;

  int nHeader = (int)header.size(); int nValues = (int)totalValues.size();
  std::ofstream xSecLatex(("./"+fileName+".tex").c_str());

  std::string bosonlabel, tablabel;
  if(fileName=="xSection_wplusEta"){
    bosonlabel="$W^{+}\\rightarrow\\mu^{+}\\nu$";
    tablabel="wplus";
  }else if(fileName=="xSection_wminEta"){
    bosonlabel="$W^{-}\\rightarrow\\mu^{-}\\nu$";
    tablabel="wminus";
  }else{
    tablabel="asymmetry";
  }

  xSecLatex << "\\begin{table}[tb!]" << std::endl;
  xSecLatex << "\\begin{center}" << std::endl;
  xSecLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) xSecLatex<<" c "; 
  xSecLatex<<" }\\hline \\hline " << std::endl;

  xSecLatex<< " & ";
  for(int k=0; k < nHeader; k++){
    xSecLatex<<std::setw(6)<<std::right;
    if(k< nHeader-1) xSecLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) xSecLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nValues; j++){
    xSecLatex<<column.at(j)<<" &";
    for(int l=0;l < nHeader; l++){
      if(std::isnan((totalValues.at(j)).at(l)) || std::isinf((totalValues.at(j)).at(l)))(totalValues.at(j)).at(l)=0.;
      xSecLatex<<std::setw(6)<<std::right;
      if(l< nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} & ";
      if(l==nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} \\\\"<<std::endl;
    }
  }
  
  xSecLatex << "\\hline \\hline"<< std::endl;
  xSecLatex << "\\end{tabular}" << std::endl;
  xSecLatex << "\\end{center}"  << std::endl;
  xSecLatex << "\\caption{This table contains the cross section values for "<<bosonlabel<<". The first column displays the muon pseudorapidity bins; the second and third columns are the data measured and predicted cross sections computed with \\cref{eq:masterdif};  the fourth column contains the unfolding factor; the fifth, sixth and seventh columns show the up, down and statistical uncertainties in percentage, respectively; last right four columns contain the data, reconstructed, truth and background number of events in that order.}" << std::endl;  
  xSecLatex << "\\label{tab:xsec"<<tablabel<<"}" << std::endl;  
  xSecLatex << "\\end{table}" << std::endl;
  xSecLatex.close();

  return;
}


void CrossSection::asyTable(Config config, std::vector<std::vector<double>> totalValues, std::vector<std::string> column, std::vector<std::string> header, std::string fileName){

  if(header.size()!=totalValues.at(1).size()) return;

  int nHeader = (int)header.size(); int nValues = (int)totalValues.size();
  std::ofstream xSecLatex(("./"+fileName+".tex").c_str());

  std::string bosonlabel, tablabel;
  if(fileName=="xSection_wplusEta"){
    bosonlabel="$W^{+}\\rightarrow\\mu^{+}\\nu$";
    tablabel="wplus";
  }else if(fileName=="xSection_wminEta"){
    bosonlabel="$W^{-}\\rightarrow\\mu^{-}\\nu$";
    tablabel="wminus";
  }else{
    tablabel="asymmetry";
  }

  xSecLatex << "\\begin{table}[tb!]" << std::endl;
  xSecLatex << "\\begin{center}" << std::endl;
  xSecLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) xSecLatex<<" c "; 
  xSecLatex<<" }\\hline \\hline " << std::endl;

  xSecLatex<< " & ";
  for(int k=0; k < nHeader; k++){
    xSecLatex<<std::setw(6)<<std::right;
    if(k< nHeader-1) xSecLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) xSecLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nValues; j++){
    xSecLatex<<column.at(j)<<" &";
    for(int l=0;l < nHeader; l++){
      if(std::isnan((totalValues.at(j)).at(l)) || std::isinf((totalValues.at(j)).at(l)))(totalValues.at(j)).at(l)=0.;
      xSecLatex<<std::setw(6)<<std::right;
      if(l< nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} & ";
      if(l==nHeader-1) xSecLatex<<"\\num{"<<(totalValues.at(j)).at(l)<<"} \\\\"<<std::endl;
    }
  }
  
  xSecLatex << "\\hline \\hline"<< std::endl;
  xSecLatex << "\\end{tabular}" << std::endl;
  xSecLatex << "\\end{center}"  << std::endl;
  xSecLatex << "\\caption{This table contains the cross section values for "<<bosonlabel<<". The first column displays the muon pseudorapidity bins; the second and third columns are the data measured and predicted cross sections computed with \\cref{eq:masterdif};  the fourth column contains the unfolding factor; the fifth, sixth and seventh columns show the up, down and statistical uncertainties in percentage, respectively; last right four columns contain the data, reconstructed, truth and background number of events in that order.}" << std::endl;  
  xSecLatex << "\\label{tab:xsec"<<tablabel<<"}" << std::endl;  
  xSecLatex << "\\end{table}" << std::endl;
  xSecLatex.close();

  return;
}



void CrossSection::xSecTableSum(std::vector<std::vector<double>> totalValues,std::vector<std::vector<double>> stat,
				std::vector<std::vector<double>> sys_up, std::vector<std::vector<double>> sys_down, 
				std::vector<std::string> column, std::vector<std::string> header, std::string fileName){

  //if(header.size()!=totalValues.at(1).size()) return;
  int nHeader = (int)header.size(); int nValues = (int)totalValues.size();//nvalues=numberetabins
  std::ofstream xSecLatex(("./"+fileName+".tex").c_str());

  xSecLatex << "\\begin{table}[tb!]" << std::endl;
  xSecLatex << "\\begin{center}" << std::endl;
  xSecLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) xSecLatex<<" c "; 
  xSecLatex<<" }\\hline \\hline " << std::endl;

  xSecLatex<< " & ";
  for(int k=0; k < nHeader; k++){
    xSecLatex<<std::setw(6)<<std::right;
    if(k< nHeader-1) xSecLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) xSecLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < 11; j++){
    double xsec_wplus  = (totalValues.at(0)).at(j)*1000.;
    double xsec_wmin   = (totalValues.at(1)).at(j)*1000;
    double asymmetry   = (totalValues.at(2)).at(j);

    double stat_wplus  = (stat.at(0)).at(j); stat_wplus   = stat_wplus  ;
    double stat_wmin   = (stat.at(1)).at(j); stat_wmin    = stat_wmin   ;  
    double stat_asym   = (stat.at(2)).at(j); stat_asym    = 2   *asymmetry *0.01;  

    double up_wplus    = (sys_up.at(0)).at(j);  up_wplus  = up_wplus   *xsec_wplus ;
    double up_wmin     = (sys_up.at(1)).at(j);	up_wmin   = up_wmin    *xsec_wmin;  
    double up_asym     = (sys_up.at(2)).at(j);	up_asym   = up_asym    *asymmetry  *0.01;  
						                         
    double dn_wplus    = (sys_down.at(0)).at(j); dn_wplus = dn_wplus   *xsec_wplus ;
    double dn_wmin     = (sys_down.at(1)).at(j); dn_wmin  = dn_wmin    *xsec_wmin; 
    double dn_asym     = (sys_down.at(2)).at(j); dn_asym  = dn_asym    *asymmetry  *0.01; 

    double lu_wplus   =  0.024  *xsec_wplus ;
    double lu_wmin    =  0.024  *xsec_wmin;  
    double lu_asym    =  0.024  *asymmetry  ;
    
    xSecLatex<<column.at(j)<<" &";
    xSecLatex<<"\\num{"<<xsec_wplus<<"}$\\pm\\num{"<<stat_wplus<<"} ^{+\\num{"<<up_wplus<<"}} _{-\\num{"<<dn_wplus<<"}} \\pm$\\num{"<<lu_wplus<<"} & ";
    xSecLatex<<"\\num{"<<xsec_wmin <<"}$\\pm\\num{"<<stat_wmin <<"} ^{+\\num{"<<up_wmin <<"}} _{-\\num{"<<dn_wmin <<"}} \\pm$\\num{"<<lu_wmin <<"} & ";
    xSecLatex<<"\\num[round-mode=places,round-precision=3]{"<<asymmetry <<"}$\\pm\\num[round-mode=places,round-precision=3]{"<<stat_asym <<"} ^{+\\num[round-mode=places,round-precision=3]{"<<up_asym <<"}} _{-\\num[round-mode=places,round-precision=3]{"<<dn_asym <<"}} \\pm$\\num[round-mode=places,round-precision=3]{"<<lu_asym <<"}   ";
    xSecLatex<<"\\\\"<<std::endl;
  }
  
  xSecLatex << "\\hline \\hline"<< std::endl;
  xSecLatex << "\\end{tabular}" << std::endl;
  xSecLatex << "\\end{center}"  << std::endl;
  xSecLatex << "\\caption{This table displays a summary of the charged current Drell-Yan cross sections and the muon charge asymmetry measurements. The first column is the muon pseudorapidity binning; the second and the third columns are the $W^{-}\\rightarrow\\mu^{-}\\nu$ and $W^{+}\\rightarrow\\mu^{+}\\nu$ cross sections, respectively, calculated utilising \\cref{eq:masterdif}; the fourth column is the muon charge asymmetry computed with \\cref{eq:asdif}. The numbers are presented in the following order: the value of the process, the data statistical uncertainty, the systematic uncertainties up and down, and the luminosity uncertainty.}" << std::endl;  
  xSecLatex << "\\label{tab:summary}" << std::endl;  
  xSecLatex << "\\end{table}" << std::endl;
  xSecLatex.close();

  return;
}


std::vector<std::string> CrossSection::Names(std::string option){

  std::vector<std::string> names;
  if(option=="systematics"){

    names.push_back("MUON\\_ID");//0
    names.push_back("MUON\\_MS");
    names.push_back("MUON\\_SAGITTA\\_RESBIAS");
    names.push_back("MUON\\_SAGITTA\\_RHO");
    names.push_back("MUON\\_SCALE");
    
    names.push_back("MET\\_SoftTrk\\_ResoPara");//5
    names.push_back("MET\\_SoftTrk\\_ResoPerp");
    names.push_back("MET\\_SoftTrk\\_Scale");
    names.push_back("PileUp");
    
    names.push_back("IDStat");//9
    names.push_back("IDSys");
    names.push_back("IsoStat");
    names.push_back("IsoSys");	      
    
    names.push_back("TrigStat");//13
    names.push_back("TrigSys");
    names.push_back("TTVAStat");
    names.push_back("TTVASys");
    
    names.push_back("ALPHAS");//17
    names.push_back("BEAM\\_ENERGY");
    names.push_back("PDF\\_EW");
    names.push_back("PDF");
    names.push_back("PI");
    
    names.push_back("SCALE\\_W");//22
    names.push_back("SCALE\\_Z");
    names.push_back("CHOICE\\_HERAPDF20");
    names.push_back("CHOICE\\_NNPDF30");
    
    names.push_back("PDF\\_EV1");//26
    names.push_back("PDF\\_EV2");
    names.push_back("PDF\\_EV3");
    names.push_back("PDF\\_EV4");
    names.push_back("PDF\\_EV5");
    names.push_back("PDF\\_EV6");
    names.push_back("PDF\\_EV7");
    
    names.push_back("JET\\_EtaIntercalibration\\_NonClosure");//33
    names.push_back("JET\\_JER\\_DataVsMC");
    names.push_back("JET\\_GroupedNP\\_1");
    names.push_back("JET\\_GroupedNP\\_2");
    names.push_back("JET\\_GroupedNP\\_3");
    
    names.push_back("JET\\_JER\\_EffectiveNP\\_1");//38
    names.push_back("JET\\_JER\\_EffectiveNP\\_2");
    names.push_back("JET\\_JER\\_EffectiveNP\\_3");
    names.push_back("JET\\_JER\\_EffectiveNP\\_4");
    
    names.push_back("JET\\_JER\\_EffectiveNP\\_5");//42
    names.push_back("JET\\_JER\\_EffectiveNP\\_6");
    names.push_back("JET\\_JER\\_EffectiveNP\\_7restTerm");

    names.push_back("Monte Carlo statistics");    
  }

  if(option=="xSection"){
    names.push_back("$\\sigma_{Data}$ [nb]");
    names.push_back("$\\sigma_{MC}$ [nb]");
    names.push_back("$C_{W}$");
    names.push_back("$\\delta_{sys}^{up}$");
    names.push_back("$\\delta_{sys}^{down}$");
    names.push_back("$\\delta_{stat}$");
    names.push_back("$N^{Data}$");
    names.push_back("$N^{reco}$");
    names.push_back("$N^{truth}$");
    names.push_back("$N^{BG}$");
  }

  if(option=="summary"){
    names.push_back("$\\sigma_{W^{-}}$[pb]");
    names.push_back("$\\sigma_{W^{+}}$[pb]");
    names.push_back("$As_{\\mu}$");

    // names.push_back("$\\sigma_{W^{-}}$ (value$\\pm$stat$\\pm$sys$\\pm$lumi) [pb]");
    // names.push_back("$\\sigma_{W^{+}}$ (value$\\pm$stat$\\pm$sys$\\pm$lumi)[pb]");
    // names.push_back("$As_{\\mu}$ (value$\\pm$stat$\\pm$sys$\\pm$lumi)");
  }

  if(option=="xSection_inc"){
    names.push_back("$Data$");
    names.push_back("$MC$");
    names.push_back("$C_{W}$");
    names.push_back("$\\delta_{sys}^{up}$");
    names.push_back("$\\delta_{sys}^{down}$");
    names.push_back("$\\delta_{stat}$");
    names.push_back("$N^{Data}$");
    names.push_back("$N^{reco}$");
    names.push_back("$N^{truth}$");
    names.push_back("$N^{BG}$");
  }


  if(option=="Asymmetry"){
    names.push_back("$As_{\\mu,Data}$");
    names.push_back("$As_{\\mu,MC}$");
    names.push_back("$\\delta_{sys}^{up}$");
    names.push_back("$\\delta_{sys}^{down}$");
    names.push_back("$\\delta_{stat}$");
  }

  return names;
}


std::vector<std::string> CrossSection::LatexBinName(){

  std::vector<std::string> name; name.clear();
  //name.push_back(" Inclusive $\\vert\\eta^{\\mu}\\vert$ ");     	     
  name.push_back(" $ 0.0\\leq\\vert\\eta^{\\mu}\\vert<0.21$");
  name.push_back(" $0.21\\leq\\vert\\eta^{\\mu}\\vert<0.42$");
  name.push_back(" $0.42\\leq\\vert\\eta^{\\mu}\\vert<0.63$");
  name.push_back(" $0.63\\leq\\vert\\eta^{\\mu}\\vert<0.84$");
  name.push_back(" $0.84\\leq\\vert\\eta^{\\mu}\\vert<1.05$");
  name.push_back(" $1.05\\leq\\vert\\eta^{\\mu}\\vert<1.37$");
  name.push_back(" $1.37\\leq\\vert\\eta^{\\mu}\\vert<1.52$");
  name.push_back(" $1.52\\leq\\vert\\eta^{\\mu}\\vert<1.74$");
  name.push_back(" $1.74\\leq\\vert\\eta^{\\mu}\\vert<1.95$");
  name.push_back(" $1.95\\leq\\vert\\eta^{\\mu}\\vert<2.18$");
  name.push_back(" $2.18\\leq\\vert\\eta^{\\mu}\\vert<2.40$");
 
 return name;
 
}


void CrossSection::xSecPlot(Config config, TH1D *xSecHist, TH1D *xSecHistMC, std::string boson, double ylow, double yhigh){

  std::string xlabel,ylabel; double xlow=0., xhigh=0.;
  TString leglabel;

  xlabel="#||{#eta^{#mu}}";
  if(boson=="zmumu")    leglabel="Z#rightarrow#mu#mu";
  if(boson=="wplus")    leglabel="W^{+}#rightarrow#mu^{+}#nu";
  if(boson=="wminus")   leglabel="W^{-}#rightarrow#mu^{-}#nu";
  if(boson=="asymmetry") leglabel="";

  ylabel="d#sigma/d#||{#eta^{#mu}} [pb]";  
  if(boson=="asymmetry")ylabel = "Asymmetry";

  plotAxisLine(xSecHist,kBlack,kGreen,20,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(xSecHistMC,kBlue,kGreen,28,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);

  TString lumi="44.3 fb^{-1}";
  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.175,0.85,0.375,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  TH1D *hAux=(TH1D*)xSecHist->Clone();
  xSecHist->SetFillColor(kGreen);
  xSecHist->SetLineWidth(0);
  xSecHist->SetFillStyle(1001);

  leg->AddEntry(hAux, "#bf{#scale[0.85]{"" Data}}","PL");  
  leg->AddEntry(xSecHistMC, "#bf{#scale[0.85]{"" Powheg+Pythia}}","PL");  
  leg->AddEntry(xSecHist, "#bf{#scale[0.85]{ Systematic Uncertainty}}","f");  

  TPaveText *box;
  box = new TPaveText(0.2,0.865,0.375,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi);

  TPaveText *box2;
  box2 = new TPaveText(0.445,0.895,0.495,0.975,"NDC");
  box2->SetBorderSize(0); box2->SetTextSize(0.035); box2->SetFillColor(0);
  box2->AddText(leglabel);
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.2,1,1);
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.015);
  pad->SetBottomMargin(0.1125);

  xSecHist->SetTitle("");
  xSecHist->GetXaxis()->SetTitleOffset(1.45);

  xSecHist->DrawCopy("e2");
  hAux->DrawCopy("p samese hist");
  xSecHistMC->DrawCopy("samese p");

  pad->Update();
  pad->Modified();
  box->Draw();
  box2->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = xSecHist->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[xSecHist->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< xSecHist->GetXaxis()->GetNbins()+1; i++)  xbins[i]=xSecHist->GetXaxis()->GetBinLowEdge(i+1);
  xbins[xSecHist->GetXaxis()->GetNbins()]=xSecHist->GetXaxis()->GetBinUpEdge(xSecHist->GetXaxis()->GetNbins());

  TH1D* ratiodataMC   = new TH1D("ratiodataMC","",nx,xbins);
  
  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){
    
    double data=xSecHist->GetBinContent(b);
    if(data==0) continue;
    double mc=xSecHistMC->GetBinContent(b);
    if(mc==0) continue;
    
    double dataratio=data/mc;
    double err=xSecHist->GetBinError(b)/mc;
    ratiodataMC->SetBinContent(b,1./dataratio);
    ratiodataMC->SetBinError(b,err);
  }

  ratioSettings(ratiodataMC, 0.74, 1.24,"","Theory/Data",0.1125,0.1125,0.35,0.15,28,kBlue,0.75,5);
  ratiodataMC->DrawCopy("p");

  TLine *line;
  line = new TLine(0.,1.,2.4,1.);
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();

  c->cd();
  c->Print((config.OutputFileDir + "xSection/"+ config.DataYears +"/xSec_"+boson+".pdf").c_str());  
  delete c;

  return;
}

void CrossSection::plotAxisLine(TH1D* hist, int lineColor, int markerColor,
			      int markerStyle, double markerSize,
			      TString title, TString xlabel, TString ylabel, bool xRange,
			      double xlow, double xhigh, bool yRange, double ylow, double yhigh)
{
  hist->SetLineColor(lineColor);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerColor(lineColor);
  hist->SetMarkerSize(markerSize);
  hist->SetTitle(title);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);

  if(xRange==true)
    hist->GetXaxis()->SetRangeUser(xlow,xhigh);
  if(yRange==true)
    hist->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}

void CrossSection::ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
				 double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
				 int markerstyle, int color, double markersize, int divisions)
{
  hist->SetMinimum(min);
  hist->SetMaximum(max);
  hist->GetXaxis()->SetTitle(xlabel);

  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->SetLabelSize(xlabelsize);
  hist->GetYaxis()->SetLabelSize(ylabelsize);
  hist->GetYaxis()->SetTitleOffset(ytitleof);
  hist->GetYaxis()->SetTitleSize(ytitlesize);
  hist->SetLineColor(color);

  hist->Sumw2();
  hist->SetMarkerStyle(markerstyle);
  hist->SetMarkerSize(markersize);
  hist->SetMarkerColor(color);

  hist->SetStats(0);
  hist->GetYaxis()->SetNdivisions(divisions);

  return;
}

void CrossSection::setstyle(){

  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(000000);
  gStyle->SetOptTitle(000000);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelFont(22,"x");
  gStyle->SetLabelFont(22,"y");

  gStyle->SetTitleSize(0.035,"x");
  gStyle->SetTitleSize(0.035,"y");
  gStyle->SetTitleOffset(1.35,"x");
  gStyle->SetTitleOffset(1.35,"y");

  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.18);

  int font = 42;
  gStyle->SetTextFont(font);

  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");

  return;
}

void CrossSection::finalize(){

  fout->cd();
  // hEtaBinXsec->Write("hEtaBinXsec");
  // hEtaBinXsecP->Write("hEtaBinXsecN");
  // hEtaBinAsymm->Write("hEtaBinAsymm");
  fout->Close();

  return;
}

#endif
