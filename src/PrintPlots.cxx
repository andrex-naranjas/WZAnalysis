//Plotting Class includes
#ifndef PRINTPLOTS_CXX
#define PRINTPLOTS_CXX

#include "PrintPlots.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>


PrintPlots::PrintPlots()
{
}

PrintPlots::~PrintPlots(){}

void PrintPlots::initialize(Config config)
{

  //ouput directory according year
  std::string dirYear="";
  if(config.DataYears=="2015+2016") dirYear="2015p2016/";
  if(config.DataYears=="2017") dirYear="2017/";
  if(config.DataYears=="2018") dirYear="2018/";

  //zmumu
  if(config.WZSelection=="zmumu"){
    fdata  = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_data15p16_z_set_insitu.root");
    fmcor  = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z_nonewPU.root");
    fmc    = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z.root");
    fmccal = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z_set_insitu.root");
    fmcfull= new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z_set_insitu_resolresponse.root");

    //fmcfull= new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z_set_insitu.root");

    wzdata  = new TFile((config.OutputFileDir+"Files/" +dirYear+ "data15p16_z.root").c_str());
    wzmcor  = new TFile((config.OutputFileDir+"Files/" +dirYear+ "zmumu_z_nominal_nonewPU.root").c_str());
    wzmc    = new TFile((config.OutputFileDir+"Files/" +dirYear+ "zmumu_z_nominal.root").c_str());
    wzmccal = new TFile((config.OutputFileDir+"Files/" +dirYear+ "zmumu_z_nominal_set_insitu.root").c_str());
    wzmcfull= new TFile((config.OutputFileDir+"Files/" +dirYear+ "zmumu_z_nominal_set_insitu_resolresponse.root").c_str());

 }
  //wplus
  if(config.WZSelection=="wplus"){
    fdata  = new TFile("/data/morales/atlas/RecoilFiles2018/RecoilHistos/RecoilCalibHistos_data15p16.root");
    fmc    = new TFile("/data/morales/atlas/RecoilFiles2018/RecoilCalibHistos_wplus.root");
    fmccal = new TFile("/data/morales/atlas/RecoilFiles2018/RecoilHistos/RecoilCalibHistos_wplus_calib.root");
  }
  //wminus
  if(config.WZSelection=="wminus"){
   
    fdata  = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_data15p16_wminus_set_insitu_resolresponse.root");
    fmcor  = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_wminmunu_wminus_nonewPU.root");
    fmc    = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_wminmunu_wminus.root");
    fmccal = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_wminmunu_wminus_set_insitu.root");
    fmcfull= new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_wminmunu_wminus_set_insitu_resolresponse.root");
    //fmcfull= new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_z_set_insitu.root");

    wzdata  = new TFile((config.OutputFileDir+"Files/" +dirYear+ "data15p16_wminus.root").c_str());
    wzmcor  = new TFile((config.OutputFileDir+"Files/" +dirYear+ "wminmunu_wminus_nominal_nonewPU.root").c_str());
    wzmc    = new TFile((config.OutputFileDir+"Files/" +dirYear+ "wminmunu_wminus_nominal.root").c_str());
    wzmccal = new TFile((config.OutputFileDir+"Files/" +dirYear+ "wminmunu_wminus_nominal_set_insitu.root").c_str());
    wzmcfull= new TFile((config.OutputFileDir+"Files/" +dirYear+ "wminmunu_wminus_nominal_set_insitu_resolresponse.root").c_str());

    //recoil backgrounds
    if(config.RecoilBG=="True"){
      wtaunu  = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_wmintaunu_wminus_set_insitu_resolresponse.root");
      zmumu   = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_zmumu_wminus_set_insitu_resolresponse.root");
      ztautau = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_ztautau_wminus_set_insitu_resolresponse.root");
      ttbar   = new TFile("~/ATLAS/ControlPlots2018/Plots/RecoilHistos/RecoilCalibHistos_ttbar_wminus_set_insitu_resolresponse.root");
      
      pwtaunu  = new TFile((config.OutputFileDir+"Files/" +dirYear+ "wmintaunu_wminus_nominal_set_insitu_resolresponse.root").c_str());
      pzmumu   = new TFile((config.OutputFileDir+"Files/" +dirYear+     "zmumu_wminus_nominal_set_insitu_resolresponse.root").c_str());
      pztautau = new TFile((config.OutputFileDir+"Files/" +dirYear+   "ztautau_wminus_nominal_set_insitu_resolresponse.root").c_str());
      pttbar   = new TFile((config.OutputFileDir+"Files/" +dirYear+     "ttbar_wminus_nominal_set_insitu_resolresponse.root").c_str());
    }

  }                
  return;
}

void PrintPlots::execute(Config config){

  std::vector<std::string> kine, nplots;
  kine.push_back("mu");
  kine.push_back("Sumet"); kine.push_back("ptvtrue"); 
  kine.push_back("upar_rec"); kine.push_back("uperp_rec"); 
  kine.push_back("upar_ptv"); kine.push_back("upar_ptv_n");
  kine.push_back("ut"); kine.push_back("utphi"); 
  kine.push_back("met");
  kine.push_back("metphi"); 
  kine.push_back("ux");
  kine.push_back("uy");
  nplots.push_back("dataMC"); nplots.push_back("MC");
  
  for(int j=0; j<kine.size();j++) 
    for(int k=0; k<nplots.size();k++)
      ClosurePlot(config,kine[j],nplots[k],"calib");

  for(int j=0; j<kine.size();j++)
    if(config.RecoilBG=="True" && config.WZSelection!="zmumu") WStackPlot(config,kine[j],"calib");
    
  kine.clear();
  kine.push_back("met");
  kine.push_back("rmass");

  if(config.WZSelection=="wminus" || config.WZSelection=="wplus"){ 
    kine.push_back("met_recoil");
    kine.push_back("rmass_recoil");}
  kine.push_back("pt_boson_reco");//uT
  kine.push_back("pt_boson_reco_p");//uT
  
  for(int j=0; j<kine.size();j++)
    for(int k=0; k<nplots.size();k++)
      ClosurePlot(config,kine[j],nplots[k],"control");

  kine.push_back("SumET");
  for(int j=0; j<kine.size();j++)
    if(config.RecoilBG=="True" && config.WZSelection!="zmumu") WStackPlot(config,kine[j],"control");
  
  std::vector<std::string> profile;
  profile.push_back("ux_phi"); profile.push_back("uy_phi");
  profile.push_back("ux_set"); profile.push_back("uy_set");
  profile.push_back("ux_set"); profile.push_back("uy_set");
  profile.push_back("bias_set"); profile.push_back("uperp_set");
  profile.push_back("bias_pt"); profile.push_back("uperp_pt");
  for(int j=0; j<profile.size();j++) 
    for(int k=0; k<nplots.size();k++)
      ProfilePlot(config,profile[j],nplots[k]);

  profile.clear();
  profile.push_back("bias_set"); profile.push_back("uperp_set");
  profile.push_back("bias_pt"); profile.push_back("uperp_pt");
  for(int j=0; j<profile.size();j++) 
    for(int k=0; k<nplots.size();k++)
      WidthPlot(config,profile[j],nplots[k]);
  
  std::vector<std::string> rr;
  rr.push_back("resp_set_mean"); rr.push_back("resol_set_mean");
  rr.push_back("resp_set_width"); rr.push_back("resol_set_width");
  rr.push_back("resp_pt_mean"); rr.push_back("resol_pt_mean");
  rr.push_back("resp_pt_width"); rr.push_back("resol_pt_width");
  for(int j=0; j<rr.size();j++){
      //ResolResponse(config,rr[j],"dataMC");
}
      //for(int k=0; k<nplots.size();k++)
    			     
  return;
}

void PrintPlots::ClosurePlot(Config config, std::string kine, std::string nplots, std::string option)
{
  TH1D *hPlotData, *hPlotMCor, *hPlotMC, *hPlotMCmu, *hPlotMCfull;
 
  if(option=="calib"){
    hPlotData   = (TH1D*)fdata   ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCor   = (TH1D*)fmcor   ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMC     = (TH1D*)fmc     ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCmu   = (TH1D*)fmccal  ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCfull = (TH1D*)fmcfull ->GetObjectUnchecked(("h_"+kine).c_str());
  }else if(option=="control"){
    hPlotData   = (TH1D*)wzdata  ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCor   = (TH1D*)wzmcor  ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMC     = (TH1D*)wzmc    ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCmu   = (TH1D*)wzmccal ->GetObjectUnchecked(("h_"+kine).c_str());
    hPlotMCfull = (TH1D*)wzmcfull->GetObjectUnchecked(("h_"+kine).c_str());
  }

  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0.;
  TString leglabel;

  if(kine=="mu"){xlabel="#LT#mu#GT"; xlow=0; xhigh=50;}
  if(kine=="Sumet"){xlabel="#sum E_{T}^{PFO} [GeV]"; xlow=0.; xhigh=2000;}
  if(kine=="ptvtrue"){xlabel="p_{T}^{V} [GeV]"; xlow=0.; xhigh=200;}
  if(kine=="upar_rec"){xlabel="u_{par}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="upar_ptv"){xlabel="u_{par}^{V}+p_{T}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="upar_ptv_n"){xlabel="u_{par}^{V}-p_{T}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="uperp_rec"){xlabel="u_{perp}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="ut"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200.;}
  if(kine=="ux"){xlabel="u_{x} [GeV]"; xlow=-100.; xhigh=100.;}
  if(kine=="uy"){xlabel="u_{y} [GeV]"; xlow=-100.; xhigh=100.;}
  if(kine=="utphi"){xlabel="#phi(u_{T}) [rad]"; xlow=-5.; xhigh=5.;}
  if(kine=="met"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200.;}
  if(kine=="metphi"){xlabel="#phi(E_{T}^{miss})"; xlow=-5.; xhigh=5.;}

  if(kine=="met" && option =="control"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="rmass"){xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=150000.;}
  if(kine=="rmass" && config.WZSelection=="zmumu"){xlabel="m_{ll} [GeV]"; xlow=66000.; xhigh=116000.;}
  if(kine=="met_recoil" && option =="control"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="rmass_recoil"){xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=150000.;}
  if(kine=="pt_boson_reco"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="pt_boson_reco_p"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200*1000.;}

  if(config.WZSelection=="zmumu") leglabel="MC Z#rightarrow#mu#mu";
  if(config.WZSelection=="wminus" || config.WZSelection=="wminus") leglabel="MC W#rightarrow#mu#nu";

  if(nplots=="dataMC") ratiolabel="Data/MC"; 
  if(nplots=="MC") ratiolabel="Calib/NoCalib";

  if(option=="calib"){ hPlotData->Rebin(5); hPlotMCor->Rebin(5); hPlotMC->Rebin(5); hPlotMCmu->Rebin(5); hPlotMCfull->Rebin(5);}

  if(config.NormPlots=="True" /*&& kine!="rmass"*/){
    hPlotData->Scale(1./hPlotData->Integral());
    hPlotMCor->Scale(1./hPlotMCor->Integral());
    hPlotMC->Scale(1./hPlotMC->Integral());
    hPlotMCmu->Scale(1./hPlotMCmu->Integral());
    hPlotMCfull->Scale(1./hPlotMCfull->Integral());
    ylabel="A.U.";
  }else
    ylabel="Events";
    
  plotAxisLine(hPlotData,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,4.);
  plotAxisLine(hPlotMCor,kGray+2,kGray+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(hPlotMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(hPlotMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(hPlotMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2015+2016") lumi="36.1";
  if(config.DataYears=="2018")      lumi="59.9";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.6,0.85,0.8,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  if(nplots=="dataMC")
    leg->AddEntry(hPlotData,"#bf{#scale[0.85]{Data "+naming+"}}","PL");
  leg->AddEntry(hPlotMCor,"#bf{#scale[0.85]{"+leglabel+", Original}}","PL");  
  leg->AddEntry(hPlotMC,"#bf{#scale[0.85]{"+leglabel+", <#mu>}}","PL");
  leg->AddEntry(hPlotMCmu,"#bf{#scale[0.85]{"+leglabel+",<#mu>+SET+Insitu}}","PL");  
  leg->AddEntry(hPlotMCfull,"#bf{#scale[0.85]{"+leglabel+", Full}}","PL");  

  TPaveText *box;
  box = new TPaveText(0.625,0.865,0.8,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("13 TeV, "+lumi+" fb^{-1}");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);

  hPlotMC->SetTitle("");
  hPlotMC->GetXaxis()->SetTitleOffset(1.45);
  hPlotMC->DrawCopy("hist p");  
  if(nplots=="dataMC") hPlotData->DrawCopy("p sames hist");
  hPlotMCor->DrawCopy("hist p sames");
  hPlotMCfull->DrawCopy("hist p sames");
  hPlotMCmu->DrawCopy("p sames hist");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = hPlotData->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[hPlotData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hPlotData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hPlotData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hPlotData->GetXaxis()->GetNbins()]=hPlotData->GetXaxis()->GetBinUpEdge(hPlotData->GetXaxis()->GetNbins());

  TH1D* ratiodataMCor   = new TH1D("ratiodataMCor","",nx,xbins);
  TH1D* ratiodataMC   = new TH1D("ratiodataMC","",nx,xbins);
  TH1D* ratiodataMCmu = new TH1D("ratiodataMCmu","",nx,xbins);
  TH1D* ratiodataMCfull = new TH1D("ratiodataMCfull","",nx,xbins);
  TH1D* ratioMCmu = new TH1D("ratioMCmu","",nx,xbins);

  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){

    double data=hPlotData->GetBinContent(b);
    //if(data==0) continue;
    double mcor=hPlotMCor->GetBinContent(b);
    if(mcor==0) continue;
    double mc=hPlotMC->GetBinContent(b);
    if(mc==0) continue;
    double mcmu=hPlotMCmu->GetBinContent(b);
    if(mcmu==0) continue;
    double mcfull=hPlotMCfull->GetBinContent(b);
    if(mcfull==0) continue;

    double dataratioor=data/mcor;
    double dataratio=data/mc;
    double dataratiocalib=data/mcmu;
    double dataratiocalibfull=data/mcfull;
    double mcratio=mcmu/mc;

    double err=hPlotData->GetBinError(b)/mc;
    double error=hPlotData->GetBinError(b)/mc;
    double errmu=hPlotData->GetBinError(b)/mcmu;
    double errmc=hPlotMCmu->GetBinError(b)/mc;
    double errfull=hPlotMCfull->GetBinError(b)/mc;
    
    ratiodataMCor->SetBinContent(b,dataratioor);
    ratiodataMCor->SetBinError(b,err);

    ratiodataMC->SetBinContent(b,dataratio);
    ratiodataMC->SetBinError(b,err);

    ratiodataMCmu->SetBinContent(b,dataratiocalib);
    ratiodataMCmu->SetBinError(b,errmu);

    ratiodataMCfull->SetBinContent(b,dataratiocalibfull);
    ratiodataMCfull->SetBinError(b,0);

    ratioMCmu->SetBinContent(b,mcratio);
    ratioMCmu->SetBinError(b,errmc);

  }

  plotAxisLine(ratiodataMCor,kCyan,kCyan,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  plotAxisLine(ratiodataMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratioMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  ratioSettings(ratiodataMCor  ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kGray+2,0.5,10); 
  ratioSettings(ratiodataMC    ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kCyan+3,0.5,10);
  ratioSettings(ratiodataMCmu  ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 
  ratioSettings(ratiodataMCfull,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kViolet,0.5,10); 
  ratioSettings(ratioMCmu      ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 

  if(nplots=="dataMC")ratiodataMC->DrawCopy("p");
  if(nplots=="dataMC")ratiodataMCfull->DrawCopy("p same");
  if(nplots=="dataMC")ratiodataMCor->DrawCopy("p same");
  if(nplots=="dataMC")ratiodataMCmu->DrawCopy("p same");
  if(nplots=="MC")ratioMCmu->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,1.,xhigh,1.);
  //line = new TLine(40000.,1.,150000.,1.);//remove
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  std::string boson;
  if(config.WZSelection=="zmumu") boson="z";
  if(config.WZSelection=="wminus") boson="wminus";
  if(config.WZSelection=="wplus") boson="wplus";

  c->Print(("./ClosurePlots/closure_"+boson+"_"+kine+"_"+nplots+"_"+config.DataYears+".pdf").c_str());
  
  delete hPlotData; delete hPlotMCor; delete hPlotMC; delete hPlotMCmu; delete hPlotMCfull;
  return;
}


void PrintPlots::ProfilePlot(Config config, std::string profile, std::string nplots)
{
  TProfile *p_PlotData, *p_PlotMC, *p_PlotMCor, *p_PlotMCmu, *p_PlotMCfull;

  p_PlotData   = (TProfile*)fdata   ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCor   = (TProfile*)fmcor   ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMC     = (TProfile*)fmc     ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCmu   = (TProfile*)fmccal  ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCfull = (TProfile*)fmcfull ->GetObjectUnchecked(("p_"+profile).c_str());

  if(profile=="ux_set" || profile=="uy_set" || profile=="uperp_set" || profile=="bias_set")
    {p_PlotData->Rebin(40);  p_PlotMCor->Rebin(40); p_PlotMC->Rebin(40); p_PlotMCmu->Rebin(40); p_PlotMCfull->Rebin(40);}

  if(profile=="uperp_pt" || profile=="bias_pt")
    {p_PlotData->Rebin(10);  p_PlotMCor->Rebin(10); p_PlotMC->Rebin(10); p_PlotMCmu->Rebin(10); p_PlotMCfull->Rebin(10);}

  // if(profile=="ux_phi" || profile=="uy_phi")
  //   {p_PlotData->Rebin(2);  p_PlotMCor->Rebin(2); p_PlotMC->Rebin(2); p_PlotMCmu->Rebin(2); p_PlotMCmu->Rebin(2);}
  
  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0., ylow=0., yhigh=0.;
  bool yRange=false;
  TString leglabel;
  if(profile=="ux_phi"){xlabel="#phi^{lepton} [rad]"; ylabel="#LT u_{x}#GT [GeV]"; xlow=-5.; xhigh=5;}
  if(profile=="uy_phi"){xlabel="#phi^{lepton} [rad]"; ylabel="#LT u_{y}#GT [GeV]"; xlow=-5.; xhigh=5;}
  if(profile=="ux_set"){xlabel="#sum E_{T}^{PFO} [GeV]";  ylabel="#LT u_{x}#GT [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=20;}
  if(profile=="uy_set"){xlabel="#sum E_{T}^{PFO} [GeV]"; ylabel="#LT u_{y}#GT [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=10;}

  if(profile=="uperp_set"){xlabel="#sum E_{T}^{PFO} [GeV]";  ylabel="#LT u_{T}#GT [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=20;}
  if(profile=="bias_set"){xlabel="#sum E_{T}^{PFO} [GeV]"; ylabel="#LT u_{||} + p_{T}^{V} #GT [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=40;}

  if(profile=="uperp_pt"){xlabel="p_{T}^{V} [GeV]";  ylabel="#LT u_{T}#GT [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-2.5; yhigh=5;}
  if(profile=="bias_pt"){xlabel="p_{T}^{V} [GeV]"; ylabel="#LT u_{||} + p_{T} #GT [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-7.5; yhigh=40;}

  if(config.WZSelection=="zmumu") leglabel="MC Z#rightarrow#mu#mu";
  if(config.WZSelection=="wminus" || config.WZSelection=="wminus") leglabel="MC W#rightarrow#mu#nu";
  if(nplots=="dataMC") ratiolabel="Data-MC [GeV]";
  if(nplots=="MC") ratiolabel="Calib-NoCalib [GeV]";

  plotAxisLine(p_PlotData,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMCor,kGray+2,kGray+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  yRange=false;

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2015+2016") lumi="36.1";
  if(config.DataYears=="2018")      lumi="59.9";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.15,0.85,0.45,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  if(nplots=="dataMC")
    leg->AddEntry(p_PlotData,"#bf{#scale[0.85]{Data "+naming+"}}","PL");

  leg->AddEntry(p_PlotMCor,"#bf{#scale[0.85]{"+leglabel+", Original}}","PL");
  leg->AddEntry(p_PlotMC,"#bf{#scale[0.85]{"+leglabel+", <#mu>}}","PL");
  leg->AddEntry(p_PlotMCmu,"#bf{#scale[0.85]{"+leglabel+", <#mu>+SET+Insitu}}","PL");
  leg->AddEntry(p_PlotMCfull,"#bf{#scale[0.85]{"+leglabel+", Full}}","PL");
    
  TPaveText *box;
  box = new TPaveText(0.2,0.865,0.375,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("13 TeV, "+lumi+" fb^{-1}"); 
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.02);

  p_PlotMC->SetTitle("");
  p_PlotMC->GetXaxis()->SetTitleOffset(1.45);
  p_PlotMC->DrawCopy("hist p");  
  if(nplots=="dataMC") p_PlotData->DrawCopy("hist p sames");
  p_PlotMCor->DrawCopy("hist p sames");
  p_PlotMCmu->DrawCopy("hist p sames");
  p_PlotMCfull->DrawCopy("hist p sames");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = p_PlotData->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[p_PlotData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< p_PlotData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=p_PlotData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[p_PlotData->GetXaxis()->GetNbins()]=p_PlotData->GetXaxis()->GetBinUpEdge(p_PlotData->GetXaxis()->GetNbins());

  TH1D* diffdataMC     = new TH1D("diffdataMC","",nx,xbins);
  TH1D* diffdataMCor     = new TH1D("diffdataMCor","",nx,xbins);
  TH1D* diffdataMCmu   = new TH1D("diffdataMCmu","",nx,xbins);
  TH1D* diffdataMCfull = new TH1D("diffdataMCfull","",nx,xbins);
  TH1D* diffMCmu = new TH1D("diffMCmu","",nx,xbins);

  for(int b=1; b< (2+(diffdataMC->GetNbinsX())); b++){
    double data=p_PlotData->GetBinContent(b);
    double mcor=p_PlotMCor->GetBinContent(b);
    double mc=p_PlotMC->GetBinContent(b);
    double mcmu=p_PlotMCmu->GetBinContent(b);
    double mcfull=p_PlotMCfull->GetBinContent(b);
    double errd=p_PlotData->GetBinError(b);
    
    diffdataMCor->SetBinContent(b,data-mcor);
    diffdataMCor->SetBinError(b,errd);

    diffdataMC->SetBinContent(b,data-mc);
    diffdataMC->SetBinError(b,errd);

    diffdataMCmu->SetBinContent(b,data-mcmu);
    diffdataMCmu->SetBinError(b,errd);

    diffdataMCfull->SetBinContent(b,data-mcfull);
    diffdataMCfull->SetBinError(b,errd);

    diffMCmu->SetBinContent(b,mcmu-mc);
    diffMCmu->SetBinError(b,errd);
  }

  plotAxisLine(diffdataMCor,kGray+2,kGray+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(diffdataMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(diffdataMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  plotAxisLine(diffdataMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(diffMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  ratioSettings(diffdataMCor,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kGray+2,0.5,10); 
  ratioSettings(diffdataMC,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kCyan+3,0.5,10); 
  ratioSettings(diffdataMCmu,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 
  ratioSettings(diffdataMCfull,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kViolet,0.5,10); 
  ratioSettings(diffMCmu,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 

  if(nplots=="dataMC") diffdataMC->DrawCopy("hist p");
  if(nplots=="dataMC") diffdataMCor->DrawCopy("hist p same");
  if(nplots=="dataMC") diffdataMCmu->DrawCopy("hist p same");
  if(nplots=="dataMC") diffdataMCfull->DrawCopy("hist p same");
  if(nplots=="MC")     diffMCmu->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,0.,xhigh,0.);
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  std::string boson;
  if(config.WZSelection=="zmumu") boson="z";
  if(config.WZSelection=="wminus") boson="wminus";
  if(config.WZSelection=="wplus") boson="wplus";

  c->Print(("./ClosurePlots/closure_"+boson+"_"+profile+"_"+nplots+"_"+config.DataYears+"_profile.pdf").c_str());

  delete diffdataMC; delete diffdataMCmu;  
  delete p_PlotData; delete p_PlotMCor; delete p_PlotMC; delete p_PlotMCmu; delete p_PlotMCfull;
  return;
}


void PrintPlots::WidthPlot(Config config, std::string profile, std::string nplots)
{
  TProfile *p_PlotData, *p_PlotMC, *p_PlotMCor, *p_PlotMCmu, *p_PlotMCfull;

  p_PlotData   = (TProfile*)fdata   ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCor   = (TProfile*)fmcor   ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMC     = (TProfile*)fmc     ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCmu   = (TProfile*)fmccal  ->GetObjectUnchecked(("p_"+profile).c_str());
  p_PlotMCfull = (TProfile*)fmcfull ->GetObjectUnchecked(("p_"+profile).c_str());

  if(profile=="uperp_set" || profile=="bias_set")
    {p_PlotData->Rebin(40);  p_PlotMCor->Rebin(40); p_PlotMC->Rebin(40); p_PlotMCmu->Rebin(40); p_PlotMCfull->Rebin(40);}

  if(profile=="uperp_pt" || profile=="bias_pt")
    {p_PlotData->Rebin(10);  p_PlotMCor->Rebin(10); p_PlotMC->Rebin(10); p_PlotMCmu->Rebin(10); p_PlotMCfull->Rebin(10);}

  unsigned int nx = p_PlotData->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[p_PlotData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< p_PlotData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=p_PlotData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[p_PlotData->GetXaxis()->GetNbins()]=p_PlotData->GetXaxis()->GetBinUpEdge(p_PlotData->GetXaxis()->GetNbins());

  TH1D* h_PlotData   = new TH1D("h_PlotData","",nx,xbins);
  TH1D* h_PlotMCor   = new TH1D("h_PlotMCor","",nx,xbins);
  TH1D* h_PlotMC     = new TH1D("h_PlotMC","",nx,xbins);
  TH1D* h_PlotMCmu   = new TH1D("h_PlotMCmu","",nx,xbins);
  TH1D* h_PlotMCfull = new TH1D("h_PlotMCfull","",nx,xbins);

  for(int b=1; b< (2+(h_PlotData->GetNbinsX())); b++){
    double data=p_PlotData->GetBinError(b);//*std::sqrt(std::fabs(p_PlotData->GetBinContent(b)));
    double mcor=p_PlotMCor->GetBinError(b);//*std::sqrt(std::fabs(p_PlotMCor->GetBinContent(b)));
    double mc=p_PlotMC->GetBinError(b);//*std::sqrt(std::fabs(p_PlotMC->GetBinContent(b)));
    double mcmu=p_PlotMCmu->GetBinError(b);//*std::sqrt(std::fabs(p_PlotMCmu->GetBinContent(b)));
    double mcfull=p_PlotMCfull->GetBinError(b);//*std::sqrt(std::fabs(p_PlotMCfull->GetBinContent(b)));
    std::cout<<data<<"   "<<mc<<"   "<<mcmu<<"   "<<mcfull<<std::endl;
    h_PlotData   ->SetBinContent(b,data);
    h_PlotMCor   ->SetBinContent(b,mcor);
    h_PlotMC     ->SetBinContent(b,mc);
    h_PlotMCmu   ->SetBinContent(b,mcmu);    
    h_PlotMCfull ->SetBinContent(b,mcfull);
  }


  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0., ylow=0., yhigh=0.;
  bool yRange=false;
  TString leglabel;
  if(profile=="uperp_set"){xlabel="#sum E_{T}^{PFO} [GeV]";  ylabel="#sigma (u_{T}) [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=100;}
  if(profile=="bias_set"){xlabel="#sum E_{T}^{PFO} [GeV]"; ylabel="#sigma (u_{||} + p_{T}) [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=100;}
  if(profile=="uperp_pt"){xlabel="p_{T}^{V} [GeV]";  ylabel="#sigma (u_{T}) [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-2.5; yhigh=100;}
  if(profile=="bias_pt"){xlabel="p_{T}^{V} [GeV]"; ylabel="#sigma (u_{||} + p_{T}) [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-7.5; yhigh=100;}

  if(config.WZSelection=="zmumu") leglabel="MC Z#rightarrow#mu#mu";
  if(config.WZSelection=="wminus" || config.WZSelection=="wminus") leglabel="MC W#rightarrow#mu#nu";
  if(nplots=="dataMC") ratiolabel="Data/MC";
  if(nplots=="MC") ratiolabel="Calib/NoCalib";

  plotAxisLine(h_PlotData,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(h_PlotMCor,kGray+2,kGray+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(h_PlotMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(h_PlotMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(h_PlotMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,yRange,ylow,yhigh);
  yRange=false;

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2018")      lumi="59.9";
  if(config.DataYears=="2015+2016") lumi="36.1";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.15,0.85,0.45,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  if(nplots=="dataMC")
    leg->AddEntry(h_PlotData,"#bf{#scale[0.85]{Data "+naming+"}}","PL");
  leg->AddEntry(h_PlotMCor,"#bf{#scale[0.85]{"+leglabel+" Original}}","PL");
  leg->AddEntry(h_PlotMC,"#bf{#scale[0.85]{"+leglabel+"<#mu>}}","PL");
  leg->AddEntry(h_PlotMCmu,"#bf{#scale[0.85]{"+leglabel+"<#mu>+SET+Insitu}}","PL");
  leg->AddEntry(h_PlotMCfull,"#bf{#scale[0.85]{"+leglabel+" Full}}","PL");

  
  TPaveText *box;
  box = new TPaveText(0.2,0.865,0.375,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("13 TeV, "+lumi+" fb^{-1}"); 
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.02);

  h_PlotMC->SetTitle("");
  h_PlotMC->GetXaxis()->SetTitleOffset(1.45);
  h_PlotMC->DrawCopy("hist p");
  //if(nplots=="dataMC")
  h_PlotData->DrawCopy("hist p sames");
  h_PlotMCor->DrawCopy("hist p sames");
  h_PlotMCmu->DrawCopy("hist p sames");
  h_PlotMCfull->DrawCopy("hist p sames");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  TH1D* ratiodataMCor   = new TH1D("ratiodataMCor","",nx,xbins);
  TH1D* ratiodataMC     = new TH1D("ratiodataMC","",nx,xbins);
  TH1D* ratiodataMCmu   = new TH1D("ratiodataMCmu","",nx,xbins);
  TH1D* ratiodataMCfull = new TH1D("ratiodataMCfull","",nx,xbins);
  TH1D* ratioMCmu = new TH1D("ratioMCmu","",nx,xbins);

  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){
    double data=h_PlotData->GetBinContent(b);
    double mc=h_PlotMC->GetBinContent(b);
    if(mc==0) continue;
    double mcor=h_PlotMCor->GetBinContent(b);
    if(mcor==0) continue;
    double mcmu=h_PlotMCmu->GetBinContent(b);
    if(mcmu==0) continue;
    double mcfull=h_PlotMCfull->GetBinContent(b);
    if(mcfull==0) continue;
    double errd=h_PlotData->GetBinError(b);
    
    ratiodataMCor->SetBinContent(b,data/mcor);
    ratiodataMCor->SetBinError(b,errd);

    ratiodataMC->SetBinContent(b,data/mc);
    ratiodataMC->SetBinError(b,errd);

    ratiodataMCmu->SetBinContent(b,data/mcmu);
    ratiodataMCmu->SetBinError(b,errd);

    ratiodataMCfull->SetBinContent(b,data/mcfull);
    ratiodataMCfull->SetBinError(b,0);

    ratioMCmu->SetBinContent(b,mcmu/mc);
    ratioMCmu->SetBinError(b,errd);
  }

  plotAxisLine(ratiodataMCor  ,kGray+2,kGray+2,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMC    ,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMCmu  ,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratiodataMCfull,kViolet,kViolet,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(ratioMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  ratioSettings(ratiodataMCor  ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kGray+2,0.5,10); 
  ratioSettings(ratiodataMC    ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kCyan+3,0.5,10); 
  ratioSettings(ratiodataMCmu  ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 
  ratioSettings(ratiodataMCfull,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kViolet,0.5,10); 
  ratioSettings(ratioMCmu      ,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 

  if(nplots=="dataMC")ratiodataMCor->DrawCopy("hist p");
  if(nplots=="dataMC") ratiodataMC->DrawCopy("hist p same");
  if(nplots=="dataMC") ratiodataMCmu->DrawCopy("hist p same");
  if(nplots=="dataMC") ratiodataMCfull->DrawCopy("hist p same");
  if(nplots=="MC")ratioMCmu->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,1.,xhigh,1.);
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  std::string boson;
  if(config.WZSelection=="zmumu") boson="z";
  if(config.WZSelection=="wminus") boson="wminus";
  if(config.WZSelection=="wplus") boson="wplus";

  profile=profile+"_width";

  c->Print(("./ClosurePlots/closure_"+boson+"_"+profile+"_"+nplots+"_"+config.DataYears+"_profile.pdf").c_str());

  delete ratiodataMC; delete ratiodataMCmu;  
  delete p_PlotData; delete p_PlotMCor; delete p_PlotMC; delete p_PlotMCmu; delete p_PlotMCfull;
  return;
}


void PrintPlots::ResolResponse(Config config, std::string rr, std::string nplots)
{

  // TH1D *p_PlotData, *p_PlotMC, *p_PlotMCmu;
   
  // std::string prof;
  // if(rr=="resp_set_mean"  || rr=="resp_pt_mean" || rr=="resp_set_width"  || rr=="resp_pt_width")
  //   prof = "respcorr";
  
  // if(rr=="resol_set_mean" || rr=="resol_pt_mean" || rr=="resol_set_width" || rr=="resol_pt_width")
  //   prof = "resolcorr";
  
  // TH3D *h3d_data, *h3d_mc, *h3d_mccal;
  // TH2D *h2d_data, *h2d_mc, *h2d_mccal;
  
  // h3d_data  = (TH3D*) fdata->Get((prof+"_vs_set_vs_ptv").c_str());      
  // h3d_mc    = (TH3D*) fmc->Get((prof+"_vs_set_vs_ptv").c_str());   
  // h3d_mccal = (TH3D*) fmccal->Get((prof+"_vs_set_vs_ptv").c_str());
  
  // TProfile2D *rr_data_2d  = (TProfile2D*)h3d_data ->Project3DProfile("p2_data_yz");    
  // TProfile2D *rr_mc_2d    = (TProfile2D*)h3d_mc   ->Project3DProfile("p2_mc_yz");
  // TProfile2D *rr_mccal_2d = (TProfile2D*)h3d_mccal->Project3DProfile("p2_mccal_yz");
  
  // TString option;
  // if(rr=="resp_set_width"  || rr=="resp_pt_width" || rr=="resol_set_width" || rr=="resol_pt_width"){
  //   option = "C=E";
  // }else{
  //   option = "";
  // }
  
  // if(prof=="resolcorr"){ 
  //   rr_data_2d->SetErrorOption("s"); rr_mc_2d->SetErrorOption("s"); rr_mccal_2d->SetErrorOption("s");}
  
  // h2d_data  = (TH2D*)rr_data_2d  ->ProjectionXY("h2d_data",option);
  // h2d_mc    = (TH2D*)rr_mc_2d    ->ProjectionXY("h2d_mc",option);
  // h2d_mccal = (TH2D*)rr_mccal_2d ->ProjectionXY("h2d_mccal",option);
  
    
  // if(rr=="resp_set_mean" || rr=="resp_set_width"){
  //   p_PlotData  = (TH1D*)h2d_data ->ProjectionX();//"pfydata",1,h2d_data->GetNbinsY(),"C=E");
  //   p_PlotMC    = (TH1D*)h2d_mc   ->ProjectionX();//"pfymc",1,h2d_mc->GetNbinsY(),"C=E");
  //   p_PlotMCmu  = (TH1D*)h2d_mccal->ProjectionX();//"pfymccal",1,h2d_mccal->GetNbinsY(),"C=E");
  //   //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  // }
  
  // if(rr=="resol_set_mean" || rr=="resol_set_width"){
  //   p_PlotData  = (TH1D*)h2d_data ->ProjectionX();//"pfydata",1,h2d_data->GetNbinsY(),"C=E");
  //   p_PlotMC    = (TH1D*)h2d_mc   ->ProjectionX();//"pfymc",1,h2d_mc->GetNbinsY(),"C=E");
  //   p_PlotMCmu  = (TH1D*)h2d_mccal->ProjectionX();//"pfymccal",1,h2d_mccal->GetNbinsY(),"C=E");
  //   // p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  // }
  
  // if(rr=="resp_pt_mean" || rr=="resp_pt_width"){
  //   p_PlotData  = (TH1D*)h2d_data ->ProjectionY();//"pfxdata",1,h2d_data->GetNbinsX(),"C=E");
  //   p_PlotMC    = (TH1D*)h2d_mc   ->ProjectionY();//"pfxmc",1,h2d_mc->GetNbinsX(),"C=E");
  //   p_PlotMCmu  = (TH1D*)h2d_mccal->ProjectionY();//"pfxmccal",1,h2d_mccal->GetNbinsX(),"C=E");
  //     //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  // }
  
  // if(rr=="resol_pt_mean" || rr=="resol_pt_width"){
  //   p_PlotData  = (TH1D*)h2d_data ->ProjectionY();//"pfxdata",1,h2d_data->GetNbinsX(),"C=E");
  //   p_PlotMC    = (TH1D*)h2d_mc   ->ProjectionY();//"pfxmc",1,h2d_mc->GetNbinsX(),"C=E");
  //   p_PlotMCmu  = (TH1D*)h2d_mccal->ProjectionY();//"pfxmccal",1,h2d_mccal->GetNbinsX(),"C=E");
  //   //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  // }
  
    
  TProfile *p_PlotData, *p_PlotMC, *p_PlotMCmu;
   
  std::string prof;
  if(rr=="resp_set_mean" || rr=="resp_set_width")
    prof = "h_bias_vs_set";

  if(rr=="resol_set_mean" || rr=="resol_set_width")
    prof = "h_uperp_pfo_vs_set";

  if(rr=="resp_pt_mean" || rr=="resp_pt_width")
    prof = "h_bias_vs_ptv";

  if(rr=="resol_pt_mean" || rr=="resol_pt_width")
    prof = "h_uperp_pfo_vs_pt";
  
  TH2D *h2d_mc, *h2d_mccal, *h2d_data;  
    
  TString option;
  if(rr=="resp_set_width"  || rr=="resp_pt_width" || rr=="resol_set_width" || rr=="resol_pt_width"){
    option = "C=E";
  }else{
    option = "";
  }
  
  //based on simple 2D histos
  h2d_data  = (TH2D*)fdata->Get((prof+"").c_str());	
  h2d_mc    = (TH2D*)fmc->Get((prof+"").c_str());   
  h2d_mccal = (TH2D*)fmccal->Get((prof+"").c_str());
  
  if(rr=="resp_set_mean" || rr=="resp_set_width"){
    p_PlotData  = (TProfile*)h2d_data ->ProfileY("pfydata",1,h2d_data->GetNbinsY(),"");
    p_PlotMC    = (TProfile*)h2d_mc   ->ProfileY("pfymc",1,h2d_mc->GetNbinsY(),"");
    p_PlotMCmu  = (TProfile*)h2d_mccal->ProfileY("pfymccal",1,h2d_mccal->GetNbinsY(),"");
    //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  }
  
  if(rr=="resol_set_mean" || rr=="resol_set_width"){
    p_PlotData  = (TProfile*)h2d_data ->ProfileY("pfydata",1,h2d_data->GetNbinsY(),"");
    p_PlotMC    = (TProfile*)h2d_mc   ->ProfileY("pfymc",1,h2d_mc->GetNbinsY(),"");
    p_PlotMCmu  = (TProfile*)h2d_mccal->ProfileY("pfymccal",1,h2d_mccal->GetNbinsY(),"");
    //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  }
  
  if(rr=="resp_pt_mean" || rr=="resp_pt_width"){
    
    p_PlotData  = (TProfile*)h2d_data ->ProfileY("pfxdata",1,h2d_data->GetNbinsY(),"");
    p_PlotMC    = (TProfile*)h2d_mc   ->ProfileY("pfxmc",1,h2d_mc->GetNbinsY(),"");
    p_PlotMCmu  = (TProfile*)h2d_mccal->ProfileY("pfxmccal",1,h2d_mccal->GetNbinsY(),"");
    //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  }
    
  if(rr=="resol_pt_mean" || rr=="resol_pt_width"){
    p_PlotData  = (TProfile*)h2d_data ->ProfileY("pfxdata",1,h2d_data->GetNbinsY(),"");
    p_PlotMC    = (TProfile*)h2d_mc   ->ProfileY("pfxmc",1,h2d_mc->GetNbinsY(),"");
    p_PlotMCmu  = (TProfile*)h2d_mccal->ProfileY("pfxmccal",1,h2d_mccal->GetNbinsY(),"");
    //p_PlotData->SetErrorOption("s"); p_PlotMC->SetErrorOption("s"); p_PlotMCmu->SetErrorOption("s");
  }
  
  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0., ylow=0., yhigh=0.;
  bool yRange=false;
  TString leglabel;
  // if(rr=="resp_pt_mean" || rr=="resp_pt_width"){ylabel="response [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-10; yhigh=40; p_PlotData->Rebin(5);  p_PlotMC->Rebin(5); p_PlotMCmu->Rebin(5);}
  // if(rr=="resp_set_mean" || rr=="resp_set_width"){ylabel="response [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-7.5; yhigh=50;}
  // if(rr=="resol_pt_mean" || rr=="resol_pt_width"){ylabel="resol [GeV]"; xlow=0.; xhigh=200.; yRange=true; ylow=-7.5; yhigh=50;}
  // if(rr=="resol_set_mean"|| rr=="resol_set_mean"){ylabel="resol [GeV]"; xlow=0.; xhigh=2000.; yRange=true; ylow=-5.5; yhigh=10; p_PlotData->Rebin(20);  p_PlotMC->Rebin(20); p_PlotMCmu->Rebin(20);}

  p_PlotData->Rebin(5);  p_PlotMC->Rebin(5); p_PlotMCmu->Rebin(5);

  for(int b=1; b< (2+p_PlotData->GetNbinsX()); b++){
  std::cout<<p_PlotData->GetBinContent(b)<<"   "<<p_PlotMC->GetBinContent(b)<<"    "<<p_PlotMCmu->GetBinContent(b)<<std::endl;

 }


  if(config.WZSelection=="zmumu") leglabel="MC Z#rightarrow#mu#mu";
  if(config.WZSelection=="wminus" || config.WZSelection=="wminus") leglabel="MC W#rightarrow#mu#nu";
  if(nplots=="dataMC") ratiolabel="Data-MC [GeV]";
  if(nplots=="MC") ratiolabel="Calib-NoCalib [GeV]";

  plotAxisLine(p_PlotData,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,yRange,ylow,yhigh);
  plotAxisLine(p_PlotMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,yRange,ylow,yhigh);

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2018")      lumi="59.9";
  if(config.DataYears=="2015+2016") lumi="36.1";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.65,0.85,0.85,0.7,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  if(nplots=="dataMC")
    leg->AddEntry(p_PlotData,"#bf{#scale[0.85]{Data "+naming+"}}","PL");

  leg->AddEntry(p_PlotMCmu,"#bf{#scale[0.85]{"+leglabel+" Calibrated}}","PL");
  
  leg->AddEntry(p_PlotMC,"#bf{#scale[0.85]{"+leglabel+"}}","PL");
  
  TPaveText *box;
  box = new TPaveText(0.65,0.865,0.825,0.965,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("13 TeV, "+lumi+" fb^{-1}"); 
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.01);

  p_PlotMC->SetTitle("");
  p_PlotMC->GetXaxis()->SetTitleOffset(1.45);
  p_PlotMC->DrawCopy("e");  
  if(nplots=="dataMC") p_PlotData->DrawCopy("e sames");
  p_PlotMCmu->DrawCopy("e sames");

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = p_PlotData->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[p_PlotData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< p_PlotData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=p_PlotData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[p_PlotData->GetXaxis()->GetNbins()]=p_PlotData->GetXaxis()->GetBinUpEdge(p_PlotData->GetXaxis()->GetNbins());

  TH1D* diffdataMC   = new TH1D("diffdataMC","",nx,xbins);
  TH1D* diffdataMCmu = new TH1D("diffdataMCmu","",nx,xbins);
  TH1D* diffMCmu = new TH1D("diffMCmu","",nx,xbins);

  for(int b=1; b< (2+(diffdataMC->GetNbinsX())); b++){
    double data=p_PlotData->GetBinContent(b);
    double mc=p_PlotMC->GetBinContent(b);
    double mcmu=p_PlotMCmu->GetBinContent(b);
    double errd=p_PlotData->GetBinError(b);
    
    diffdataMC->SetBinContent(b,data-mc);
    diffdataMC->SetBinError(b,errd);

    diffdataMCmu->SetBinContent(b,data-mcmu);
    diffdataMCmu->SetBinError(b,errd);

    diffMCmu->SetBinContent(b,mcmu-mc);
    diffMCmu->SetBinError(b,errd);
  }

  plotAxisLine(diffdataMC,kCyan+3,kCyan+3,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(diffdataMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);
  plotAxisLine(diffMCmu,kBlue,kBlue,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,0);

  ratioSettings(diffdataMC,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kCyan+3,0.5,10); 
  ratioSettings(diffdataMCmu,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 
  ratioSettings(diffMCmu,-0.5,0.5,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,23,kBlue,0.5,10); 

  if(nplots=="dataMC")diffdataMC->DrawCopy("p ");
  if(nplots=="dataMC") diffdataMCmu->DrawCopy("p same");
  if(nplots=="MC")diffMCmu->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,0.,xhigh,0.);
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  std::string boson;
  if(config.WZSelection=="zmumu")  boson="z";
  if(config.WZSelection=="wminus") boson="wminus";
  if(config.WZSelection=="wplus")  boson="wplus";

  c->Print(("./ClosurePlots/closure_"+boson+"_"+rr+"_"+nplots+"_"+config.DataYears+"_rr.pdf").c_str());

  delete diffdataMC; delete diffdataMCmu;
  delete p_PlotData; delete p_PlotMC; delete p_PlotMCmu;
  return;
}


void PrintPlots::WStackPlot(Config config, std::string kine, std::string option)
{
  TH1D *hData, *hWmunu, *hWtaunu, *hZmumu, *hZtautau, *hTtbar, *sum;

  if(option=="calib"){
    hData    = (TH1D*)fdata  ->GetObjectUnchecked(("h_"+kine).c_str());
    hWmunu   = (TH1D*)fmcfull->GetObjectUnchecked(("h_"+kine).c_str());
    sum      = (TH1D*)hWmunu->Clone("sum");  
    hWtaunu  = (TH1D*)wtaunu ->GetObjectUnchecked(("h_"+kine).c_str());
    hZmumu   = (TH1D*)zmumu  ->GetObjectUnchecked(("h_"+kine).c_str());
    hZtautau = (TH1D*)ztautau->GetObjectUnchecked(("h_"+kine).c_str());
    hTtbar   = (TH1D*)ttbar  ->GetObjectUnchecked(("h_"+kine).c_str());

  }else if(option=="control"){
    hData    = (TH1D*)wzdata  ->GetObjectUnchecked(("h_"+kine).c_str());
    hWmunu   = (TH1D*)wzmcfull->GetObjectUnchecked(("h_"+kine).c_str());
    sum      = (TH1D*)hWmunu->Clone("sum");
    hWtaunu  = (TH1D*)pwtaunu ->GetObjectUnchecked(("h_"+kine).c_str());
    hZmumu   = (TH1D*)pzmumu  ->GetObjectUnchecked(("h_"+kine).c_str());
    hZtautau = (TH1D*)pztautau->GetObjectUnchecked(("h_"+kine).c_str());
    hTtbar   = (TH1D*)pttbar  ->GetObjectUnchecked(("h_"+kine).c_str());
  }

  if(option=="calib" && kine!="Sumet"){ hData->Rebin(4); hWmunu->Rebin(4); hWtaunu->Rebin(4); hZmumu->Rebin(4); hZtautau->Rebin(4); hTtbar->Rebin(4);}
  if(option=="calib" && kine=="Sumet"){ hData->Rebin(10); hWmunu->Rebin(10); hWtaunu->Rebin(10); hZmumu->Rebin(10); hZtautau->Rebin(10); hTtbar->Rebin(10);}

  std::string xlabel,ylabel,ratiolabel; double xlow=0., xhigh=0., yhigh=0.;
  bool yRange=false;
  TString leglabel;

  if(config.NormPlots=="True"){
    hData   ->Scale(1./hData   ->Integral());
    hWmunu  ->Scale(1./hWmunu  ->Integral());
    hWtaunu ->Scale(1./hWtaunu ->Integral());
    hZmumu  ->Scale(1./hZmumu  ->Integral());
    hZtautau->Scale(1./hZtautau->Integral());
    hTtbar  ->Scale(1./hTtbar  ->Integral());
    ylabel="A.U.";
  }else
    ylabel="Events";

  if(kine=="mu"){xlabel="#LT#mu#GT"; xlow=0; xhigh=50;}
  if(kine=="Sumet"){xlabel="#sum E_{T}^{PFO} [GeV]"; xlow=0.; xhigh=2000;}
  if(kine=="ptv"){xlabel="p_{T}^{V} [GeV]"; xlow=0.; xhigh=200;}
  if(kine=="upar_rec"){xlabel="u_{par}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="upar_ptv"){xlabel="u_{par}^{V}+p_{T}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="upar_ptv_n"){xlabel="u_{par}^{V}-p_{T}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="uperp_rec"){xlabel="u_{perp}^{V} [GeV]"; xlow=-200.; xhigh=200.;}
  if(kine=="ut"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200.;}
  if(kine=="ux"){xlabel="u_{x} [GeV]"; xlow=-100.; xhigh=100.;}
  if(kine=="uy"){xlabel="u_{y} [GeV]"; xlow=-100.; xhigh=100.;}
  if(kine=="utphi"){xlabel="#phi(u_{T}) [rad]"; xlow=-5.; xhigh=5.; yRange=true; yhigh=9500000.;}
  if(kine=="met"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200.;}
  if(kine=="metphi"){xlabel="#phi(E_{T}^{miss})"; xlow=-5.; xhigh=5.; yRange=true; yhigh=9500000.;}

  if(kine=="met" && option =="control"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="rmass"){xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=150000.;}
  if(kine=="SumET"){xlabel="#sum E_{T}^{PFO} [GeV]"; xlow=0.; xhigh=2000*1000;}
  if(kine=="met_recoil" && option =="control"){xlabel="E_{T}^{miss} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="rmass_recoil"){xlabel="m_{T} [GeV]"; xlow=40000.; xhigh=150000.;}
  if(kine=="pt_boson_reco"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200*1000.;}
  if(kine=="pt_boson_reco_p"){xlabel="u_{T} [GeV]"; xlow=0.; xhigh=200*1000.;}

  plotAxisLine(hData,kBlack,kBlack,8,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hData ->SetFillColor(kWhite);
  plotAxisLine(hWmunu,kBlack,kBlack,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hWmunu->SetFillColor(kWhite);
  plotAxisLine(hWtaunu,kBlack,kBlack,34,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hWtaunu->SetFillColor(kAzure);
  plotAxisLine(hZmumu,kBlack,kBlack,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hZmumu->SetFillColor(kOrange+1);
  plotAxisLine(hZtautau,kBlack,kBlack,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hZtautau->SetFillColor(kGreen+2);
  plotAxisLine(hTtbar,kBlack,kBlack,23,0.5,"",xlabel.c_str(),ylabel.c_str(),true,xlow,xhigh,false,0.,yhigh);
  hTtbar->SetFillColor(kRed+2);

  THStack *hs = new THStack("hs","");
  hs ->Add(hTtbar);
  hs ->Add(hZtautau);
  hs ->Add(hZmumu);
  hs ->Add(hWtaunu);
  hs ->Add(hWmunu); 

  sum->Add(hWtaunu); sum->Add(hZmumu); sum->Add(hZtautau); sum->Add(hTtbar);

  if(config.WZSelection=="zmumu") leglabel="MC ";
  if(config.WZSelection=="wminus" || config.WZSelection=="wminus") leglabel="MC W#rightarrow#mu#nu";

  ratiolabel="Data/MC"; 

  TString lumi;
  if(config.DataYears=="2015")      lumi="3.2";
  if(config.DataYears=="2016")      lumi="32.8";
  if(config.DataYears=="2017")      lumi="32.8";
  if(config.DataYears=="2018")      lumi="59.9";
  if(config.DataYears=="2015+2016") lumi="36.1";

  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.65,0.825,0.85,0.675,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  leg->AddEntry(hData,"#bf{#scale[0.85]{Data "+naming+"}}","P");
  leg->AddEntry(hWmunu,"#bf{#scale[0.85]{MC W#rightarrow#mu#nu}}","f");
  leg->AddEntry(hWtaunu,"#bf{#scale[0.85]{MC W#rightarrow#tau#nu}}","f");  
  leg->AddEntry(hZmumu,"#bf{#scale[0.85]{Z#rightarrow#mu#mu}}","f");
  leg->AddEntry(hZtautau,"#bf{#scale[0.85]{Z#rightarrow#tau#tau}}","f");  
  leg->AddEntry(hTtbar,"#bf{#scale[0.85]{ttbar}}","f");  
  
  TPaveText *box;
  box = new TPaveText(0.65,0.845,0.825,0.945,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("13 TeV, "+lumi+" fb^{-1}");
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.15,1,1);
  pad->SetTicks(1,1); 
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.035);

  hData->SetTitle("");
  hData->GetXaxis()->SetTitleOffset(1.45);
  //sum->Draw();
  hs->Draw("hist");
  hData->DrawCopy("pe X0 SAME");
  hs->GetXaxis()->SetTitle(xlabel.c_str());
  hs->GetYaxis()->SetTitle(ylabel.c_str());
  hs->GetYaxis()->SetTitleOffset(1.55);
  // hs->Draw("hist");
  // hs->GetYaxis()->SetLimits(0.,yhigh);

  pad->Update();
  pad->Modified();
  box->Draw();
  leg->Draw();
  c->cd();

  //ratio panel
  TPad *padr = new TPad("padr","padr",0,0.0,1,0.2);
  padr->SetTopMargin(0);
  padr->SetBottomMargin(0.2);
  padr->SetLeftMargin(0.1125);
  padr->SetRightMargin(0.03);
  padr->SetGridy();
  padr->Draw();
  padr->cd();

  unsigned int nx = hData->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[hData->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hData->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hData->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hData->GetXaxis()->GetNbins()]=hData->GetXaxis()->GetBinUpEdge(hData->GetXaxis()->GetNbins());

  TH1D* ratiodataMC = new TH1D("ratiodataMC","",nx,xbins);

  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){

    double data=hData->GetBinContent(b);
    double mc=sum->GetBinContent(b);
    if(mc==0) continue;

    double dataratio=data/mc;
    double err=hData->GetBinError(b)/mc;

    ratiodataMC->SetBinContent(b,dataratio);
    ratiodataMC->SetBinError(b,err);
  }

  ratioSettings(ratiodataMC,0.95,1.05,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,34,kBlack,0.5,10);
  ratiodataMC->DrawCopy("");

  TLine *line;
  line = new TLine(xlow,1.,xhigh,1.);
  //line = new TLine(40000.,1.,150000.,1.);//remove
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  c->cd();
  std::string boson;
  if(config.WZSelection=="zmumu") boson="z";
  if(config.WZSelection=="wminus") boson="wminus";
  if(config.WZSelection=="wplus") boson="wplus";

  c->Print(("./StackPlots/stack_"+boson+"_"+kine+"_"+config.DataYears+".pdf").c_str());
  
  //  delete hData; delete hWmunu; delete hWtaunu; delete hZmumu; delete hZtautau; delete hTtbar; delete sum;
  delete hs;
  return;
}



void PrintPlots::plotAxisLine(TH1D* hist, int lineColor, int markerColor,
			      int markerStyle, double markerSize,
			      TString title, TString xlabel, TString ylabel, bool xRange,
			      double xlow, double xhigh, bool yRange, double ylow, double yhigh)
{
  hist->SetLineColor(lineColor);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerColor(markerColor);
  hist->SetMarkerSize(markerSize);
  hist->SetTitle(title);
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  std::cout<<ylow<<"    "<<yhigh<<std::endl;
  if(xRange==true)
    hist->GetXaxis()->SetRangeUser(xlow,xhigh);
  if(yRange==true)
    hist->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}


void PrintPlots::profAxisLine(TProfile* prof, int lineColor, int markerColor,
			      int markerStyle, double markerSize,
			      TString title, TString xlabel, TString ylabel, bool xRange,
			      double xlow, double xhigh, bool yRange, double ylow, double yhigh)
{
  prof->SetLineColor(lineColor);
  prof->SetMarkerStyle(markerStyle);
  prof->SetMarkerColor(markerColor);
  prof->SetMarkerSize(markerSize);
  prof->SetTitle(title);
  prof->GetXaxis()->SetTitle(xlabel);
  prof->GetYaxis()->SetTitle(ylabel);
  if(xRange==true)
    prof->GetXaxis()->SetRangeUser(xlow,xhigh);
  if(yRange==true)
    prof->GetYaxis()->SetRangeUser(ylow,yhigh);

  return;
}


void PrintPlots::ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
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


void PrintPlots::profRatioSettings(TProfile* prof, double min, double max, TString xlabel, TString ylabel,
			       double xlabelsize, double ylabelsize, double ytitleof, double ytitlesize,
			       int markerstyle, int color, double markersize, int divisions)
{
  prof->SetMinimum(min);
  prof->SetMaximum(max);
  prof->GetXaxis()->SetTitle(xlabel);
  prof->GetYaxis()->SetTitle(ylabel);
  prof->GetXaxis()->SetLabelSize(xlabelsize);
  prof->GetYaxis()->SetLabelSize(ylabelsize);
  prof->GetYaxis()->SetTitleOffset(ytitleof);
  prof->GetYaxis()->SetTitleSize(ytitlesize);
  prof->SetLineColor(color);

  prof->Sumw2();
  prof->SetMarkerStyle(markerstyle);
  prof->SetMarkerSize(markersize);
  prof->SetMarkerColor(color);

  prof->SetStats(0);
  prof->GetYaxis()->SetNdivisions(divisions);

  return;
}

void PrintPlots::setstyle(){

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

#endif
