//This class is to sum diferent numbers row by row in several text files,
//handles mass slices and writes a latex cutflow table
//big update 20/03/2019
//author: Andres Ramirez Morales <andres.ramirez.morales@cern.ch>
#ifndef CUTFLOWMAKER_CXX
#define CUTFLOWMAKER_CXX

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "CutFlowMaker.h"


CutFlowMaker::CutFlowMaker()
{
}

CutFlowMaker::~CutFlowMaker(){}


void CutFlowMaker::initialize(Config config){

  if(config.WZSelection=="zmumu"){
    wzchannel="z";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus";}
  
  std::string data="", total="", dirInclusive="", puname="";
  if(config.DataYears=="2015+2016") {dirYear="2015p2016/"; data="data15p16";}
  if(config.DataYears=="2017") {dirYear="2017/"; data="data17";}
  if(config.DataYears=="2018") {dirYear="2018/"; data="data18";}

  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}

  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
  std::string  systematic="_nominal";

  fileDir=config.OutputFileDir;
   
  //load files
  if(config.WZSelection=="wplus"){
    fsig   = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "wplusmunu" + "_" + wzchannel + systematic + puname + total+".root").c_str());
    fwtau  = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "wplustaunu"+ "_" + wzchannel + systematic + puname + total+".root").c_str());
  }else if(config.WZSelection=="wminus"){
    fsig   = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "wminmunu" + "_" + wzchannel + systematic + puname + total+".root").c_str());
    fwtau  = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "wmintaunu"+ "_" + wzchannel + systematic + puname + total+".root").c_str());  
  }
  fzmu   = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "zmumu"     + "_" + wzchannel + systematic + puname + total+".root").c_str());
  fztau  = new TFile((config.OutputFileDir +"Files/"+ dirYear + dirInclusive + "ztautau"+ "_" + wzchannel + systematic + puname + total+".root").c_str());
  ftop   = new TFile((config.OutputFileDir +"Files/"+ dirYear +  "top"     + "_" + wzchannel + systematic + puname + ".root").c_str());
  fdibo  = new TFile((config.OutputFileDir +"Files/"+ dirYear +  "diboson" + "_" + wzchannel + systematic + puname + ".root").c_str());
  fdata  = new TFile((config.OutputFileDir +"Files/"+ dirYear +  data  + "_" + wzchannel +".root").c_str());

  if(!fsig->IsOpen() || !fwtau->IsOpen() || !fzmu->IsOpen() || !fztau->IsOpen() || !ftop->IsOpen() ||  !fdibo->IsOpen() || !fdata->IsOpen()){
    std::cout<<"Files missing for the cutflow, bye!"<<std::endl; exit(10);}
    
  return;
}

void CutFlowMaker::execute(){

  FetchValues("");//peak(inclusive) region
  latexConvertor(totalFlow, names,"peak");

  totalFlow.clear(); names.clear();

  FetchValues("_high");
  latexConvertor(totalFlow, names,"high");

  return;
}


void CutFlowMaker::FetchValues(std::string option){

  std::vector<double> data, signal, wtaunu, zmumu, ztautau, top, diboson;
  data.clear(); signal.clear(); wtaunu.clear(); zmumu.clear(); ztautau.clear(); top.clear(); diboson.clear();

  TH1D *hData, *hSignal, *hWtaunu, *hZmumu, *hZtau, *hTop, *hDibo;
  hData     = (TH1D*)fdata ->Get(("h_cutflow"+option).c_str());
  hSignal   = (TH1D*)fsig  ->Get(("h_cutflow"+option).c_str());
  hWtaunu   = (TH1D*)fwtau ->Get(("h_cutflow"+option).c_str());
  hZmumu    = (TH1D*)fzmu  ->Get(("h_cutflow"+option).c_str());
  hZtau     = (TH1D*)fztau ->Get(("h_cutflow"+option).c_str());
  hTop      = (TH1D*)ftop  ->Get(("h_cutflow"+option).c_str());
  hDibo     = (TH1D*)fdibo ->Get(("h_cutflow"+option).c_str());

  int nx = hData->GetXaxis()->GetNbins();
  int lastCuts = 0;//exclude the last cuts of the hist if needed
  for(int k=1; k<(nx+1)-lastCuts; k++){//
    if(k==5 || k==6) continue;//these are useless cuts(?)
    if(hSignal->GetBinContent(k)==0) continue;
    signal.push_back(hSignal->GetBinContent(k));
    wtaunu.push_back(hWtaunu->GetBinContent(k));
    zmumu.push_back(hZmumu->GetBinContent(k));
    ztautau.push_back(hZtau->GetBinContent(k));
    top.push_back(hTop->GetBinContent(k));
    diboson.push_back(hDibo->GetBinContent(k));
  }

  for(int k=1; k<(nx+1)-lastCuts; k++){//for 2017 run number cuts are not displayed in the h_cutflow
    if(hData->GetBinContent(k)==0) continue;
    data.push_back(hData->GetBinContent(k));
  }
  
  totalFlow.push_back(data); totalFlow.push_back(signal); totalFlow.push_back(wtaunu); totalFlow.push_back(zmumu);
  totalFlow.push_back(ztautau); totalFlow.push_back(top); totalFlow.push_back(diboson);
  
  //names
  names.push_back("data");
  if(wzchannel=="wplus"){
    names.push_back("wplusmunu"); names.push_back("wplustaunu");
  }else if(wzchannel=="wminus"){
    names.push_back("wminmunu"); names.push_back("wmintaunu");}
  names.push_back("zmumu");
  names.push_back("ztautau");
  names.push_back("top");
  names.push_back("diboson");

  return;
}


void CutFlowMaker::latexConvertor(std::vector<std::vector</*long int*/double>> totalValues, std::vector<std::string> names, std::string option){
    
  std::vector<std::string> cutName; cutName.clear(); cutName=Names(option);

  if(names.size()!=totalValues.size()) return;

  int nNames = (int)names.size(); int nValues = 0;
  nValues = (int)totalValues.at(1).size()-1;
  std::ofstream cutFlowLatex((fileDir + "CutFlow/"+dirYear+"All_sum_"+wzchannel+"_"+option+".tex").c_str());

  double sumTotalMC=0.;  
  for(int k=1;k<nNames;k++) sumTotalMC+=totalValues.at(k).at(nValues-1);

  //begin of the table
  cutFlowLatex << "\\documentclass[12pt]{book}" << std::endl;
  cutFlowLatex << "\\usepackage{lscape}" << std::endl;
  cutFlowLatex << "\\setlength{\\topmargin}{1.25 in}" << std::endl;

  cutFlowLatex << "\\usepackage{siunitx}"<<std::endl;
  cutFlowLatex << "\\sisetup{exponent-product=\\ensuremath{{}\\cdot{}}}"<<std::endl;
  cutFlowLatex << "\\sisetup{round-mode=places,round-precision=2}"<<std::endl;
  
  cutFlowLatex << "\\begin{document}" << std::endl;
  cutFlowLatex << "\\begin{landscape}" << std::endl;

  cutFlowLatex << "\\begin{table}[h!]" << std::endl;
  //cutFlowLatex << "\\begin{center}" << std::endl;
  cutFlowLatex << "\\scriptsize{" << std::endl;
  cutFlowLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nNames; k++) cutFlowLatex<<" c "; 
  cutFlowLatex<<" }\\hline \\hline " << std::endl;

  cutFlowLatex << " & ";
  for(int k=0; k < nNames; k++){
    cutFlowLatex<<std::setw(4)<<std::right;
    names.at(k)=processName(names.at(k));
    if(k< nNames-1) cutFlowLatex<<" "<<names.at(k)<<" &";
    if(k==nNames-1) cutFlowLatex<<" "<<names.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nValues; j++){
    cutFlowLatex<<std::setw(12)<<std::right;
    cutFlowLatex<<cutName.at(j)<<" & ";
    for(int l=0;l < nNames; l++){
      cutFlowLatex<<std::setw(12)<<std::right;
      if(l< nNames-1){
	cutFlowLatex<<"  \\num{"<<(totalValues.at(l)).at(j)<<"}";
	if(j==0 && true){ cutFlowLatex<<"  (NA)       ";}
	if(j!=0 && true){ cutFlowLatex<<std::setprecision(4)<<"  ("<<((totalValues.at(l)).at(j)/(totalValues.at(l)).at(j-1))*100.<<"\\%)  ";}
	cutFlowLatex<<std::right;
	if(j!=0 && std::ceil(((totalValues.at(l)).at(j)/(totalValues.at(l)).at(j-1))*100000)==100000) cutFlowLatex<<"  ";//very ugly way to have nice formated table
	cutFlowLatex<<" & ";
      }
      if(l==nNames-1){
	cutFlowLatex<<"  \\num{"<<(totalValues.at(l)).at(j)<<"}";
	if(j==0) cutFlowLatex<<" (NA) ";
	if(j!=0) cutFlowLatex<<"  ("<<std::setprecision(4)<<((totalValues.at(l)).at(j)/(totalValues.at(l)).at(j-1))*100.<<"\\%) ";
	cutFlowLatex<<" \\\\"<<std::endl;
      }
    }
  }
  
  cutFlowLatex << "\\hline\\hline" << std::endl;

  cutFlowLatex << "Total Efficiency & ";
  for(int k=0; k < nNames; k++){
    cutFlowLatex<<std::setw(4)<<std::right;
    if(k< nNames-1) cutFlowLatex<<std::setprecision(2)<<"  "<<((totalValues.at(k)).at(nValues-1)/(totalValues.at(k)).at(0))*100<<"\\%  &";
    if(k==nNames-1) cutFlowLatex<<std::setprecision(2)<<"  "<<((totalValues.at(k)).at(nValues-1)/(totalValues.at(k)).at(0))*100<<"\\%  \\\\"<<std::endl;
  }

  cutFlowLatex << "\\hline\\hline" << std::endl;

  cutFlowLatex << "ProcessMC/TotalMC & ";
  for(int k=0; k < nNames; k++){
    cutFlowLatex<<std::setw(4)<<std::right;
    if(k==0) cutFlowLatex<<std::setprecision(4)<<"  NA  &";
    if(k< nNames-1 && k!=0) cutFlowLatex<<std::setprecision(4)<<"  "<<((totalValues.at(k)).at(nValues-1)/sumTotalMC)*100<<"\\%  &";
    if(k==nNames-1) cutFlowLatex<<std::setprecision(4)<<"  "<<((totalValues.at(k)).at(nValues-1)/sumTotalMC)*100<<"\\%  \\\\"<<std::endl;
  }

  cutFlowLatex << "  Data/ProcessMC & ";
  for(int k=0; k < nNames; k++){
    cutFlowLatex<<std::setw(4)<<std::right;
    if(k==0) cutFlowLatex<<std::setprecision(4)<<"  NA  &";
    if(k< nNames-1 && k!=0) cutFlowLatex<<std::setprecision(4)<<"  "<<((totalValues.at(k)).at(nValues-1)/(totalValues.at(0)).at(nValues-1))*100<<"\\%  &";
    if(k==nNames-1) cutFlowLatex<<std::setprecision(4)<<"  "<<((totalValues.at(k)).at(nValues-1)/(totalValues.at(0)).at(nValues-1))*100<<"\\%  \\\\"<<std::endl;
  }

  cutFlowLatex << "\\hline" << std::endl;
  cutFlowLatex << " Data/TotalMC &  "<<std::setprecision(4)<<((totalValues.at(0)).at(nValues-1)/sumTotalMC)*100.<<"\\%  \\\\"<<std::endl;

  cutFlowLatex << "\\hline \\hline" << std::endl;
  cutFlowLatex << "\\end{tabular}" << std::endl;
  cutFlowLatex << "}" << std::endl;//scriptsize
  //cutFlowLatex << "\\end{center}" << std::endl;
  cutFlowLatex << "\\caption{Cutflow \\label{cut}}" << std::endl;  
  cutFlowLatex << "\\end{table}" << std::endl;
  cutFlowLatex << "\\end{landscape}" << std::endl;
  cutFlowLatex << "\\end{document}" << std::endl;

  cutFlowLatex.close();

  std::cout<<"CutFlowMaker done! Please find the cutflow files in /yourDirectory/CutFlow/DataYears/"<<std::endl;
  
  return;

}


std::string CutFlowMaker::processName(std::string inputName){

  std::string stylizedName;

  if(inputName=="data")  stylizedName="Data";

  if(inputName=="wminmunu")  stylizedName="$W^{-}\\rightarrow\\mu\\nu$";
  if(inputName=="wmintaunu") stylizedName="$W^{-}\\rightarrow\\tau\\nu$";

  if(inputName=="wplusmunu")  stylizedName="$W^{+}\\rightarrow\\mu\\nu$";
  if(inputName=="wplustaunu") stylizedName="$W^{+}\\rightarrow\\tau\\nu$";

  if(inputName=="zmumu")  stylizedName="$Z\\rightarrow\\mu\\mu$";
  if(inputName=="ztautau") stylizedName="$Z\\rightarrow\\tau\\tau$";

  if(inputName=="top")  stylizedName="Top";
  if(inputName=="diboson") stylizedName="Diboson";

  return stylizedName;
}


std::vector<std::string> CutFlowMaker::Names(std::string option){

  std::vector<std::string> names;

  names.push_back("Initial                     ");                  
  names.push_back("GRL                         ");		    
  names.push_back("Good Calo                   ");		    
  names.push_back("PriVtx                      ");		    
  names.push_back("Trigger                     ");		      
  names.push_back("MU\\_N $30>=1$                  ");		      
  names.push_back("Trigger match                 ");		    
  names.push_back("JetCleaning: Loose Bad      ");		    
  names.push_back("MU\\_N == 1                  ");		    
  names.push_back("EL\\_N == 0                  ");		    
  names.push_back("$E_{T}^{miss} > 25$ GeV     ");		    
  names.push_back("$m_{T} > 55$ GeV            ");		    
  names.push_back("Save                        ");		    
  names.push_back("NOBADMUON                   ");		    
  if(option=="peak"){						    
  names.push_back("Lepton Veto                 ");		    
  names.push_back("Charge  Selection           ");		     
  names.push_back("$\\vert\\eta\\vert^{\\mu} < 2.4$");		     
  names.push_back("$p_{T}^{\\mu} > 30$ GeV      ");
  names.push_back("$E_{T}^{miss} > 25$ GeV     ");
  names.push_back("$m_{T} > 55$ GeV            ");
  }else if(option=="high"){
  names.push_back("Lepton Veto                 ");
  names.push_back("Charge  Selection           ");
  names.push_back("$\\vert\\eta\\vert^{\\mu} < 2.4$");
  names.push_back("$p_{T}^{\\mu} > 55$ GeV      ");
  names.push_back("$E_{T}^{miss} > 55$ GeV     ");
  names.push_back("$m_{T} > 110$ GeV           ");
  }

  return names;
}


//The functions sumFiles and LoadFiles are not used anymore, but they are nice code ;)
void CutFlowMaker::sumFiles(std::vector<std::vector<std::string>> samplesRun, std::string option){

  for(int iSamples=0; iSamples < (int)samplesRun.size();iSamples++){

    std::ifstream inputFile[50];
    std::ofstream cutFlowText(( fileDir + "CutFlow/"+dirYear+samplesRun[iSamples][0]+"_Total_"+wzchannel+"_"+option+".txt").c_str());

    double sum[50][50], temp=0;
    int contador=0;
    std::string nominal="_nominal";
    if(samplesRun[iSamples][0]=="data18" || samplesRun[iSamples][0]=="data17" || samplesRun[iSamples][0]=="data16" || samplesRun[iSamples][0]=="data15") nominal="";
    
    for(int iMassRange=0; iMassRange < (int)samplesRun[iSamples].size(); iMassRange++){
      
      inputFile[iMassRange].open((fileDir + "CutFlow/"+dirYear+samplesRun[iSamples][iMassRange]+"_"+wzchannel+nominal+"_"+option+".txt").c_str());

      int l = 0;
      bool openFile = false;
      if(inputFile[iMassRange].is_open()){
	openFile=true;
	while(inputFile[iMassRange].good()){
	  inputFile[iMassRange]>>temp;
	  sum[iMassRange][l]=temp;
	  //std::cout<<temp<<"    "<<contador<<std::endl;
	  if(iMassRange==0) contador++;
	  l++;	  
	}
      }else std::cout<<"Unable to open the file"<<std::endl;
      
      if(openFile) inputFile[iMassRange].close();
     
    }//mass range loop

    std::cout<<((samplesRun[iSamples][0]+"_Total_"+wzchannel).c_str())<<std::endl;
    std::cout<<"Number of cutFlow entries:   "<<contador-1<< "     Number of mass slices/samples for process:   "<<(int)samplesRun[iSamples].size()<<std::endl;
    temp=0;
    for(int ll=0;ll<contador-1;ll++){
      for(int kk = 0;kk<(int)samplesRun[iSamples].size();kk++){
	//std::cout<<sum[kk][ll]<<std::endl;
	temp=temp+sum[kk][ll];
      }
      cutFlowText<<temp<<std::endl;//std::cout<<temp<<std::endl;
      temp=0;  
    }  
    cutFlowText.close();

  }//samples loop
  
  return;
}


void CutFlowMaker::loadValues(std::vector<std::vector<std::string>> samplesRun, std::string option){

  for(int iSamples=0; iSamples < (int)samplesRun.size();iSamples++){ 
    
    std::ifstream sumFile((fileDir + "CutFlow/"+dirYear+samplesRun[iSamples][0]+"_Total_"+wzchannel+"_"+option+".txt").c_str());
    int k=0; /*long int*/double value=0;
    std::vector</*long int*/double> tempVec; tempVec.clear();
    bool openFile = false;
    if(sumFile.is_open()){
      openFile=true;
      while(sumFile.good()){
	sumFile>>value;
	tempVec.push_back(value);
	k++; value=0.;
      }
    }else std::cout<<"Unable to open the file"<<std::endl;
        
    totalFlow.push_back(tempVec); tempVec.clear();
    if(openFile) sumFile.close();
    
    names.push_back(samplesRun[iSamples][0]);
    
  }//samples loop
  
  return;  
}

#endif
