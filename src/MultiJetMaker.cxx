//This class is to calculate the multijet background
//using the fraction fit methodology
//author: Andres Ramirez Morales <andres.ramirez.morales@cern.ch>
#ifndef MULTIJETMAKER_CXX
#define MULTIJETMAKER_CXX

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "MultiJetMaker.h"
#include "TCanvas.h"

MultiJetMaker::MultiJetMaker()
{
}

MultiJetMaker::~MultiJetMaker(){}


void MultiJetMaker::initialize(Config config){

  if(config.WZSelection=="zmumu"){
    wzchannel="z";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus";}
  
  std::string year="";
  if(config.DataYears=="2015+2016") {dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017") {dirYear="2017/"; year="17";}
  if(config.DataYears=="2018") {dirYear="2018/"; year="18";}

  std::string dirInclusive="", total="";

  std::string systematic = "";
  if(config.MultijetVariation=="Nominal") systematic="25";
  else if(config.MultijetVariation=="30") systematic="30"; //(config.MultijetVariation=="d0up") systematic="30";
  else if(config.MultijetVariation=="20") systematic="20";
  else if(config.MultijetVariation=="xSecup") systematic="xSecup";
  else if(config.MultijetVariation=="xSecdown") systematic="xSecdown";

  fdata0 = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "data" + year + "_" + wzchannel + "_multijet_signal.root").c_str());					    
  fmc0   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "background_"   + wzchannel + "_nominal_nonewPU_multijet_signal" + total +"_" + systematic+".root").c_str());

  fdata0_p = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "data" + year + "_" + wzchannel + "_multijet_signal.root").c_str());					    
  fmc0_p   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "background_"   + wzchannel + "_nominal_nonewPU_multijet_signal" + total + "_"+ systematic+".root").c_str());

  fdata1 = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "data" + year + "_" + wzchannel + "_multijet_fr1.root").c_str());					    
  fmc1   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "background_"   + wzchannel + "_nominal_nonewPU_multijet_fr1" + total + "_"+ systematic+".root").c_str());
  
  fdata2 = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "data" + year + "_" + wzchannel + "_multijet_fr2.root").c_str());					    
  fmc2   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "background_"   + wzchannel + "_nominal_nonewPU_multijet_fr2" + total + "_"+ systematic+".root").c_str());
  
  fout   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "multijet_" + wzchannel + total + "_"+ systematic+".root").c_str(), "RECREATE");
  
  if(wzchannel=="wminus") scale = 1;//debug
  else if (wzchannel=="wplus") scale = 1;

  if(!fdata0->IsOpen() || !fmc0->IsOpen() || !fdata0_p->IsOpen() || !fmc0_p->IsOpen() || !fdata1->IsOpen() || !fmc1->IsOpen() || !fdata2->IsOpen() || !fmc2->IsOpen() || !fout->IsOpen()){
    std::cout<<"The input is incomplete for MultijetMaker, try again. Ciao!"<<std::endl;
    exit(1);}

  return;
}


void MultiJetMaker::execute(Config config){

  std::vector<std::string> kine; kine.clear();
  kine.push_back("rmass");
  kine.push_back("met");
  kine.push_back("d0sig");
  kine.push_back("pt");
  kine.push_back("eta");
  kine.push_back("iso30");
  kine.push_back("eta_d");
  kine.push_back("phi");
  kine.push_back("phi_met");
  kine.push_back("delta_phi");
  kine.push_back("mu");

  std::vector<std::string> MJdisct; MJdisct.clear();
  MJdisct.push_back("0"); MJdisct.push_back("1"); MJdisct.push_back("2");
  MJdisct.push_back("3"); MJdisct.push_back("4"); MJdisct.push_back("5");

  //eta binning, for multijet estimation
  std::vector<std::string> eta_bin; eta_bin.clear();
  eta_bin.push_back("");      
  eta_bin.push_back("_eta1");  
  eta_bin.push_back("_eta2");  
  eta_bin.push_back("_eta3");  
  eta_bin.push_back("_eta4");  
  eta_bin.push_back("_eta5");  
  eta_bin.push_back("_eta6");  
  eta_bin.push_back("_eta7");  
  eta_bin.push_back("_eta8");  
  eta_bin.push_back("_eta9");  
  eta_bin.push_back("_eta10"); 
  eta_bin.push_back("_eta11");

  //final fraction
  finalFractionAv.clear();

  //loop to fill the normalisation fraction vectors
  OriginalMJFracMET.clear();  OriginalMJFracMWT.clear();
  FinalMJFracMET.clear();     FinalMJFracMWT.clear();
  for(int k=0; k<(int)kine.size(); k++)
    for(int j=0;j<(int)eta_bin.size();j++)
      MultiJetNorm(config,kine[k],"a",eta_bin[j]);

  //get the shapes and the final MC+MJ distributions
  TH1D* htemp; data_eta_bin.clear(); mc_eta_bin.clear(); multijet_eta_bin.clear();
  for(int k=0; k<(int)kine.size(); k++){
    for(int j=0;j<(int)eta_bin.size();j++){
      htemp = MultiJetShape(config,kine[k],/*MJdisct[l]*/"1",eta_bin[j]);
      MultijetFinal(config,j,htemp,kine[k],eta_bin[j]);
    }
  }


  for(int j=0;j<(int)eta_bin.size();j++){
    std::cout<<FinalMJFracMET.at(j)<<"  "<<OriginalMJFracMET.at(j)<<"   ";
    std::cout<<status_met.at(j)<<"   "<<FinalMJFracMWT.at(j)<<"  ";
    std::cout<<OriginalMJFracMWT.at(j)<<"  "<< status_mwt.at(j)<<"  ";
    std::cout<<FinalMJFracMET.at(j)/FinalMJFracMWT.at(j)<<"  "<<0.5*(FinalMJFracMET.at(j)+FinalMJFracMWT.at(j))<<" ";
    std::cout<<(FinalMJFracMET.at(j)+FinalMJFracMWT.at(j))/(FinalMJFracMET.at(0)+FinalMJFracMWT.at(0))<<" "<<j<<" ";
    std::cout<<"   "<<data_eta_bin.at(j)<<"    "<<mc_eta_bin.at(j)<<"     "<<multijet_eta_bin.at(j)<<std::endl;
  }

  //dummy histograms to store the final numbers, ugly way!
  unsigned int nx = hData_eta->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[hData_eta->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hData_eta->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hData_eta->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hData_eta->GetXaxis()->GetNbins()]=hData_eta->GetXaxis()->GetBinUpEdge(hData_eta->GetXaxis()->GetNbins());

  TH1D* hDummy1 = new TH1D("data_final_numbers_eta"    ,     "",nx,xbins);
  TH1D* hDummy2 = new TH1D("mc_final_numbers_eta"      ,     "",nx,xbins);
  TH1D* hDummy3 = new TH1D("multijet_final_numbers_eta",     "",nx,xbins);
  TH1D* hDummy4 = new TH1D("fraction_final_numbers_eta",     "",nx,xbins);
  TH1D* hDummy5 = new TH1D("fractionerror_final_numbers_eta","",nx,xbins);

  for(int k=0;k<hDummy1->GetNbinsX()+1;k++){//k=0 inclsive eta bin
    hDummy1->SetBinError(k,0); hDummy1->SetBinContent(k,data_eta_bin.at(k));
    hDummy2->SetBinError(k,0); hDummy2->SetBinContent(k,mc_eta_bin.at(k));      
    hDummy3->SetBinError(k,0); hDummy3->SetBinContent(k,multijet_eta_bin.at(k));
    hDummy4->SetBinError(k,0); hDummy4->SetBinContent(k,finalFractionAv.at(k));         
    hDummy5->SetBinError(k,0); hDummy5->SetBinContent(k,finalFractionErrorAv.at(k));//fit error
  }

  fout->cd();
  hDummy1->Write();
  hDummy2->Write();
  hDummy3->Write();
  hDummy4->Write();
  hDummy5->Write();

  return;
}


TH1D* MultiJetMaker::MultiJetShape(Config config, std::string kine,
				   std::string MJdisct, std::string eta_bin){

  TH1D *hDataP, *hBGP, *hShapeP, *hShapeP_norm, *hDummy;
  hDataP = (TH1D*)fdata0 ->GetObjectUnchecked(("h_"+kine+"_"+MJdisct+eta_bin).c_str());
  hBGP   = (TH1D*)fmc0   ->GetObjectUnchecked(("h_"+kine+"_"+MJdisct+eta_bin).c_str()); hBGP->Scale(scale);
  hShapeP= (TH1D*)hDataP ->Clone(("shape_"+kine+"_antid0"+eta_bin).c_str());
  hShapeP->Add(hBGP,-1.);

  hShapeP_norm = (TH1D*)hShapeP->Clone(("shape_"+kine+"_norm_antid0"+eta_bin).c_str());
  hShapeP_norm->Scale(1./hShapeP_norm->Integral());

  fout->cd();
  hShapeP->Write(("shape_"+kine+"_antid0"+eta_bin).c_str());

  hDummy=(TH1D*)hShapeP->Clone("dummy");

  return hDummy;
}


void MultiJetMaker::MultiJetNorm(Config config, std::string kine, std::string MJdisct,std::string eta_bin){

  //signal
  hData       = (TH1D*)fdata0 ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str());
  hBG         = (TH1D*)fmc0   ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str()); hBG->Scale(scale);
  hDataAntid0 = (TH1D*)fdata0 ->GetObjectUnchecked(("h_"+kine+"_1"+eta_bin).c_str());            
  hBGAntid0   = (TH1D*)fmc0   ->GetObjectUnchecked(("h_"+kine+"_1"+eta_bin).c_str()); hBGAntid0->Scale(scale);
  hMultiSR    = (TH1D*)hDataAntid0->Clone(("multi_signal_"+kine+eta_bin).c_str());
  hMultiSR->Add(hBGAntid0,-1.);
  TCanvas *c = new TCanvas("canvas","canvas",550,450);
  c->cd();
  hMultiSR->Draw();
  if(kine=="eta") c->Print(("./multi_signal_"+kine+eta_bin+"_"+wzchannel+"_initial.pdf").c_str());
  delete c;

  //FR1
  hData1       = (TH1D*)fdata1 ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str());
  hBG1         = (TH1D*)fmc1   ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str()); hBG1->Scale(scale);  
  hData1Antid0 = (TH1D*)fdata1 ->GetObjectUnchecked(("h_"+kine+"_1"+eta_bin).c_str());                       
  hBG1Antid0   = (TH1D*)fmc1   ->GetObjectUnchecked(("h_"+kine+"_1"+eta_bin).c_str()); hBG1Antid0->Scale(scale);
  hMultiFR1    = (TH1D*)hData1Antid0->Clone(("multi_fr1_"+kine+eta_bin+"_initial").c_str());
  hMultiFR1->Add(hBG1Antid0,-1.);

  hData_nod0   = (TH1D*)fdata0 ->GetObjectUnchecked(("h_"+kine+"_0"+eta_bin).c_str());
  hBG_nod0     = (TH1D*)fmc0   ->GetObjectUnchecked(("h_"+kine+"_0"+eta_bin).c_str()); hBG_nod0->Scale(scale);
 
  fout->cd();
  hData       ->Write(("data_signal_d0_"+kine+eta_bin).c_str());
  hBG         ->Write(("bg_signal_d0_"+kine+eta_bin).c_str());
  hDataAntid0 ->Write(("data_signal_antid0_"+kine+eta_bin).c_str());
  hBGAntid0   ->Write(("bg_signal_antid0_"+kine+eta_bin).c_str());
  hMultiSR    ->Write(("multi_signal_"+kine+eta_bin+"_initial").c_str());

  hData1      ->Write(("data_fr1_d0_"+kine+eta_bin).c_str());
  hBG1        ->Write(("bg_fr1_d0_"+kine+eta_bin).c_str());
  hData1Antid0->Write(("data_fr1_antid0_"+kine+eta_bin).c_str());
  hBG1Antid0  ->Write(("bg_fr1_antid0_"+kine+eta_bin).c_str());
  hMultiFR1   ->Write(("multi_fr1_"+kine+eta_bin+"_initial").c_str());

  hData_nod0  ->Write();
  hBG_nod0    ->Write();

  //Template fit
  if(kine=="rmass" || kine=="met"){
    
    SingleNormalisation(hData1, hMultiFR1, hBG1, MultiFrac1, MCFrac1, MultiErr1, MCErr1,kine,eta_bin);
    //SingleNormalisation(hData2, hMultiFR2, hBG2, MultiFrac2, MCFrac2, MultiErr2, MCErr2);
        
    //if(status!=0){MCFrac1=1.; MultiFrac1=1.;}
    hData1_norm = (TH1D*)hData1->Clone(("data_fr1_d0_"+kine+"_norm"+eta_bin).c_str());
    hData1_norm->Scale(1./hData1_norm->Integral());
    
    hBG1_norm      = (TH1D*)hBG1->Clone(("bg_fr1_d0_"+kine+"_norm"+eta_bin).c_str());
    hBG1_norm->Scale(1./hBG1_norm->Integral());
    hBG1_norm_scaled = (TH1D*)hBG1_norm->Clone(("bg_fr1_d0_"+kine+"_norm_scaled"+eta_bin).c_str());
    hBG1_norm_scaled->Scale(MCFrac1);
    
    hMultiFR1_norm = (TH1D*)hMultiFR1->Clone(("multi_fr1_"+kine+"_norm"+eta_bin).c_str());
    hMultiFR1_norm->Scale(1./hMultiFR1_norm->Integral());
    hMultiFR1_norm_scaled = (TH1D*)hMultiFR1_norm->Clone(("multi_fr1_"+kine+"_norm_scaled"+eta_bin).c_str());
    hMultiFR1_norm_scaled->Scale(MultiFrac1);
    
    h_sum_fr1_norm_scaled     = (TH1D*)hBG1_norm_scaled->Clone(("sum_fr1_"+kine+"_norm_scaled"+eta_bin).c_str());
    h_sum_fr1_norm_scaled->Add(hMultiFR1_norm_scaled);
    
    //normalisation (template fit)
    //This is to get the actual T fraction
    double NormMC=hData1->Integral()/hBG1->Integral();
    //if(status!=0) NormMC=1.;
    hBG1_scaled    = (TH1D*)hBG1->Clone(("bg_fr1_d0_"+kine+"_scaled"+eta_bin).c_str());
    hBG1_scaled ->Scale(MCFrac1*NormMC);
    
    //This is to get the actual T fraction
    double NormMultijet = hData1->Integral()/hMultiFR1->Integral();
    //if(status!=0) NormMultijet=1.;
    hMultiFR1_scaled = (TH1D*)hMultiFR1->Clone(("multi_fr1_"+kine+"_scaled"+eta_bin).c_str());
    hMultiFR1_scaled->Scale(MultiFrac1*NormMultijet);
    
    h_sum_fr1_scaled          = (TH1D*)hBG1_scaled->Clone(("sum_fr1_"+kine+"_scaled"+eta_bin).c_str());
    h_sum_fr1_scaled->Add(hMultiFR1_scaled);

    h_sum_fr1 = (TH1D*)hBG1->Clone(("sum_fr1_"+kine+"_sum"+eta_bin).c_str());
    h_sum_fr1->Add(hMultiFR1);

    fout->cd();
    hData1_norm   ->Write();
    hBG1_norm     ->Write();
    hMultiFR1_norm->Write();
    hBG1_norm_scaled->Write();
    hMultiFR1_norm_scaled->Write();
    h_sum_fr1_norm_scaled->Write();//scaled

    hBG1_scaled->Write();
    hMultiFR1_scaled->Write();
    h_sum_fr1_scaled->Write();//scaled
    h_sum_fr1->Write();//not scaled

    if(kine=="rmass"){
      OriginalMJFracMWT.push_back(MultiFrac1); 
      FinalMJFracMWT.push_back(MultiFrac1*NormMultijet);
      FinalMJFracErrorMWT.push_back(MultiErr1*NormMultijet);
      status_mwt.push_back(status);
    }else if(kine=="met"){
      OriginalMJFracMET.push_back(MultiFrac1);
      FinalMJFracMET.push_back(MultiFrac1*NormMultijet);
      FinalMJFracErrorMET.push_back(MultiErr1*NormMultijet);
      status_met.push_back(status);
    }
    
    std::cout<<"****************        "<<kine<<"       *****************"<<std::endl;
    std::cout<<MultiFrac1<<"    FRACTIONS     "<<MCFrac1<<std::endl;
    std::cout<<MultiFrac1*NormMultijet<<"    FRACTIONS     "<<MCFrac1*NormMC<<std::endl; 
    std::cout<<NormMultijet<<"    NORMAL...     "<<NormMC<<std::endl;
  }

  return;
}


void MultiJetMaker::SingleNormalisation(TH1D *data, TH1D *mc0, TH1D *mc1,
					double &fracMu, double &fracMC,
					double &errMu,  double  &errMC,
					std::string kine, std::string eta_bin){

  //data->Rebin(1); mc0->Rebin(1); mc1->Rebin(1);
  //data->Scale(1./data->Integral()); mc0->Scale(1./data->Integral()); mc1->Scale(1./data->Integral());

  //find low statistics bins
  std::vector<int> exclude; exclude.clear();
  for(int k=1; k<mc0->GetXaxis()->GetNbins();k++)
    if(mc0->GetBinContent(k)<10) exclude.push_back(k);
  
  TObjArray *mc = new TObjArray(2);
  mc->Add(mc0);
  mc->Add(mc1);
  TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
  //fit->ReleaseRangeX();
  //fit->SetRangeX(1,29);
  fit->SetRangeX(1,data->GetXaxis()->GetNbins());
  // for(int k=0; k<(int)exclude.size(); k++)
  //   fit->ExcludeBin(exclude.at(k));

  status = fit->Fit();               // perform the fit
  std::cout << "fit status: " << status << std::endl;
  if (status == 0) {                    // check on fit status
    result = (TH1D*) fit->GetPlot();
    result_norm = (TH1D*)result->Clone(("fit_result_norm_"+kine+eta_bin).c_str());
    result_norm->Scale(1./result_norm->Integral());
    fout->cd();
    result->Write(("fit_result_"+kine+eta_bin).c_str());
    result_norm->Write();
    // data->Draw("Ep");
    // result->Draw("same");
    fit->GetResult(0,fracMu,errMu);
    fit->GetResult(1,fracMC,errMC);
  }else if (status==4){
    fit->GetResult(0,fracMu,errMu);
    fit->GetResult(1,fracMC,errMC);    
  }

  return;
}


void MultiJetMaker::MultijetFinal(Config config,int it, TH1D* h, std::string kine, std::string eta_bin){

  TH1D *hMultiSR_Final    = (TH1D*)h->Clone(("h_multijet_signal_final_"+kine+eta_bin).c_str());
  TH1D *hData_Final       = (TH1D*)fdata0_p ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str());
  TH1D *hMultijetMC_Final = (TH1D*)fmc0_p   ->GetObjectUnchecked(("h_"+kine+"_i"+eta_bin).c_str()); hMultijetMC_Final->Scale(scale);
  TH1D *hMC_Final         = (TH1D*)hMultijetMC_Final->Clone(("hMC_Final"+kine+eta_bin).c_str());

  if(kine=="eta"){
    data_eta_bin.push_back(hData_Final->Integral(0,hData_Final->GetNbinsX()+1));
    mc_eta_bin.push_back(hMultijetMC_Final->Integral(0,hMultijetMC_Final->GetNbinsX()+1));
  }

  //decide which fraction to use  
  double fmet=1.,fmwt=1.,fav=1.;
  double emet=0.,emwt=0.,eav=0;
  if(status_met.at(it)==0 || status_met.at(it)==4){fmet = FinalMJFracMET.at(it); emet=FinalMJFracErrorMET.at(it);}
  if(status_mwt.at(it)==0 || status_mwt.at(it)==4){fmwt = FinalMJFracMWT.at(it); emwt=FinalMJFracErrorMWT.at(it);}
  //fav = 0.5*(fmet + fmwt);
  if(status_met.at(it)==0 && status_mwt.at(it)==0){fav = 0.5*(fmet + fmwt); eav = 0.5*(emet + emwt);}
  if(status_met.at(it)==0 && status_mwt.at(it)!=0){fav = fmet;		    eav = emet;}		   
  if(status_met.at(it)!=0 && status_mwt.at(it)==0){fav = fmwt;              eav = emwt;}             

  double metOK=false, mwtOK=false;
  if(fmwt>2.5 && fmwt<4.5) mwtOK = true;
  if(fmet>2.5 && fmet<4.5) metOK = true;

  if(status_met.at(it)!=0 && status_mwt.at(it)!=0){
    if(mwtOK && metOK){fav = 0.5*(fmet + fmwt); fav = 0.5*(fmet + fmwt);}
    else if(mwtOK)    {fav = fmwt;	       fav = fmwt;}	      
    else if(metOK)    {fav = fmet;              fav = fmet;}             
    else if(!mwtOK && !metOK){fav = 0.5*(FinalMJFracMET.at(0)+FinalMJFracMWT.at(0)); eav = 0.5*(FinalMJFracErrorMET.at(0)+FinalMJFracErrorMWT.at(0));}
  }

  finalFractionAv.push_back(fav);
  finalFractionErrorAv.push_back(eav);

  //normalise the signal region MJ shape with the T (coming from the fit)
  std::cout<<hMultiSR->Integral(0,hMultiSR->GetNbinsX()+1)<<"    "<<hMultiSR_Final->Integral(0,hMultiSR_Final->GetNbinsX()+1)<<"   "<<  fav <<"    "<<kine<<"    debug"<<std::endl;
  hMultiSR_Final->Scale(fav);
  if(kine=="eta") multijet_eta_bin.push_back(hMultiSR_Final->Integral(0,hMultiSR_Final->GetNbinsX()+1));

  //Add MC+MJ
  hMultijetMC_Final->Add(hMultiSR_Final);

  //full distribution
  if(kine=="eta"){
    if(it==1) {hMC_eta  = (TH1D*)hMultijetMC_Final->Clone("eta_combined_final"); hMC_eta_multi  = (TH1D*)hMultiSR_Final->Clone("eta_combined_multi_final"); }
    if(it>1)  {hMC_eta->Add(hMultijetMC_Final); hMC_eta_multi->Add(hMultiSR_Final);}
    if(it==11){
      hData_eta = (TH1D*)fdata0_p ->GetObjectUnchecked(("h_"+kine+"_i").c_str());
      fout->cd();
      hData_eta->Write("eta_data_final");
      hMC_eta->Write();
      hMC_eta_multi->Write();      
      //hData_eta->Rebin(); hMC_eta->Rebin(); hMC_eta_multi->Rebin();
    }
  }

  fout->cd();
  hMultiSR_Final    ->Write(("mj_final_"+kine+eta_bin).c_str());
  hData_Final       ->Write(("data_final_"+kine+eta_bin).c_str()); 
  hMC_Final         ->Write(("mc_final_"+kine+eta_bin).c_str()); 
  hMultijetMC_Final ->Write(("mcmj_final_"+kine+eta_bin).c_str()); 

  hMultiSR_Final    = NULL;
  hData_Final       = NULL;
  hMC_Final         = NULL;
  hMultijetMC_Final = NULL;

  return;
}

void MultiJetMaker::finalize(){
  fout->Close();
  return;
}

#endif
