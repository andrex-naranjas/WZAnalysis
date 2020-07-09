//author: Andres Ramirez Morales <andres.ramirez.morales@cern.ch>
#ifndef MULTIJETRESULTS_CXX
#define MULTIJETRESULTS_CXX

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "MultiJetResults.h"


MultiJetResults::MultiJetResults()
{
}

MultiJetResults::~MultiJetResults(){}


void MultiJetResults::initialize(Config config){

  if(config.WZSelection=="zmumu"){
    wzchannel="z";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus";}
  
  std::string year="";
  if(config.DataYears=="2015+2016") {dirYear="2015p2016/"; year="15p16";}
  if(config.DataYears=="2017") {dirYear="2017/"; year="17";}
  if(config.DataYears=="2017") {dirYear="2018/"; year="18";}

  std::string dirInclusive="", total="";
  std::string systematic  = "nominal";

  fdata0 = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "data" + year + "_" + wzchannel + "_multijet_signal.root").c_str());					    
  fmc0   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "background_"   + wzchannel + "_nominal_nonewPU_multijet_signal" + total+"_25.root").c_str());    

  fout   = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "multijet_final.root").c_str(), "RECREATE");

  //Check if the files are open/exist
  if(!fdata0->IsOpen() || !fmc0->IsOpen() || !fout->IsOpen()){
    std::cout<<"Multijet results enabled, but not enough files, Ciao!"<<std::endl;
    exit(1);
  }

  setstyle();

  return;
}


void MultiJetResults::execute(Config config){

  //Plots
  Plots(config);
  fout->Close();
  //latex tables
  std::vector<std::string> eta_bin_name; eta_bin_name.clear();
  eta_bin_name=LatexBinName();
  MultiTable("multijet_wplus", eta_bin_name, "wplus");
  MultiTable("multijet_wminus", eta_bin_name,"wminus");
  TableCompare("multijet_compare", eta_bin_name);

  TableSys("multijet_wplus_sys", eta_bin_name, "wplus","d0sig");
  TableSys("multijet_wminus_sys", eta_bin_name,"wminus","d0sig");

  TableSys("multijet_wplus_sys", eta_bin_name, "wplus","xSec");
  TableSys("multijet_wminus_sys", eta_bin_name,"wminus","xSec");

  TableSys("multijet_wplus_sys", eta_bin_name, "wplus","stat");
  TableSys("multijet_wminus_sys", eta_bin_name,"wminus","stat");


  return;
}

void MultiJetResults::Plots(Config config){

  std::vector<std::string> kine; kine.clear();
  kine.push_back("rmass");
  kine.push_back("met");
  kine.push_back("d0sig");
  kine.push_back("pt");
  kine.push_back("eta");
  //kine.push_back("iso30"); 
  kine.push_back("eta_d");
  kine.push_back("phi");
  kine.push_back("phi_met");
  kine.push_back("delta_phi");
  kine.push_back("mu");


  //eta binning
  std::vector<std::string> eta_bin_label; eta_bin_label.clear();
  std::vector<std::string> eta_bin; eta_bin.clear();
  eta_bin.push_back("");         eta_bin_label.push_back(" Inclusive #||{#eta}");   
  eta_bin.push_back("_eta1");  	 eta_bin_label.push_back(" 0.00#leq#||{#eta}<0.21");
  // eta_bin.push_back("_eta2");  	 eta_bin_label.push_back(" 0.21#leq#||{#eta}<0.42");
  // eta_bin.push_back("_eta3");  	 eta_bin_label.push_back(" 0.42#leq#||{#eta}<0.63");
  // eta_bin.push_back("_eta4");  	 eta_bin_label.push_back(" 0.63#leq#||{#eta}<0.84");
  // eta_bin.push_back("_eta5");  	 eta_bin_label.push_back(" 0.84#leq#||{#eta}<1.05");
  // eta_bin.push_back("_eta6");  	 eta_bin_label.push_back(" 1.05#leq#||{#eta}<1.37");
  // eta_bin.push_back("_eta7");  	 eta_bin_label.push_back(" 1.37#leq#||{#eta}<1.52");
  // eta_bin.push_back("_eta8");  	 eta_bin_label.push_back(" 1.52#leq#||{#eta}<1.74");
  // eta_bin.push_back("_eta9");  	 eta_bin_label.push_back(" 1.74#leq#||{#eta}<1.95");
  // eta_bin.push_back("_eta10"); 	 eta_bin_label.push_back(" 1.95#leq#||{#eta}<2.18");
  // eta_bin.push_back("_eta11");	 eta_bin_label.push_back(" 2.18#leq#||{#eta}<2.40");

  std::vector<std::string> channel; channel.clear();
  channel.push_back("wplus"); channel.push_back("wminus"); 

  std::vector<std::string> systematic; systematic.clear();
  systematic.push_back("25");//=="Nominal") 
  systematic.push_back("30");//=="30")    
  systematic.push_back("20");//=="25")  
  systematic.push_back("xSecup");   
  systematic.push_back("xSecdown");

  for(int h=0; h<(int)systematic.size();h++){//sys loop
    for(int i=0;i<(int)channel.size();i++){//loop over channels
      fIN = new TFile((config.OutputFileDir+"Files/"+ dirYear + "Multijet/" + dirInclusive + "multijet_"+channel[i]+"_"+systematic[h]+".root").c_str());;
      //fIN = new TFile(("/data/morales/atlas/r21ControlPlots/Files/2017/Multijet/background_wminus_nominal_nonewPU_multijet_signal_original.root"));

      if(systematic[h]=="25"){      
	for(int j=0;j<(int)kine.size();j++){
	  for(int k=0;k<(int)eta_bin.size();k++){	
	    //signal
	    hMultiSR=NULL;
	    hData       = (TH1D*)fIN ->GetObjectUnchecked(("data_signal_d0_"+kine[j]+eta_bin[k]).c_str());
	    hBG         = (TH1D*)fIN ->GetObjectUnchecked(("bg_signal_d0_"+kine[j]+eta_bin[k]).c_str());
	    //hData       = (TH1D*)fdata0 ->GetObjectUnchecked("h_d0sig_i");//original plot
	    //hBG         = (TH1D*)fIN ->GetObjectUnchecked("h_d0sig_i");//orinal plot
	    hDataAntid0 = (TH1D*)fIN ->GetObjectUnchecked(("data_signal_antid0_"+kine[j]+eta_bin[k]).c_str());
	    hBGAntid0   = (TH1D*)fIN ->GetObjectUnchecked(("bg_signal_antid0_"+kine[j]+eta_bin[k]).c_str());
	    hMultiSR    = (TH1D*)fIN ->GetObjectUnchecked(("shape_"+kine[j]+"_antid0"+eta_bin[k]).c_str());//("multi_signal_"+kine[j]+eta_bin[k]+"_initial").c_str());
	        
	    //FR1
	    hData1       = (TH1D*)fIN ->GetObjectUnchecked(("data_fr1_d0_"+kine[j]+eta_bin[k]).c_str());
	    hBG1         = (TH1D*)fIN ->GetObjectUnchecked(("bg_fr1_d0_"+kine[j]+eta_bin[k]).c_str());
	    hData1Antid0 = (TH1D*)fIN ->GetObjectUnchecked(("data_fr1_antid0_"+kine[j]+eta_bin[k]).c_str());
	    hBG1Antid0   = (TH1D*)fIN ->GetObjectUnchecked(("bg_fr1_antid0_"+kine[j]+eta_bin[k]).c_str());
	    hMultiFR1    = (TH1D*)fIN ->GetObjectUnchecked(("multi_fr1_"+kine[j]+eta_bin[k]).c_str());
	    
	    hData_nod0   = (TH1D*)fdata0 ->GetObjectUnchecked(("h_"+kine[j]+"_0"+eta_bin[k]).c_str());
	    hBG_nod0     = (TH1D*)fmc0   ->GetObjectUnchecked(("h_"+kine[j]+"_0"+eta_bin[k]).c_str());
	    
	    //fit results
	    hData1_norm           = (TH1D*)fIN ->GetObjectUnchecked(("data_fr1_d0_"+kine[j]+"_norm"+eta_bin[k]).c_str());    
	    hBG1_norm             = (TH1D*)fIN ->GetObjectUnchecked(("bg_fr1_d0_"+kine[j]+"_norm"+eta_bin[k]).c_str());
	    hMultiFR1_norm        = (TH1D*)fIN ->GetObjectUnchecked(("multi_fr1_"+kine[j]+"_norm"+eta_bin[k]).c_str());
	    hBG1_norm_scaled      = (TH1D*)fIN ->GetObjectUnchecked(("bg_fr1_d0_"+kine[j]+"_norm_scaled"+eta_bin[k]).c_str());    
	    hMultiFR1_norm_scaled = (TH1D*)fIN ->GetObjectUnchecked(("multi_fr1_"+kine[j]+"_norm_scaled"+eta_bin[k]).c_str());      
	    h_sum_fr1_norm_scaled = (TH1D*)fIN ->GetObjectUnchecked(("sum_fr1_"+kine[j]+"_norm_scaled"+eta_bin[k]).c_str());    
	    
	    hBG1_scaled           = (TH1D*)fIN ->GetObjectUnchecked(("bg_fr1_d0_"+kine[j]+"_scaled"+eta_bin[k]).c_str());
	    hMultiFR1_scaled      = (TH1D*)fIN ->GetObjectUnchecked(("multi_fr1_"+kine[j]+"_scaled"+eta_bin[k]).c_str());      
	    h_sum_fr1_scaled      = (TH1D*)fIN ->GetObjectUnchecked(("sum_fr1_"+kine[j]+"_scaled"+eta_bin[k]).c_str());  
	    h_sum_fr1             = (TH1D*)fIN ->GetObjectUnchecked(("sum_fr1_"+kine[j]+"_sum"+eta_bin[k]).c_str());
	    
	    //shape
	    // hShapeP      = (TH1D*)fIN ->GetObjectUnchecked(("shape_"+kine[j]+"_antid0"+eta_bin[k]).c_str());
	    // hShapeP_norm = (TH1D*)fIN ->GetObjectUnchecked(("shape_"+kine[j]+"_norm_antid0"+eta_bin[k]).c_str());
	    
	    //Plots
	    //control plots/shape plots
	    //hBG->Scale(1.04);  hBG_nod0->Scale(1.02);
	    ComparePlot(config,kine[j],"DataMC","Data/Pred.", hData, hBG,NULL,NULL, "FullSelection",eta_bin[k],eta_bin_label[k], 0, 1e7,false,channel[i]);
	    ComparePlot(config,kine[j],"DataMC","Data/Pred.", hData_nod0,hBG_nod0,NULL,NULL, "FullSelection_nod0",eta_bin[k], eta_bin_label[k], 1e3, 1e9,true,channel[i]);

	    ComparePlot(config,kine[j],"DataMC","Data/Pred.",hDataAntid0,hBGAntid0,NULL,NULL,"Signal_antid0_multi",eta_bin[k],eta_bin_label[k],0,1e6,false,channel[i]);
	    
	    double ymax=2.5e5;
	    if(kine[j]=="eta") ymax=1.25e6;
	    else if(kine[j]=="phi" || kine[j]=="phi_met") ymax=7e4;
	    if(eta_bin[i]!="") ymax=ymax/1.;
	    SinglePlot(config,kine[j],"DataMCMulti",hDataAntid0,hBGAntid0,hMultiSR,NULL,"Signal_antid0_multi",eta_bin[k],eta_bin_label[k],0.1,ymax,false,channel[i]);
	    SinglePlot(config,kine[j],"DataMCMulti",hData, hBG,hMultiSR,NULL,"Signal_d0_multi",eta_bin[k],eta_bin_label[k],0,1e7,false,channel[i]);
	    
	    //normalisation (template fit)
	    if(kine[j]=="rmass" || kine[j]=="met"){
	      
	      result_norm = (TH1D*)fIN ->GetObjectUnchecked(("fit_result_norm_"+kine[j]+eta_bin[k]).c_str());
	      result      = (TH1D*)fIN ->GetObjectUnchecked(("fit_result_"+kine[j]+eta_bin[k]).c_str());
	      
	      double ymax=0.;
	      if(kine[j]=="rmass") ymax=1.e7;
	      if(kine[j]=="met")   ymax=1.75e7;
	      //if(config.WZSelection=="wplus")    ymax=2*ymax;
	      if(eta_bin[i]=="_eta1") ymax=ymax/1.;
	      //hBG->Scale(1.04);	      hBG1->Scale(1.04);
	      SinglePlot(config,kine[j],"DataMC",hData,  hBG, NULL,NULL,"Signal_d0",eta_bin[k],eta_bin_label[k],0,ymax,false,channel[i]);
	      SinglePlot(config,kine[j],"DataMC",hData1, hBG1,NULL,NULL,"FR1_d0",eta_bin[k],eta_bin_label[k],0,ymax,false,channel[i]);
	      
	      //hBG->Scale(0.96);	      hBG1->Scale(0.96);
	      if(kine[j]=="rmass") ymax=1.5e7;
	      if(kine[j]=="met")   ymax=5e8;
	      if(config.WZSelection=="wplus") ymax=2*ymax;
	      if(eta_bin[i]=="_eta1") ymax=ymax/1.;
	      SinglePlot(config,kine[j],"DataMCMulti",hDataAntid0, hBGAntid0,hMultiSR,NULL,"Signal_antid0",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      SinglePlot(config,kine[j],"DataMCMulti",hData1Antid0, hBG1Antid0,hMultiFR1,NULL,"FR1_antid0",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);    
	      if(kine[j]=="rmass") ymax=1.5e9;
	      if(kine[j]=="met")   ymax=1e11;
	      if(config.WZSelection=="wplus") ymax=2*ymax;
	      if(eta_bin[i]!="") ymax=ymax/1.;
	      SinglePlot(config,kine[j],"DataMCMulti",hData1, hBG1,hMultiFR1,NULL,"FR1_d0_multi",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      
	      double ynormmax=0.;
	      if(kine[j]=="rmass") ynormmax=1.5e4;
	      if(kine[j]=="met")   ynormmax=5e4;
	      //results
	      
	      // if(result==NULL) std::cin.get();
	      if(result!=NULL && result_norm!=NULL){
		SinglePlot(config,kine[j],"DataComp",hData1, result, hBG1, hMultiFR1,"FR1_antid0",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
		SinglePlot(config,kine[j],"DataCompNorm",hData1_norm,result_norm,hBG1_norm,hMultiFR1_norm,"FR1_antid0_norm",eta_bin[k],eta_bin_label[k],0,ynormmax,true,channel[i]);
		SinglePlot(config,kine[j],"DataCompNorm",hData1_norm,result_norm,hBG1_norm_scaled,hMultiFR1_norm_scaled,"FR1_antid0_norm_scaled",eta_bin[k],eta_bin_label[k],0,ynormmax,true,channel[i]);
	      }
	    
	      if(result!=NULL)
		ComparePlot(config,kine[j],"DataFit","Data/Fit",hData1, result,NULL,NULL,"FR1_antid0",eta_bin[k],eta_bin_label[k],0, ymax,true,channel[i]);
	      
	      ComparePlot(config,kine[j],"DataFitSumNorm","Data/Fit",hData1_norm,h_sum_fr1_norm_scaled,NULL,NULL,"FR1_antid0_norm_sum", eta_bin[k],eta_bin_label[k], 0,ynormmax,true,channel[i]);

	      if(h_sum_fr1!=NULL && hBG1!=NULL && hMultiFR1!=NULL) {h_sum_fr1->Scale(1.06); hBG1->Scale(1.06); hMultiFR1->Scale(1.06);}
	      ComparePlot(config,kine[j],"DataSumNoScaled","Data/Pred.",hData1, h_sum_fr1,hBG1,hMultiFR1,"FR1_antid0_sum", eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      ComparePlot(config,kine[j],"DataSum","Data/Pred.",hData1, h_sum_fr1_scaled,hBG1_scaled,hMultiFR1_scaled,"FR1_antid0_scaled_sum", eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      ComparePlot(config,kine[j],"DataSumOnly","Data/Pred.",hData1, h_sum_fr1_scaled,NULL,NULL,"FR1_antid0_sum",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      
	      //final template fit plot
	      if(result!=NULL)
		ComparePlot(config,kine[j],"DataFitFinal","Data/Fit",hData1,result,hBG1_scaled,hMultiFR1_scaled,"FR1_antid0",eta_bin[k],eta_bin_label[k],0,ymax,true,channel[i]);
	      
	    }//if met && mwt
	    
	  }//eta bin loop
	  
	  //final plots      
	  hMultiSR_Final    = (TH1D*)fIN ->GetObjectUnchecked(("mj_final_"+kine[j]+eta_bin[0]).c_str());
	  hData_Final       = (TH1D*)fIN ->GetObjectUnchecked(("data_final_"+kine[j]+eta_bin[0]).c_str()); 
	  hMC_Final         = (TH1D*)fIN ->GetObjectUnchecked(("mc_final_"+kine[j]+eta_bin[0]).c_str()); 
	  hMultijetMC_Final = (TH1D*)fIN ->GetObjectUnchecked(("mcmj_final_"+kine[j]+eta_bin[0]).c_str()); 
	  
	  //write to the multijet final file
	  fout->cd();
	  hMultiSR_Final   ->Write((kine[j]+"_"+channel[i]).c_str());
	  
	  double ymax=1e7;
	  if(kine[j]=="eta") ymax=5.5e7;
	  if(config.WZSelection=="wplus") ymax=2*ymax;
	  if(eta_bin[i]!="") ymax=ymax/1.;
	  ComparePlot(config,kine[j],"Final","Data/Pred",hData_Final,hMultijetMC_Final,hMultiSR_Final,NULL,"SR_d0",eta_bin[0],eta_bin_label[0],0,ymax,false,channel[i]);
		
	}//kine loop
      
	hMC_eta        = (TH1D*)fIN ->GetObjectUnchecked(("eta_combined_final"));
	hMC_eta_multi  = (TH1D*)fIN ->GetObjectUnchecked(("eta_combined_multi_final")); 
	hData_eta      = (TH1D*)fIN ->GetObjectUnchecked(("eta_data_final"));
	
	//write to the multijet final file
	fout->cd();
	hMC_eta_multi->Write(("eta_bin_"+channel[i]).c_str());
      
	double ymax=1e7;
	if(config.WZSelection=="wplus") ymax=1.2*ymax;
	ymax=1e8;
	ComparePlot(config,"eta","FinalEtaCombined","Data/Pred.",hData_eta,hMC_eta,hMC_eta_multi,NULL,"SR_d0",eta_bin[0],eta_bin_label[0],0,ymax,false,channel[i]);
      
      }//if sys nominal


      //for tables
      data_numbers     = (TH1D*)fIN ->GetObjectUnchecked(("data_final_numbers_eta"    ));
      mc_numbers       = (TH1D*)fIN ->GetObjectUnchecked(("mc_final_numbers_eta"      ));
      multijet_numbers = (TH1D*)fIN ->GetObjectUnchecked(("multijet_final_numbers_eta"));
      fraction_numbers = (TH1D*)fIN ->GetObjectUnchecked(("fraction_final_numbers_eta"));
      fractionerror_numbers = (TH1D*)fIN ->GetObjectUnchecked(("fractionerror_final_numbers_eta"));
      
      if(systematic[h]=="25"){
	for(int k=0;k<data_numbers->GetNbinsX()+1;k++){     
	  if(channel[i]=="wminus"){
	    data_eta_bin_wminus.push_back(data_numbers->GetBinContent(k));
	    mc_eta_bin_wminus.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wminus.push_back(multijet_numbers->GetBinContent(k));
	    finalFractions_wminus.push_back(fraction_numbers->GetBinContent(k));
	    finalErrorFractions_wminus.push_back(fractionerror_numbers->GetBinContent(k));
	  }
	  if(channel[i]=="wplus"){
	    data_eta_bin_wplus.push_back(data_numbers->GetBinContent(k));
	    mc_eta_bin_wplus.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wplus.push_back(multijet_numbers->GetBinContent(k));
	    finalFractions_wplus.push_back(fraction_numbers->GetBinContent(k));
	    finalErrorFractions_wplus.push_back(fractionerror_numbers->GetBinContent(k));
	  }      
	}
      }

      if(systematic[h]=="20"){
	for(int k=0;k<data_numbers->GetNbinsX()+1;k++){     
	  if(channel[i]=="wminus"){
	    mc_eta_bin_wminus_d0down.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wminus_d0down.push_back(multijet_numbers->GetBinContent(k));
	  }
	  if(channel[i]=="wplus"){
	    mc_eta_bin_wplus_d0down.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wplus_d0down.push_back(multijet_numbers->GetBinContent(k));
	  }      
	}
      }

      if(systematic[h]=="30"){
	for(int k=0;k<data_numbers->GetNbinsX()+1;k++){     
	  if(channel[i]=="wminus"){	 
	    mc_eta_bin_wminus_d0up.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wminus_d0up.push_back(multijet_numbers->GetBinContent(k));
	  }
	  if(channel[i]=="wplus"){	 
	    mc_eta_bin_wplus_d0up.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wplus_d0up.push_back(multijet_numbers->GetBinContent(k));
	  }      
	}
      }

      if(systematic[h]=="xSecup"){
	for(int k=0;k<data_numbers->GetNbinsX()+1;k++){     
	  if(channel[i]=="wminus"){	 
	    mc_eta_bin_wminus_xSecup.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wminus_xSecup.push_back(multijet_numbers->GetBinContent(k));
	  }
	  if(channel[i]=="wplus"){	 
	    mc_eta_bin_wplus_xSecup.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wplus_xSecup.push_back(multijet_numbers->GetBinContent(k));
	  }      
	}
      }

      if(systematic[h]=="xSecdown"){
	for(int k=0;k<data_numbers->GetNbinsX()+1;k++){     
	  if(channel[i]=="wminus"){	 
	    mc_eta_bin_wminus_xSecdown.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wminus_xSecdown.push_back(multijet_numbers->GetBinContent(k));
	  }
	  if(channel[i]=="wplus"){	 
	    mc_eta_bin_wplus_xSecdown.push_back(mc_numbers->GetBinContent(k));
	    multijet_eta_bin_wplus_xSecdown.push_back(multijet_numbers->GetBinContent(k));
	  }      
	}
      }
      
    }//channel for loop
  }//sys loop
  
  return;
}


void MultiJetResults::MultiTable(std::string fileName, std::vector<std::string> eta_bin_name,
				 std::string channel){

  std::vector<double> finalFractionAv;
  std::vector<double> data_eta_bin;
  std::vector<double> mc_eta_bin;
  std::vector<double> multijet_eta_bin;
  std::string labelboson;
  
  if(channel=="wminus"){
    labelboson="$W^{-}\\rightarrow\\mu^{-}\\nu$";
    multijet_eta_bin = multijet_eta_bin_wminus;
    mc_eta_bin	     = mc_eta_bin_wminus;
    data_eta_bin     = data_eta_bin_wminus;
    finalFractionAv  = finalFractions_wminus;
  }else if (channel=="wplus"){
    labelboson="$W^{+}\\rightarrow\\mu^{+}\\nu$";
    multijet_eta_bin =   multijet_eta_bin_wplus;
    mc_eta_bin	     =   mc_eta_bin_wplus;
    data_eta_bin     =   data_eta_bin_wplus;
    finalFractionAv  =   finalFractions_wplus;
  }
  
  std::vector<std::string> header; header.clear();
  header.push_back("$T$");
  //header.push_back("$T(m_{T})$");
  //header.push_back("$T(E_{T}^{miss})/T(m_{T})$"); header.push_back("$(T(E_{T}^{miss})+T(m_{T}))/2$");
  //header.push_back("Average ratio");
  header.push_back("$N_{Data}$");
  header.push_back("$N_{MC}$");
  header.push_back("$N_{Multijet}$");
  header.push_back("$N_{Multijet}/N_{Data}$");

  int nHeader = (int)header.size(); int nEtaBins = (int)data_eta_bin.size();
  std::ofstream multiLatex(("./"+fileName+"_T.tex").c_str());

  multiLatex << "\\begin{table}[tb!]" << std::endl;
  multiLatex << "\\begin{center}" << std::endl;
  multiLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++){
    if(k==0)multiLatex<<" C{1.25cm} "; 
    else if(k==1 || k==2) multiLatex<<" C{2cm} "; 
    else multiLatex<<" c ";}
  multiLatex<<" }\\hline \\hline " << std::endl;

  //header row
  multiLatex<< "Eta bin & ";
  for(int k=0; k < nHeader; k++){
    multiLatex<<std::setw(4)<<std::right;
    if(k< nHeader-1) multiLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) multiLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nEtaBins; j++){
    if(multijet_eta_bin.at(j)<0) multijet_eta_bin.at(j)=0;
    //multiLatex<<std::setw(6)<<std::right;
    multiLatex<<"  "<<eta_bin_name.at(j)    <<" &";
    multiLatex<<"  \\num{"<<finalFractionAv.at(j) <<"} &";  /*std::fixed<<std::setprecision(2)<<*/
    multiLatex<<"  \\num{"<<data_eta_bin.at(j)    <<"} &";
    multiLatex<<"  \\num{"<<mc_eta_bin.at(j)      <<"} &";
    multiLatex<<"  \\num{"<<multijet_eta_bin.at(j)<<"} &";
    multiLatex<<"  \\num{"<<100.*(multijet_eta_bin.at(j)/data_eta_bin.at(j))<<"}  \\\\ \\hline"<<std::endl;
  }
  
  multiLatex << "\\hline" << std::endl;
  multiLatex << "\\end{tabular}" << std::endl;
  multiLatex << "\\end{center}"  << std::endl;
  multiLatex << "\\caption{This table displays the results of the multijet calculation in each muon pseudorapidity bin for the "<<labelboson<<" selection. The first column indicates the muon pseudorapidity bin. The second column contains the $T$ values corresponding to the multijet normalisation. The third, fourth, and fifth columns are the number of events found in each bin for data, Monte Carlo simulation, and multijet respectively. The last column is the multijet fraction of events w.r.t. the data number of events, given in percentage.}" << std::endl;  
  multiLatex << "\\label{tab:T"<<channel<<"}" << std::endl;  
  multiLatex << "\\end{table}"<< std::endl;
  multiLatex.close();

  return;
}


void MultiJetResults::TableSys(std::string fileName, std::vector<std::string> eta_bin_name,
			       std::string channel, std::string systematic){

  std::vector<double> data;         
  std::vector<double> nominal_multi;
  std::vector<double> up_multi;
  std::vector<double> down_multi;
  std::string labelboson;
  std::string labelsys;

  if(channel=="wminus"){
    labelboson="$W^{-}\\rightarrow\\mu^{-}\\nu$";
    data          = data_eta_bin_wminus;
    nominal_multi = multijet_eta_bin_wminus;	  
  }else if(channel=="wplus"){
    labelboson="$W^{+}\\rightarrow\\mu^{+}\\nu$";
    data          = data_eta_bin_wplus;
    nominal_multi = multijet_eta_bin_wplus;
  }
  
  if(channel=="wminus" && systematic=="d0sig"){
    labelsys=" $d_{0}$ significance weight variation by $\\pm 0.5$ standard deviation ";
    up_multi    = multijet_eta_bin_wminus_d0up;  
    down_multi  = multijet_eta_bin_wminus_d0down;    
  }else if(channel=="wplus" && systematic=="d0sig"){
    labelsys=" $d_{0}$ significance weight variation by $\\pm 0.5$ standard deviation ";
    up_multi    = multijet_eta_bin_wplus_d0up;   
    down_multi  = multijet_eta_bin_wplus_d0down;     
  }else if(channel=="wminus" && systematic=="xSec"){
    labelsys=" Monte Carlo cross section variation of $\\pm 0.5\\%$";
    up_multi    = multijet_eta_bin_wminus_xSecup;   
    down_multi  = multijet_eta_bin_wminus_xSecdown; 
  }else if(channel=="wplus" && systematic=="xSec"){
    labelsys=" Monte Carlo cross section variation of $\\pm 0.5\\%$";
    up_multi    = multijet_eta_bin_wplus_xSecup;   
    down_multi  = multijet_eta_bin_wplus_xSecdown; 
  }else if(channel=="wminus" && systematic=="stat"){
    labelsys=" Statistical uncertainty ";
    up_multi    = finalErrorFractions_wminus;   
    down_multi  = finalErrorFractions_wminus; 
  }else if(channel=="wplus" && systematic=="stat"){
    labelsys=" Statistical uncertainty";
    up_multi    = finalErrorFractions_wplus;   
    down_multi  = finalErrorFractions_wplus; 
  }
 
  std::vector<std::string> header; header.clear();
  //header.push_back("$N_{Data}$");
  header.push_back("$N_{Multijet}^{nominal}$");
  header.push_back("$N_{Multijet}^{up}$");
  header.push_back("$N_{Multijet}^{down}$");
  header.push_back("$\\delta^{up}$");
  header.push_back("$\\delta^{down}$");

  int nHeader = (int)header.size(); int nEtaBins = (int)data.size();
  std::ofstream multiLatex(("./"+fileName+"_"+systematic+".tex").c_str());

  multiLatex << "\\begin{table}[tb!]" << std::endl;
  multiLatex << "\\begin{center}" << std::endl;
  multiLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) {
    if(k==4 || k==5) multiLatex<<" C{1.5cm} ";
    else multiLatex<<" C{2.25cm} ";
  }
  multiLatex<<" }\\hline \\hline " << std::endl;

  //header row
  multiLatex<< "Eta bin & ";
  for(int k=0; k < nHeader; k++){
    multiLatex<<std::setw(4)<<std::right;
    if(k< nHeader-1) multiLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) multiLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  for(int j=0; j < nEtaBins; j++){
    //multiLatex<<std::setw(6)<<std::right;
    multiLatex<<"  "<<eta_bin_name.at(j)  <<" &";
    //multiLatex<<"  \\num{"<<data.at(j)          <<"} &";  /*std::fixed<<std::setprecision(2)<<*/
    if(nominal_multi.at(j)<0 ){
      nominal_multi.at(j)=0;
      up_multi.at(j)= 0;
      down_multi.at(j) = 0;
    }else if(up_multi.at(j)<0 || down_multi.at(j)<0 ){ 
      up_multi.at(j)= 3 *nominal_multi.at(j);
      down_multi.at(j) = 3 *nominal_multi.at(j);
    }

    if(systematic!="stat"){
      multiLatex<<"  \\num{"<<nominal_multi.at(j) <<"} &";
      multiLatex<<"  \\num{"<<up_multi.at(j)   <<"} &";
      multiLatex<<"  \\num{"<<down_multi.at(j) <<"} &";
    }else{
      multiLatex<<"  \\num{"<<nominal_multi.at(j) <<"} &";
      multiLatex<<"  \\num{"<<up_multi.at(j)*nominal_multi.at(j)+nominal_multi.at(j)  <<"} &";
      multiLatex<<"  \\num{"<<down_multi.at(j)*nominal_multi.at(j)+nominal_multi.at(j) <<"} &";
    }

    if(nominal_multi.at(j)==0){nominal_multi.at(j)=1; up_multi.at(j)= 3; down_multi.at(j) = 3;}
    if(systematic!="stat"){
      multiLatex<<"  \\num{"<<std::fabs(100.*(nominal_multi.at(j)-up_multi.at(j))/nominal_multi.at(j)   )<<"} &";                      
      multiLatex<<"  \\num{"<<std::fabs(100.*(nominal_multi.at(j)-down_multi.at(j))/nominal_multi.at(j) )<<"} \\\\ \\hline"<<std::endl;
    }else{
      multiLatex<<"  \\num{"<<std::fabs(100.*(nominal_multi.at(j)-(up_multi.at(j)*nominal_multi.at(j)+nominal_multi.at(j)) ) /nominal_multi.at(j)  ) <<"} &";
      multiLatex<<"  \\num{"<<std::fabs(100.*(nominal_multi.at(j)-(down_multi.at(j)*nominal_multi.at(j)+nominal_multi.at(j))) /nominal_multi.at(j) )<<"} \\\\ \\hline"<<std::endl;
    }
  }

  
  multiLatex << " \\hline" << std::endl;
  multiLatex << "\\end{tabular}" << std::endl;
  multiLatex << "\\end{center}"  << std::endl;
  multiLatex << "\\caption{This table displays the systematics uncertainties arising from the "<<labelsys<<" in the "<<labelboson<<" selection. The first column contains the muon pseudorapidity bins; the second column is the nominal number of multijet background events; the third and fourth columns are the up and down varied number of multijet events; the fifth and sixth columns show the up and down variation computed with \\cref{eq:binsys}, in percentage.}" << std::endl;  
  multiLatex << "\\label{tab:"<<systematic<<channel<<"}" << std::endl;  
  multiLatex << "\\end{table}" << std::endl;
  multiLatex.close();

  return;
}


void MultiJetResults::TableCompare(std::string fileName, std::vector<std::string> eta_bin_name){
  
  std::vector<std::string> header; header.clear();
  header.push_back("$N_{Multijet}(W^{-})$");
  //header.push_back("$N_{Multijet}/N_{Data}$");
  header.push_back("$N_{Multijet}(W^{+})$");
  //header.push_back("$N_{Multijet}/N_{Data}$");
  header.push_back("$\\Delta_{Multijet}$");


  int nHeader = (int)header.size(); int nEtaBins = (int)data_eta_bin_wplus.size();
  std::ofstream multiLatex(("./"+fileName+".tex").c_str());

  std::cout<<"dalliiii "<<nHeader<<"   "<<nEtaBins<<std::endl;

  multiLatex << "\\begin{table}[tb!]" << std::endl;
  multiLatex << "\\begin{center}" << std::endl;
  multiLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++) multiLatex<<" C{3.5cm} "; 
  multiLatex<<" }\\hline \\hline " << std::endl;

  //header row
  multiLatex<< "Eta bin & ";
  for(int k=0; k < nHeader; k++){
    multiLatex<<std::setw(4)<<std::right;
    if(k< nHeader-1) multiLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) multiLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }
 
  for(int j=0; j < nEtaBins; j++){
    //multiLatex<<std::setw(6)<<std::right;

    if(multijet_eta_bin_wminus.at(j)<0) multijet_eta_bin_wminus.at(j)=0;
    if(multijet_eta_bin_wplus.at(j)<0)  multijet_eta_bin_wplus.at(j)=0;

    multiLatex<<"  "<<eta_bin_name.at(j)<<" &";
    multiLatex<<"  \\num{"<<multijet_eta_bin_wminus.at(j) <<"} &";  /*std::fixed<<std::setprecision(2)<<*/
    //multiLatex<<"  \\num{"<<100.*(multijet_eta_bin_wminus.at(j)/data_eta_bin_wminus.at(j)) <<"} &";
    multiLatex<<"  \\num{"<<multijet_eta_bin_wplus.at(j)      <<"} &";
    //multiLatex<<"  \\num{"<<100.*(multijet_eta_bin_wplus.at(j)/data_eta_bin_wplus.at(j))<<"} &";
    multiLatex<<"  \\num[explicit-sign=+]{"<<100.*(1.-multijet_eta_bin_wplus.at(j)/multijet_eta_bin_wminus.at(j)) <<"} \\\\ \\hline"<<std::endl;
  }
  
  multiLatex << " \\hline" << std::endl;
  multiLatex << "\\end{tabular}" << std::endl;
  multiLatex << "\\end{center}" << std::endl;
  multiLatex << "\\caption{This table compares the yields, in both positive and negative $W$ boson selection channels, coming from the multijet calculation. The first column indicates the muon pseudorapidity bin; the second and third columns correspond to the multijet events coming from the $W^{-}\\rightarrow\\mu^{-}\\nu$ and $W^{-}\\rightarrow\\mu^{-}\\nu$ selections, respectively; the fourth column is the percentage difference between the positive channel w.r.t. the negative channel multijet events, defined as $\\Delta_{Multijet}=(N_{Multijet}(W^{+})-N_{Multijet}(W^{-}))/N_{Multijet}(W^{+})$.}"<< std::endl;  
  multiLatex << "\\label{tab:mjcomp}"<< std::endl;  
  multiLatex << "\\end{table}" << std::endl;
  multiLatex.close();

  return;
}


void MultiJetResults::finalize(){

  return;
}


void MultiJetResults::ComparePlot(Config config, std::string kine,
				  std::string option, std::string ratiolabel, 
				  TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4, 
				  std::string cuts, std::string eta_bin, std::string label,
				  double ylow, double yhigh, bool ylog, std::string channel){
  
  std::string xlabel,ylabel; double xlow=0., xhigh=0.;
  TString leglabel="";

  if(kine=="rmass") xlabel="m_{T} [GeV]";
  if(kine=="d0sig") xlabel="d0sig";
  if(kine=="pt")    xlabel="p_{T} [GeV]";
  if(kine=="eta")   xlabel="#||{#eta^{#mu}}";
  if(kine=="met")   xlabel="E_{T}^{miss} [GeV]";
  if(kine=="iso30") xlabel= "iso30";

  if(channel=="zmumu")  leglabel="Z#rightarrow#mu#mu";
  if(channel=="wplus")  leglabel="W^{+}#rightarrow#mu^{+}#nu";
  if(channel=="wminus") leglabel="W^{-}#rightarrow#mu^{-}#nu";

  //  xlabel="d0sig"; leglabel="W^{-}#rightarrow#mu^{-}#nu"; for original plot
  ylabel="Events";
  if(option=="DataCompNorm" || option=="DataFitSumNorm") ylabel="Arbitrary Units";

  if(ylog) ylow = 1;
  if(ylog && (option=="DataCompNorm" || option=="DataFitSumNorm") ) ylow = 1e-10;
  if(kine=="d0sig" && ylog) ylow = 1e2;

  //divide bin eta binwidth to smooth the plot
  if(kine=="eta"){
    for(int k=1;k<hist1->GetXaxis()->GetNbins()+1;k++){
      hist1->SetBinContent(k,(hist1   ->GetBinContent(k)/hist1   ->GetBinWidth(k)));
      hist2->SetBinContent(k,(hist2   ->GetBinContent(k)/hist2   ->GetBinWidth(k)));
      if(hist3!=NULL)
	hist3->SetBinContent(k,(hist3   ->GetBinContent(k)/hist3   ->GetBinWidth(k)));
      if(hist4!=NULL)
	hist4->SetBinContent(k,(hist4   ->GetBinContent(k)/hist4   ->GetBinWidth(k)));
    }
  }

  plotAxisLine(hist1,kBlack,kGreen,20,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  plotAxisLine(hist2,kBlue,kGreen,28,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  if(hist3!=NULL){
    if(hist4==NULL) plotAxisLine(hist3,kRed,kGreen,21,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
    if(hist4!=NULL) plotAxisLine(hist3,kViolet+2,kGreen,21,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  }
  
  if(hist4!=NULL) 
    plotAxisLine(hist4,kRed,kGreen,20,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);

  TString lumi="44.3 fb^{-1}";
  TCanvas *c = new TCanvas("canvas","canvas",550,600);
  TLegend* leg = new TLegend(0.175,0.8,0.375,0.65,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  TString leg1, leg2, leg3, leg4;
  if(option=="DataMC" || option=="DataFit" || option=="DataFitSum" || option=="DataFitSumNorm") leg1="Data";
  if(option=="DataMC") leg2="MC Signal+BG";
  if(option=="DataFit")leg2="Fit";
  if(option=="DataFitSumNorm")leg2="Norm+weighted MC+MJ";
  if(option=="DataFitSum")leg2="Weighted MC+MJ";
  if(option=="FitSum"){leg1="Fit"; leg2="MC Signal+BG";}
  if(option=="DataSumNoScaled" || option=="DataSum")
    {leg1="Data"; leg2="MC+Multijet"; leg3="MC Signal+BG"; leg4="Multijet";}
  if(option=="DataSumOnly") {leg1="Data"; leg2="Weighted MC+Multijet";}
  if(option=="DataFitFinal")
    {leg1="Data"; leg2="Fit"; leg3="MC Signal+BG"; leg4="Multijet";}

  if(option=="Final" || option=="FinalEtaCombined")
    {leg1="Data"; leg2="MC Signal+BG+MJ"; leg3="Multijet";}

  leg->AddEntry(hist1, "#bf{#scale[0.85]{"+leg1+"}}","P");  
  leg->AddEntry(hist2, "#bf{#scale[0.85]{"+leg2+"}}","PL");
  if(hist3!=NULL && option!="Final")
    leg->AddEntry(hist3, "#bf{#scale[0.85]{"+leg3+"}}","PL");
  if(hist3!=NULL && option=="Final")
    leg->AddEntry(hist3, "#bf{#scale[0.85]{"+leg3+"}}","f");
  if(hist4!=NULL)
    leg->AddEntry(hist4, "#bf{#scale[0.85]{"+leg4+"}}","PL");

  TPaveText *box;
  box = new TPaveText(0.2,0.815,0.375,0.915,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi);

  TPaveText *box2;
  box2 = new TPaveText(0.745,0.75,0.795,0.875,"NDC");
  box2->SetBorderSize(0); box2->SetTextSize(0.035); box2->SetFillColor(0);
  box2->AddText(leglabel);
  box2->AddText((cuts+" region").c_str());
  box2->AddText(label.c_str());

  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0.2,1,1);
  if(ylog) pad->SetLogy();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.05);
  pad->SetBottomMargin(0.1125);

  hist1->SetTitle("");
  hist1->GetYaxis()->SetMaxDigits(4);
  hist1->GetXaxis()->SetTitleOffset(1.25);
  //  hist1->GetYaxis()->SetMaxDigits(3);

  hist1->DrawCopy("e2");
  hist2->DrawCopy("samese hist");
  if(hist3!=NULL && option!="Final")
    hist3->DrawCopy("samese hist");
  if(hist3!=NULL && option=="Final"){
    hist3->SetFillColor(kRed);
    hist3->DrawCopy("hist sames");
  }
  if(hist3!=NULL &&  option=="FinalEtaCombined"){
    hist3->SetFillColor(kRed);
    hist3->DrawCopy("sames hist");
  }
  if(hist4!=NULL)
    hist4->DrawCopy("samese hist");

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

  unsigned int nx = hist1->GetXaxis()->GetNbins();
  Double_t*  xbins = new Double_t[hist1->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hist1->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hist1->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hist1->GetXaxis()->GetNbins()]=hist1->GetXaxis()->GetBinUpEdge(hist1->GetXaxis()->GetNbins());

  TH1D* ratiodataMC   = new TH1D("ratiodataMC","",nx,xbins);
  
  for(int b=1; b< (2+(ratiodataMC->GetNbinsX())); b++){
    
    double data=hist1->GetBinContent(b);
    if(data==0) continue;
    double mc=hist2->GetBinContent(b);
    if(mc==0) continue;
    
    double dataratio=data/mc;
    double err=hist1->GetBinError(b)/mc;
    ratiodataMC->SetBinContent(b,dataratio);
    ratiodataMC->SetBinError(b,err);
  }

  ratioSettings(ratiodataMC, 0.74, 1.24,"",ratiolabel.c_str(),0.1125,0.1125,0.35,0.15,20,kBlack,0.75,5);
  ratiodataMC->DrawCopy("p");

  TLine *line;
  line = new TLine(hist1->GetXaxis()->GetBinLowEdge(1),1.,hist1->GetXaxis()->GetBinUpEdge(hist1->GetXaxis()->GetNbins()),1.);
  line->SetLineColor(kRed+2);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();
  padr->Update();  

  std::string years="";
  if(config.DataYears=="2017")
    years="2017";
  else if(config.DataYears=="2015+2016")
    years="2015p2016";

  c->cd();
  c->Print((config.OutputFileDir + "Plots/"+ years +"/Multijet/"+kine+"_multijet_"+option+"_"+channel+"_"+cuts+"_cuts"+eta_bin+".pdf").c_str());  
  delete c;

  return;
}


void MultiJetResults::SinglePlot(Config config, std::string kine, std::string option,
			       TH1D *hist1, TH1D *hist2, TH1D *hist3, TH1D *hist4,
			       std::string cuts, std::string eta_bin, std::string label,
				 double ylow, double yhigh, bool ylog, std::string channel){

  std::string xlabel,ylabel; double xlow=0., xhigh=0.;
  TString leglabel="";

  if(kine=="rmass") xlabel="m_{T} [GeV]";
  if(kine=="d0sig") xlabel="d0sig";
  if(kine=="pt")    xlabel="p_{T} [GeV]";
  if(kine=="eta")   xlabel="#||{#eta^{#mu}}";
  if(kine=="met")   xlabel="E_{T}^{miss} [GeV]";
  if(kine=="iso30") xlabel= "iso30";
  if(kine=="eta_d") xlabel="#||{#eta^{#mu}}";
  if(kine=="phi") xlabel="#phi^{#mu} [rad]";
  if(kine=="phi_met") xlabel="#phi^{E_{T}^{miss}} [rad]";
  if(kine=="delta_phi") xlabel="#||{#Delta#phi} [rad]";

  if(channel=="zmumu")  leglabel="Z#rightarrow#mu#mu";
  if(channel=="wplus")  leglabel="W^{+}#rightarrow#mu^{+}#nu";
  if(channel=="wminus") leglabel="W^{-}#rightarrow#mu^{-}#nu";

  ylabel="Events";
  if(option=="DataCompNorm" || option=="DataFitSumNorm") ylabel="Arbitrary Units";

  if(ylog) ylow = 1;
  if(ylog && (option=="DataCompNorm" || option=="DataFitSumNorm") ) ylow = 1e-10;

  //divide bin eta bin to smooth the plot
  // if(kine=="eta"){
  //   for(int k=1;k<hist1->GetXaxis()->GetNbins()+1;k++){
  //     hist1->SetBinContent(k,(hist1   ->GetBinContent(k)/hist1   ->GetBinWidth(k)));
  //     if(hist2!=NULL)
  // 	hist2->SetBinContent(k,(hist2   ->GetBinContent(k)/hist2   ->GetBinWidth(k)));
  //     if(hist3!=NULL)
  // 	hist3->SetBinContent(k,(hist3   ->GetBinContent(k)/hist3   ->GetBinWidth(k)));
  //     if(hist4!=NULL)
  // 	hist4->SetBinContent(k,(hist4   ->GetBinContent(k)/hist4   ->GetBinWidth(k)));
  //   }
  // }

  plotAxisLine(hist1,kBlack,kGreen,20,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  if(hist2!=NULL) 
    plotAxisLine(hist2,kBlue,kGreen,28,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  if(hist3!=NULL){
    if(hist4==NULL) plotAxisLine(hist3,kRed,kGreen,21,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
    if(hist4!=NULL) plotAxisLine(hist3,kViolet+2,kGreen,21,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);
  }

  if(hist4!=NULL) 
    plotAxisLine(hist4,kRed,kGreen,20,0.75,"",xlabel.c_str(),ylabel.c_str(),false,xlow,xhigh,true,ylow,yhigh);

  TString lumi="44.3 fb^{-1}";
  TCanvas *c = new TCanvas("canvas","canvas",550,450);
  TLegend* leg = new TLegend(0.175,0.8,0.375,0.65,"");
  leg->SetBorderSize(0); leg->SetTextSize(0.0375); leg->SetFillColor(0);
  TString naming=config.DataYears;

  TString leg1, leg2, leg3,leg4;
  if(option=="DataMC" || option=="DataFit" || option=="DataSum") leg1="Data";
  if(option=="DataMC") leg2="MC Signal+BG";
  if(option=="DataFit")leg2="Fit";
  if(option=="DataSum")leg2="Weighted sum:MC+MJ";
  if(option=="FitSum"){leg1="Fit"; leg2="Sum";}
  if(option=="DataMCMulti"){leg1="Data"; leg2="MC Signal+BG"; leg3="Multijet";}
  if(option=="DataComp" || option=="DataCompNorm"){leg1="Data"; leg2="Fit"; leg3="MC Signal+BG"; leg4="Multijet";}
  if(option=="Shape") leg1="MultiJet Shape";

  leg->AddEntry(hist1, "#bf{#scale[0.85]{"+leg1+"}}","P");
  if(option=="DataComp" || option=="DataCompNorm"){
    if(hist2!=NULL) leg->AddEntry(hist2, "#bf{#scale[0.85]{"+leg2+"}}","L"); 
    if(hist3!=NULL) leg->AddEntry(hist3, "#bf{#scale[0.85]{"+leg3+"}}","L"); 
    if(hist4!=NULL) leg->AddEntry(hist4, "#bf{#scale[0.85]{"+leg4+"}}","L"); 
  }else{
    if(hist2!=NULL) leg->AddEntry(hist2, "#bf{#scale[0.85]{"+leg2+"}}","P"); 
    if(hist3!=NULL) leg->AddEntry(hist3, "#bf{#scale[0.85]{"+leg3+"}}","P"); 
    if(hist4!=NULL) leg->AddEntry(hist4, "#bf{#scale[0.85]{"+leg4+"}}","P"); 
  }

  TPaveText *box;
  box = new TPaveText(0.2,0.815,0.375,0.915,"NDC");
  box->SetBorderSize(0); box->SetTextSize(0.035); box->SetFillColor(0);
  //box->AddText("#bf{#it{#scale[1.125]{ATLAS}}} Internal");
  box->AddText("CCDY Analysis");
  box->AddText("13 TeV, "+lumi);

  TPaveText *box2;
  box2 = new TPaveText(0.745,0.75,0.795,0.875,"NDC");
  box2->SetBorderSize(0); box2->SetTextSize(0.035); box2->SetFillColor(0);
  box2->AddText(leglabel);
  box2->AddText((cuts+" region").c_str());
  box2->AddText(label.c_str());
  c->cd();

  //main panel
  TPad *pad = new TPad("pad","pad",0,0,1,1);
  if(ylog) pad->SetLogy();
  // pad->SetLogx();
  pad->SetTicks(1,1);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1125);
  pad->SetRightMargin(0.03);
  pad->SetTopMargin(0.05);
  pad->SetBottomMargin(0.1125);

  hist1->SetTitle("");
  hist1->GetYaxis()->SetMaxDigits(4);
  hist1->GetXaxis()->SetTitleOffset(1.25);
  hist1->DrawCopy("e2");
  
  if(option=="DataComp" || option=="DataCompNorm"){
    if(hist2!=NULL) hist2->DrawCopy("hist sames");
    if(hist3!=NULL) hist3->DrawCopy("hist sames");
    if(hist4!=NULL) hist4->DrawCopy("hist sames");
  }else{
    if(hist2!=NULL) hist2->DrawCopy("samese p");
    if(hist3!=NULL) hist3->DrawCopy("samese p");
    if(hist4!=NULL) hist4->DrawCopy("samese p");
  }
  
  pad->Update();
  pad->Modified();
  box->Draw();
  box2->Draw();
  leg->Draw();
  c->cd();

  std::string years="";
  if(config.DataYears=="2017")
    years="2017";
  else if(config.DataYears=="2015+2016")
    years="2015p2016";

  c->Print((config.OutputFileDir + "Plots/"+ years +"/Multijet/"+kine+"_multijet_single_"+option+"_"+channel+"_"+cuts+"_cuts"+eta_bin+".pdf").c_str());  

  return;
}


void MultiJetResults::plotAxisLine(TH1D* hist, int lineColor, int markerColor,
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



void MultiJetResults::ratioSettings(TH1D* hist, double min, double max, TString xlabel, TString ylabel,
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


std::vector<std::string> MultiJetResults::LatexBinName(){

  std::vector<std::string> name; name.clear();

  name.push_back(" Inclusive $\\vert\\eta^{\\mu}\\vert$ ");     	     
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


void MultiJetResults::setstyle(){

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

  //TGaxis::SetMaxDigits(3);

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
