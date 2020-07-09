#include "MyWZAnalysis.h"
#include "MyWZAnalysisTruth.h"

#include "ConfigSettings.h"
#include "LoadSettings.h"
#include "SampleCaller.h"
#include "PrintPlots.h"
#include "TruthPlots.h"
#include "FileAdd.h"
#include "StackPlots.h"
#include "SysVariations.h"
#include "CwCalculator.h"
#include "CutFlowMaker.h"
#include "CrossSection.h"
#include "SFPlots.h"
#include "MultiJetMaker.h"
#include "MultiJetResults.h"

#include <string>
#include <vector>
#include <iostream>


int main(int argc, char** argv){

  if(argc!=2){
    std::cout<<"Please provide one config file, you have provided "<<argc-1<<" files/arguments"<<std::endl;
    std::cout<<"Try again. BYE!"<<std::endl;
    return 1;
  }

  //Get the cofiguration from configFile.txt
  LoadSettings *settings = new LoadSettings();
  Config config;
  settings->loadConfig(config, argv[1]);

  //lumi normalisation
  double lumiData=0.;
  if(config.DataYears=="2015"){
    lumiData=3219.56;
  }else if(config.DataYears=="2018"){
    lumiData=59937.2;
  }else if(config.DataYears=="2017"){
    lumiData=44307.4;
  }else if(config.DataYears=="2015+2016"){
    lumiData=  32988.1+3219.56;}

  //32988.1
  //3219.56
  //32988.1 2016

  //kinematical cuts
  double eta_cut = 2.4, pt_cut1 = 30000., pt_cut2 = 55000, met_cut1 = 25000., met_cut2 = 55000., mwt_cut1 = 55000., mwt_cut2 = 110000.;
  //eta_cut = 2.4, pt_cut1 = 30000., pt_cut2 = 55000, met_cut1 = 0., met_cut2 = 0., mwt_cut1 = 0., mwt_cut2 = 0.;

  SampleCaller *SampleLumi = new SampleCaller();

  std::vector<std::string> systematics; systematics.clear();
  SampleLumi->systematics(config,systematics);

  std::vector<std::vector<std::string>> samplesRun;
  std::vector<std::vector<double>> lumiRun;
  samplesRun.clear(); lumiRun.clear();

  if(config.OnlyMC!="True"){
    std::vector<std::string> data_sample; std::vector<double> data_lumi;
    SampleLumi->data(config,data_sample);
    for(int i=0;i<data_sample.size();i++) data_lumi.push_back(1.);
    samplesRun.push_back(data_sample);
    lumiRun.push_back(data_lumi);
  }
  if(config.OnlyData!="True"){
    if(config.WChannel == "True" && config.MuFlavor=="True" && config.WZSelection=="wplus"){
	std::vector<std::string> wplusmunu_sample; std::vector<double> wplusmunu_lumi;
	SampleLumi->w_plusmunu(config,lumiData,wplusmunu_sample,wplusmunu_lumi);
	samplesRun.push_back(wplusmunu_sample);
	lumiRun.push_back(wplusmunu_lumi);
    }
    if(config.WChannel == "True" && config.TauFlavor=="True" && config.WZSelection=="wplus"){
      std::vector<std::string> wplustaunu_sample; std::vector<double> wplustaunu_lumi;
      SampleLumi->w_plustaunu(config,lumiData,wplustaunu_sample,wplustaunu_lumi);
      samplesRun.push_back(wplustaunu_sample);
      lumiRun.push_back(wplustaunu_lumi);
    }
    if(config.WChannel == "True" && config.MuFlavor=="True" && config.WZSelection=="wminus"){
      std::vector<std::string> wminusmunu_sample; std::vector<double> wminusmunu_lumi;
      SampleLumi->w_minusmunu(config,lumiData,wminusmunu_sample,wminusmunu_lumi);
      samplesRun.push_back(wminusmunu_sample);
      lumiRun.push_back(wminusmunu_lumi);
    }
    if(config.WChannel == "True" && config.TauFlavor=="True" && config.WZSelection=="wminus"){
      std::vector<std::string> wminustaunu_sample; std::vector<double> wminustaunu_lumi;
      SampleLumi->w_minustaunu(config,lumiData,wminustaunu_sample,wminustaunu_lumi);
      samplesRun.push_back(wminustaunu_sample);
      lumiRun.push_back(wminustaunu_lumi);
    }
    if(config.ZChannel == "True" && config.MuFlavor=="True"){
      std::vector<std::string> zmumu_sample; std::vector<double> zmumu_lumi;
      SampleLumi->zmumu(config,lumiData,zmumu_sample,zmumu_lumi);
      samplesRun.push_back(zmumu_sample);
      lumiRun.push_back(zmumu_lumi);
    }
    if(config.ZChannel == "True" && config.TauFlavor=="True"){
      std::vector<std::string> ztautau_sample; std::vector<double> ztautau_lumi;
      SampleLumi->ztautau(config,lumiData,ztautau_sample,ztautau_lumi);
      samplesRun.push_back(ztautau_sample);
      lumiRun.push_back(ztautau_lumi);
    }
    if(config.TopChannel == "True"){
      std::vector<std::string> top_sample; std::vector<double> top_lumi;
      SampleLumi->top(config,lumiData,top_sample,top_lumi);
      samplesRun.push_back(top_sample);
      lumiRun.push_back(top_lumi);
    }
    if(config.DibosonChannel == "True"){
      std::vector<std::string> diboson_sample; std::vector<double> diboson_lumi;
      SampleLumi->diboson(config,lumiData,diboson_sample,diboson_lumi);
      samplesRun.push_back(diboson_sample);
      lumiRun.push_back(diboson_lumi);
    }
  }  

  if(config.DoOnlyPlots!="True" && config.RecoAnalysis=="True"){
    //run the analysis for every sample
    for(int iRun=0; iRun< samplesRun.size();iRun++){
      for(int iSamples=0; iSamples<samplesRun[iRun].size(); iSamples++){
  	if(samplesRun[iRun][iSamples]=="data15" || samplesRun[iRun][iSamples]=="data16"|| samplesRun[iRun][iSamples]=="data17"|| samplesRun[iRun][iSamples]=="data18" ){
	  
  	  MyWZAnalysis r(config, samplesRun[iRun][iSamples], "nominal", lumiRun[iRun][iSamples],
  			 eta_cut, pt_cut1, pt_cut2, met_cut1, met_cut2, mwt_cut1, mwt_cut2);  
  	  r.Loop();
  	  std::cout<<samplesRun[iRun][iSamples]<<"    "<<lumiRun[iRun][iSamples]<<std::endl;
  	  continue;}
	
  	for(int iSystematics = 0; iSystematics < (int)systematics.size(); ++iSystematics){
  	  MyWZAnalysis r(config, samplesRun[iRun][iSamples], systematics[iSystematics],
  			 lumiRun[iRun][iSamples],
  			 eta_cut, pt_cut1, pt_cut2, met_cut1, met_cut2, mwt_cut1, mwt_cut2);
  	  r.Loop();
	  
  	  std::cout<<samplesRun[iRun][iSamples]<<"   "<<lumiRun[iRun][iSamples]<<"     "<<systematics[iSystematics]<<std::endl;
  	}
      }
    }
  }

  //Multijet
  if(config.DoMultijet=="True"){

    if(config.OnlyMultijetResults!="True"){
      MultiJetMaker *m_multijet = new MultiJetMaker();
      m_multijet->initialize(config);
      m_multijet->execute(config);
      m_multijet->finalize();
    }
    
    MultiJetResults *m_multijet_results = new MultiJetResults();
    m_multijet_results->initialize(config);
    m_multijet_results->execute(config);
    m_multijet_results->finalize();
  }

  //cutflows
  if(config.DoCutFlow=="True"){
    CutFlowMaker *m_cuflow = new CutFlowMaker();
    m_cuflow->initialize(config);
    m_cuflow->execute();
  }

  //Fetch the systematics
  if(config.Systematics =="True" || config.SFVariations=="True") {
    SysVariations *m_sys; m_sys = new SysVariations();
    m_sys->initialize(config);
    m_sys->execute(config);
  }


  //Truth analysis: we only care about the signal MC
  samplesRun.clear(); lumiRun.clear();
  if(config.OnlyData!="True"){
    if(config.WChannel == "True" && config.MuFlavor=="True" && config.WZSelection=="wplus"){
	std::vector<std::string> wplusmunu_sample; std::vector<double> wplusmunu_lumi;
	SampleLumi->w_plusmunu(config,lumiData,wplusmunu_sample,wplusmunu_lumi);
	samplesRun.push_back(wplusmunu_sample);	lumiRun.push_back(wplusmunu_lumi);
    }
    if(config.WChannel == "True" && config.MuFlavor=="True" && config.WZSelection=="wminus"){
      std::vector<std::string> wminusmunu_sample; std::vector<double> wminusmunu_lumi;
      SampleLumi->w_minusmunu(config,lumiData,wminusmunu_sample,wminusmunu_lumi);
      samplesRun.push_back(wminusmunu_sample); lumiRun.push_back(wminusmunu_lumi);
    }
    if(config.ZChannel == "True" && config.MuFlavor=="True"){
      std::vector<std::string> zmumu_sample; std::vector<double> zmumu_lumi;
      SampleLumi->zmumu(config,lumiData,zmumu_sample,zmumu_lumi);
      samplesRun.push_back(zmumu_sample); lumiRun.push_back(zmumu_lumi);
    }
  }

 
  if(config.TruthAnalysis=="True" && config.DoOnlyPlots!="True"){
    
    std::vector<MyWZAnalysisTruth*> runTruth;//clever ;)... this fixed a difficult bug 
    //run truth analysis for every sample
    for(int iRun=0; iRun< samplesRun.size();iRun++){
      for(int iSamples=0; iSamples<samplesRun[iRun].size(); iSamples++){	
	runTruth.push_back(new MyWZAnalysisTruth(config, samplesRun[iRun][iSamples], lumiRun[iRun][iSamples],
						 0., 0., eta_cut, pt_cut1, pt_cut2, met_cut1, met_cut2, mwt_cut1, mwt_cut2));
      }
    }
    for(int k=0; k<runTruth.size();k++) runTruth[k]->Loop();  
  }
  

  if(config.DoCwFactor=="True"){
    CwCalculator *m_cw; m_cw = new CwCalculator();
    m_cw->initialize(config);
    m_cw->execute(config);
    m_cw->finalize();
  }

  if(config.DoCrossSection=="True"){
    CrossSection *m_xsection = new CrossSection();
    m_xsection->initialize(config);
    m_xsection->execute(config);
    m_xsection->finalize();
  }

  if(config.DoCalibPlots=="True"){
    //closure recoil plots
    PrintPlots *m_plots; m_plots = new PrintPlots();
    m_plots->setstyle();
    m_plots->initialize(config);
    m_plots->execute(config);
  }

  if(config.DoSFPlots=="True"){
    SFPlots *m_SFplots = new SFPlots();  
    m_SFplots->setstyle();
    m_SFplots->initialize(config);
    m_SFplots->execute(config);
  }

  if(config.DoTruthPlots=="True"){
    TruthPlots *m_truthplots = new TruthPlots();  
    m_truthplots->setstyle();
    m_truthplots->initialize(config);
    m_truthplots->execute(config);
  }
  
  if(config.DoStackPlots=="True"){ 
    //control plots
    StackPlots *m_stack; m_stack = new StackPlots();
    m_stack->initialize(config);
    m_stack->execute(config);
  }
  
  return 0;
}//end of program (run)


  // if(config.OnlyInclusive!="True" && config.DoOnlyPlots!="True"){
  //   FileAdd *m_fileAdd;
  //   m_fileAdd = new FileAdd();

  //   if(config.ZChannel=="True" && config.MuFlavor=="True"){
  //     m_fileAdd->execute(config,"zmumu",false);
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"zmumu",true);
  //   }   
  //   if(config.ZChannel=="True" && config.TauFlavor=="True"){
  //     m_fileAdd->execute(config,"ztautau",false);
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"ztautau",true);
  //   }
  //   if(config.WChannelP=="True" && config.MuFlavor=="True"){
  //     m_fileAdd->execute(config,"wplusmunu",false);
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"wplusmunu",true);
  //   }
  //   if(config.WChannelP=="True" && config.TauFlavor=="True"){
  //     m_fileAdd->execute(config,"wplustaunu",false);
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"wplustaunu",true);
  //   }
  //   if(config.WChannelN=="True" && config.MuFlavor=="True"){
  //     m_fileAdd->execute(config,"wminmunu",false); 
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"wminmunu",true);
  //   }
  //   if(config.WChannelN=="True" && config.TauFlavor=="True"){
  //     m_fileAdd->execute(config,"wmintaunu",false);
  //     if(config.HasRecoilInfo=="True")
  // 	m_fileAdd->execute(config,"wmintaunu",true);
  //   }
  // }
