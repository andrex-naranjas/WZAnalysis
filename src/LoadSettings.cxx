//Simple class to parse arguments from config file
#ifndef LOADSETTINGS_CXX
#define LOADSETTINGS_CXX

#include "LoadSettings.h"
#include "ConfigSettings.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


LoadSettings::LoadSettings()
{
}

LoadSettings::~LoadSettings(){}

void LoadSettings::loadConfig(Config& config, std::string fileName){

  std::ifstream fin(fileName);
  std::string line;
  while (getline(fin,line)) {
    std::istringstream sin(line.substr(line.find("=") + 1));
    if (line.find("NormPlots") != -1)
       sin >> config.NormPlots;
    else if (line.find("DataYears") != -1)
       sin >> config.DataYears;
    else if (line.find("DoOnlyPlots") != -1)
       sin >> config.DoOnlyPlots;
    else if (line.find("DoStackPlots") != -1)
       sin >> config.DoStackPlots;
    else if (line.find("DoTruthPlots") != -1)
       sin >> config.DoTruthPlots;
    else if (line.find("DoSFPlots") != -1)
       sin >> config.DoSFPlots;
    else if (line.find("DoCalibPlots") != -1)
       sin >> config.DoCalibPlots;
    else if (line.find("DoCutFlow") != -1)
       sin >> config.DoCutFlow;
    else if (line.find("InputFileDir") != -1)
       sin >> config.InputFileDir;
    else if (line.find("OutputFileDir") != -1)
       sin >> config.OutputFileDir;
    else if (line.find("NumberOfEvents") != -1)
       sin >> config.NumberOfEvents;
    else if (line.find("TruthAnalysis") != -1)
       sin >> config.TruthAnalysis;
    else if (line.find("RecoAnalysis") != -1)
       sin >> config.RecoAnalysis;
    else if (line.find("DoCwFactor") != -1)
       sin >> config.DoCwFactor;
    else if (line.find("DoCrossSection") != -1)
       sin >> config.DoCrossSection;
    else if (line.find("DoMultijet") != -1)
       sin >> config.DoMultijet;
    else if (line.find("MultiFitRegion") != -1)
       sin >> config.MultiFitRegion;
    else if (line.find("MultijetVariation") != -1)
       sin >> config.MultijetVariation;
    else if (line.find("OnlyMultijetResults") != -1)
       sin >> config.OnlyMultijetResults;
    else if (line.find("Dod0Cut") != -1)
       sin >> config.Dod0Cut;
    else if (line.find("Applyd0Weight") != -1)
       sin >> config.Applyd0Weight;
    else if (line.find("Nod0Shift") != -1)
       sin >> config.Nod0Shift;
    else if (line.find("xBinsCw") != -1){
      double number;
      while(sin >> number)
	config.xBinsCw.push_back(number);
    }
    else if (line.find("OnlyMC") != -1)
       sin >> config.OnlyMC;
    else if (line.find("OnlyData") != -1)
       sin >> config.OnlyData;
    else if (line.find("OnlyInclusive") != -1)
       sin >> config.OnlyInclusive;
    else if (line.find("WZSelection") != -1)
      sin >> config.WZSelection;
    else if (line.find("ZChannel") != -1)
      sin >> config.ZChannel;
    else if (line.find("WChannel") != -1)
      sin >> config.WChannel;
    else if (line.find("TopChannel") != -1)
      sin >> config.TopChannel;
    else if (line.find("DibosonChannel") != -1)
      sin >> config.DibosonChannel;
    else if (line.find("MuFlavor") != -1)
      sin >> config.MuFlavor;
    else if (line.find("TauFlavor") != -1)
      sin >> config.TauFlavor;
    else if (line.find("Systematics") != -1)
      sin >> config.Systematics;
    else if (line.find("HasRecoilInfo") != -1)
      sin >> config.HasRecoilInfo;
    else if (line.find("OnTheFlyPileUp") != -1)
      sin >> config.OnTheFlyPileUp;
    else if (line.find("SETCalibration") != -1)
      sin >> config.SETCalibration;
    else if (line.find("RecoilCalibFileDir") != -1)
      sin >> config.RecoilCalibFileDir;
    else if (line.find("PRWFile") != -1)
      sin >> config.PRWFile;
    else if (line.find("PRWName") != -1)
      sin >> config.PRWName;
    else if (line.find("InsituCorrection") != -1)
      sin >> config.InsituCorrection;
    else if (line.find("ResolResponse") != -1)
      sin >> config.ResolResponse;
    else if (line.find("TruthMatching") != -1)
      sin >> config.TruthMatching;
    else if (line.find("RecoMatching") != -1)
      sin >> config.RecoMatching;
    else if (line.find("SFVariations") != -1)
      sin >> config.SFVariations;
    else if (line.find("RecoilBG") != -1)
      sin >> config.RecoilBG;
    else if (line.find("SysInplots") != -1)      
      sin >> config.SysInplots;
    else if (line.find("PlotMulti") != -1)      
      sin >> config.PlotMulti;
    else if (line.find("PlotData") != -1)      
      sin >> config.PlotData;
    else if (line.find("PlotZmumu") != -1)      
      sin >> config.PlotZmumu;
    else if (line.find("PlotZtautau") != -1)    
      sin >> config.PlotZtautau;
    else if (line.find("PlotWplusTaunu") != -1) 
      sin >> config.PlotWplusTaunu;
    else if (line.find("PlotWplusMunu") != -1)  
      sin >> config.PlotWplusMunu;
    else if (line.find("PlotWminTaunu") != -1) 
      sin >> config.PlotWminTaunu;
    else if (line.find("PlotWminMunu") != -1)
      sin >> config.PlotWminMunu;
    else if (line.find("PlotTtbar") != -1)
      sin >> config.PlotTtbar;
    else if (line.find("PlotWt_top") != -1)
      sin >> config.PlotWt_top;
    else if (line.find("PlotWt_antitop") != -1)
      sin >> config.PlotWt_antitop;
    else if (line.find("PlotSingle_top") != -1)
      sin >> config.PlotSingle_top;
    else if (line.find("PlotSingle_antitop") != -1)
      sin >> config.PlotSingle_antitop;    
    else if (line.find("PlotZzqqll") != -1)
      sin >> config.PlotZzqqll;
    else if (line.find("PlotWwqqll") != -1)
      sin >> config.PlotWwqqll;   
    else if (line.find("PlotWwpqqmlnu") != -1)
      sin >> config.PlotWwpqqmlnu;
    else if (line.find("PlotWwplnumqq") != -1)
      sin >> config.PlotWwplnumqq;
    else if (line.find("PlotWzlnuqq") != -1)  
      sin >> config.PlotWzlnuqq;  
    else if (line.find("PlotZzllll") != -1)   
      sin >> config.PlotZzllll;   
    else if (line.find("PlotWzlnull") != -1)  
      sin >> config.PlotWzlnull;  
    else if (line.find("PlotWzlnununu") != -1)
      sin >> config.PlotWzlnununu;
    else if (line.find("PlotZzllnunu") != -1) 
      sin >> config.PlotZzllnunu; 

  }//firstwhile
}//loadConfig

#endif
