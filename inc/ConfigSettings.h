//Config Settings
#ifndef CONFIGSETTINGS_H
#define CONFIGSETTINGS_H

#include <vector>
#include <string>

struct Config{

  std::string NormPlots;
  std::string DoOnlyPlots;
  std::string DoStackPlots;
  std::string DoTruthPlots;
  std::string DoSFPlots;
  std::string DoCalibPlots;
  std::string DoCutFlow;
  std::string DataYears;
  std::string InputFileDir;
  std::string OutputFileDir;
  long int NumberOfEvents;
  std::string DoCwFactor;
  std::string DoCrossSection;
  std::string DoMultijet;
  std::string MultiFitRegion;
  std::string MultijetVariation;
  std::string OnlyMultijetResults;
  std::string Dod0Cut;
  std::string Applyd0Weight;
  std::string Nod0Shift;
  std::string TruthAnalysis;
  std::string RecoAnalysis;

  std::vector<double> xBinsCw;

  std::string OnlyMC;
  std::string OnlyData;
  std::string OnlyInclusive;
  std::string WZSelection;
  std::string ZChannel;
  std::string WChannel;
  std::string TopChannel;
  std::string DibosonChannel;
  std::string MuFlavor;
  std::string TauFlavor;
  std::string Systematics;
  std::string CutFlowName;
  std::string OnTheFlyPileUp;
  std::string HasRecoilInfo;
  std::string SETCalibration;
  std::string RecoilCalibFileDir;
  std::string PRWFile;
  std::string PRWName;
  std::string InsituCorrection;
  std::string ResolResponse;
  std::string TruthMatching;
  std::string RecoMatching;
  std::string SFVariations;
  std::string RecoilBG;
  bool SysInplots;
  bool PlotMulti;
  bool PlotData;
  bool PlotZmumu;
  bool PlotZtautau;
  bool PlotWplusTaunu;
  bool PlotWplusMunu;
  bool PlotWminTaunu;
  bool PlotWminMunu;
  bool PlotTtbar;
  bool PlotWt_top;
  bool PlotWt_antitop;
  bool PlotSingle_top;
  bool PlotSingle_antitop;    
  bool PlotZzqqll;   
  bool PlotWwqqll;   
  bool PlotWwpqqmlnu;
  bool PlotWwplnumqq;
  bool PlotWzlnuqq;  
  bool PlotZzllll;   
  bool PlotWzlnull;  
  bool PlotWzlnununu;
  bool PlotZzllnunu;

};

#endif
