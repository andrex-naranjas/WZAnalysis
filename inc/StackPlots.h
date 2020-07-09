//StackPlots Includes
#ifndef STACKPLOTS_H
#define STACKPLOTS_H

#include "ConfigSettings.h"

#include "TFile.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <THStack.h>
#include <TString.h>


class StackPlots{

 public:
  StackPlots();
  virtual ~StackPlots();

  virtual void initialize(Config config);
  virtual void execute(Config config);

 private:
  virtual void GetHistos(Config config);
  
  std::vector<std::string> kine;

  TFile *grzmumuP, *grztautauP, *grwplustaunu, *grwplusmunu, *grttbarP,  *grwt_topP, *grwt_antitopP, *grsingle_topP, *grsingle_antitopP, *grzzqqllP, *grwzqqllP, *grwwpqqmlnuP, *grwwplnumqqP, *grwzlnuqqP, *grzzllllP, *grwzlnullP, *grwzlnununuP, *grzzllnunuP;
  
  TFile *grzmumuN, *grztautauN, *grwmintaunu, *grwminmunu, *grttbarN, *grwt_topN, *grwt_antitopN, *grsingle_topN, *grsingle_antitopN, *grzzqqllN, *grwzqqllN, *grwwpqqmlnuN, *grwwplnumqqN, *grwzlnuqqN, *grzzllllN, *grwzlnullN, *grwzlnununuN, *grzzllnunuN;

  TFile *grzmumuZ, *grztautauZ, *grttbarZ, *grwt_topZ, *grwt_antitopZ, *grsingle_topZ, *grsingle_antitopZ, *grzzqqllZ, *grwzqqllZ, *grwwpqqmlnuZ, *grwwplnumqqZ, *grwzlnuqqZ, *grzzllllZ, *grwzlnullZ, *grwzlnununuZ, *grzzllnunuZ;

  TFile  *grDataP, *grSysP, *grDataN, *grSysN, *grDataZ, *grSysZ;

  std::string d0Name="";
  
};

#endif //> !STACKPLOTS_H
