//CutFlowMaker includes
#ifndef CUTFLOWMAKER_H
#define CUTFLOWMAKER_H

#include "ConfigSettings.h"

#include <string>
#include <vector>
#include <iostream>

#include <TH1.h>
#include <TF1.h>
#include <TFile.h>


class CutFlowMaker{

 public:
  CutFlowMaker();
  virtual ~CutFlowMaker();

  virtual void initialize(Config config);
  virtual void execute();
  virtual void FetchValues(std::string option);

 private:
  virtual void sumFiles(std::vector<std::vector<std::string>> samplesRun, std::string option);
  virtual void loadValues(std::vector<std::vector<std::string>> samplesRun, std::string option);
  virtual void latexConvertor(std::vector<std::vector</*long int*/double>> totalValues, std::vector<std::string> names, std::string option);

  virtual std::string processName(std::string inputName);
  virtual std::vector<std::string> Names(std::string option);

  std::string wzchannel;
  std::string dirYear;
  std::string fileDir;
  
  std::vector<std::vector</*long int*/double>> totalFlow;
  std::vector<std::string> names;

  TFile *fsig, *fwtau, *fzmu, *fztau, *ftop, *fdibo, *fdata;
  
};

#endif
