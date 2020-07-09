//SampleCaller includes
#ifndef SAMPLECALLER_H
#define SAMPLECALLER_H

#include "ConfigSettings.h"

#include <string>
#include <vector>
#include <iostream>


class SampleCaller{

 public:
  SampleCaller();
  virtual ~SampleCaller();

  virtual void data(Config config, std::vector<std::string>& samples);

  virtual void w_plusmunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void w_plustaunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void w_minusmunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void w_minustaunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void zmumu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void ztautau(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void top(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void diboson(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi);
  virtual void systematics(Config config, std::vector<std::string>& systematics);

};

#endif
