//Simple class to handle mass slices samples for W/Z analysis
#ifndef SAMPLECALLER_CXX
#define SAMPLECALLER_CXX

#include "SampleCaller.h"


SampleCaller::SampleCaller()
{
}

SampleCaller::~SampleCaller(){}

void SampleCaller::data(Config config, std::vector<std::string>& samples){

  if(config.DataYears=="2015") samples.push_back("data15");
  if(config.DataYears=="2016") samples.push_back("data16");
  if(config.DataYears=="2017") samples.push_back("data17");
  if(config.DataYears=="2018") samples.push_back("data18");
  if(config.DataYears=="2015+2016"){
    samples.push_back("data15");
    samples.push_back("data16");
  }
  return;
}

void SampleCaller::w_plusmunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("wplusmunu"); lumi.push_back(lumiData*11306.);
  if(config.OnlyInclusive=="False"){
    samples.push_back("wplusmunu_120");  lumi.push_back(lumiData*32.059);      
    samples.push_back("wplusmunu_180");  lumi.push_back(lumiData*5.0033);     
    samples.push_back("wplusmunu_250");  lumi.push_back(lumiData*1.7546);     
    samples.push_back("wplusmunu_400");  lumi.push_back(lumiData*0.31237);        
    samples.push_back("wplusmunu_600");  lumi.push_back(lumiData*0.060788);       
    samples.push_back("wplusmunu_800");  lumi.push_back(lumiData*0.017666);   
    samples.push_back("wplusmunu_1000"); lumi.push_back(lumiData*0.0072888);  
    samples.push_back("wplusmunu_1250"); lumi.push_back(lumiData*0.002507);  
    samples.push_back("wplusmunu_1500"); lumi.push_back(lumiData*0.000986); 
    samples.push_back("wplusmunu_1750"); lumi.push_back(lumiData*0.000425); 
    samples.push_back("wplusmunu_2000"); lumi.push_back(lumiData*0.00019462); 
    samples.push_back("wplusmunu_2250"); lumi.push_back(lumiData*9.3e-05); 
    samples.push_back("wplusmunu_2500"); lumi.push_back(lumiData*4.6e-05); 
    samples.push_back("wplusmunu_2750"); lumi.push_back(lumiData*2.3e-05); 
    samples.push_back("wplusmunu_3000"); lumi.push_back(lumiData*1.8e-05);  
    samples.push_back("wplusmunu_3500"); lumi.push_back(lumiData*5.0963e-06); 
    samples.push_back("wplusmunu_4000"); lumi.push_back(lumiData*1.4305e-06); 
    samples.push_back("wplusmunu_4500"); lumi.push_back(lumiData*4.0127e-07); 
    samples.push_back("wplusmunu_5000"); lumi.push_back(lumiData*1.5342e-07); 
  }
  return;
}

void SampleCaller::w_plustaunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){
  
  samples.push_back("wplustaunu");  lumi.push_back(lumiData*11306.);
  if(config.OnlyInclusive=="False"){
    samples.push_back("wplustaunu_120");   lumi.push_back(lumiData*32.059);
    samples.push_back("wplustaunu_180");   lumi.push_back(lumiData*5.003);     
    samples.push_back("wplustaunu_250");   lumi.push_back(lumiData*1.7545);     
    samples.push_back("wplustaunu_400");   lumi.push_back(lumiData*0.31237);    
    samples.push_back("wplustaunu_600");   lumi.push_back(lumiData*0.060788);   
    samples.push_back("wplustaunu_800");   lumi.push_back(lumiData*0.017666);   
    samples.push_back("wplustaunu_1000");  lumi.push_back(lumiData*0.007289);  
    samples.push_back("wplustaunu_1250");  lumi.push_back(lumiData*0.002507);  
    samples.push_back("wplustaunu_1500");  lumi.push_back(lumiData*0.000986); 
    samples.push_back("wplustaunu_1750");  lumi.push_back(lumiData*0.000425); 
    samples.push_back("wplustaunu_2000");  lumi.push_back(lumiData*0.000195); 
    samples.push_back("wplustaunu_2250");  lumi.push_back(lumiData*9.3e-05);  
    samples.push_back("wplustaunu_2500");  lumi.push_back(lumiData*4.6e-05);  
    samples.push_back("wplustaunu_2750");  lumi.push_back(lumiData*2.3e-05);  
    samples.push_back("wplustaunu_3000");  lumi.push_back(lumiData*1.8e-05);  
    samples.push_back("wplustaunu_3500");  lumi.push_back(lumiData*5.0e-06); 
    samples.push_back("wplustaunu_4000");  lumi.push_back(lumiData*1.0e-06); 
    samples.push_back("wplustaunu_4500");  lumi.push_back(lumiData*4.0127e-07); 
    samples.push_back("wplustaunu_5000");  lumi.push_back(lumiData*1.5346e-07); 
    }
  return;
}


void SampleCaller::w_minusmunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("wminmunu"); lumi.push_back(lumiData*8282.9);
  if(config.OnlyInclusive=="False"){
    samples.push_back("wminmunu_120");   lumi.push_back(lumiData*22.194);	  
    samples.push_back("wminmunu_180");   lumi.push_back(lumiData*3.2849);	   
    samples.push_back("wminmunu_250");   lumi.push_back(lumiData*1.0831);	   
    samples.push_back("wminmunu_400");   lumi.push_back(lumiData*0.1754);	   
    samples.push_back("wminmunu_600");   lumi.push_back(lumiData*0.030979);	   
    samples.push_back("wminmunu_800");   lumi.push_back(lumiData*0.008286);  
    samples.push_back("wminmunu_1000");  lumi.push_back(lumiData*0.003159);  
    samples.push_back("wminmunu_1250");  lumi.push_back(lumiData*0.001003);  
    samples.push_back("wminmunu_1500");  lumi.push_back(lumiData*0.000368); 
    samples.push_back("wminmunu_1750");  lumi.push_back(lumiData*0.000149); 
    samples.push_back("wminmunu_2000");  lumi.push_back(lumiData*6.5e-05); 
    samples.push_back("wminmunu_2250");  lumi.push_back(lumiData*3.0e-05); 
    samples.push_back("wminmunu_2500");  lumi.push_back(lumiData*1.4549e-05); 
    samples.push_back("wminmunu_2750");  lumi.push_back(lumiData*7.0e-06); 
    samples.push_back("wminmunu_3000");  lumi.push_back(lumiData*6.0e-06); 
    samples.push_back("wminmunu_3500");  lumi.push_back(lumiData*1.5975e-06); 
    samples.push_back("wminmunu_4000");  lumi.push_back(lumiData*4.721e-07);  
    samples.push_back("wminmunu_4500");  lumi.push_back(lumiData*1.4279e-07); 
    samples.push_back("wminmunu_5000");  lumi.push_back(lumiData*6.1624e-08);    
  }
  return;
}

void SampleCaller::w_minustaunu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("wmintaunu");   lumi.push_back(lumiData*8282.6);
  if(config.OnlyInclusive=="False"){
    samples.push_back("wmintaunu_120");  lumi.push_back(lumiData*22.194);	  
    samples.push_back("wmintaunu_180");  lumi.push_back(lumiData*3.2849);	   
    samples.push_back("wmintaunu_250");  lumi.push_back(lumiData*1.0831);	   
    samples.push_back("wmintaunu_400");  lumi.push_back(lumiData*0.1754);	   
    samples.push_back("wmintaunu_600");  lumi.push_back(lumiData*0.030979);	   
    samples.push_back("wmintaunu_800");  lumi.push_back(lumiData*0.008286); 
    samples.push_back("wmintaunu_1000"); lumi.push_back(lumiData*0.003159);  
    samples.push_back("wmintaunu_1250"); lumi.push_back(lumiData*0.001003); 
    samples.push_back("wmintaunu_1500"); lumi.push_back(lumiData*0.000368);
    samples.push_back("wmintaunu_1750"); lumi.push_back(lumiData*0.000149); 
    samples.push_back("wmintaunu_2000"); lumi.push_back(lumiData*6.5e-05);
    samples.push_back("wmintaunu_2250"); lumi.push_back(lumiData*3.0e-05);
    samples.push_back("wmintaunu_2500"); lumi.push_back(lumiData*1.5e-05);
    samples.push_back("wmintaunu_2750"); lumi.push_back(lumiData*7.0e-06);
    samples.push_back("wmintaunu_3000"); lumi.push_back(lumiData*6.0e-06);
    samples.push_back("wmintaunu_3500"); lumi.push_back(lumiData*2e-06);
    samples.push_back("wmintaunu_4000"); lumi.push_back(lumiData*4.7207e-07); 
    samples.push_back("wmintaunu_4500"); lumi.push_back(lumiData*1.4279e-07);
    samples.push_back("wmintaunu_5000"); lumi.push_back(lumiData*6.16194e-08);
  }
  return;
}


void SampleCaller::zmumu(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("zmumu");  lumi.push_back(lumiData*1901.1);
  if(config.OnlyInclusive=="False"){
    samples.push_back("zmumu_120");   lumi.push_back(lumiData*17.478);    
    samples.push_back("zmumu_180"); 	lumi.push_back(lumiData*2.9215);    
    samples.push_back("zmumu_250"); 	lumi.push_back(lumiData*1.0819);     
    samples.push_back("zmumu_400"); 	lumi.push_back(lumiData*0.19551);    
    samples.push_back("zmumu_600"); 	lumi.push_back(lumiData*0.037403);  
    samples.push_back("zmumu_800"); 	lumi.push_back(lumiData*0.010607);  
    samples.push_back("zmumu_1000");	lumi.push_back(lumiData*0.004259); 
    samples.push_back("zmumu_1250");	lumi.push_back(lumiData*0.001422); 
    samples.push_back("zmumu_1500");	lumi.push_back(lumiData*0.000545);
    samples.push_back("zmumu_1750");	lumi.push_back(lumiData*0.00023);
    samples.push_back("zmumu_2000");	lumi.push_back(lumiData*0.000104);
    samples.push_back("zmumu_2250");	lumi.push_back(lumiData*4.9e-05);  
    samples.push_back("zmumu_2500");	lumi.push_back(lumiData*2.4e-05);
    samples.push_back("zmumu_2750");	lumi.push_back(lumiData*1.3e-05);
    samples.push_back("zmumu_3000");	lumi.push_back(lumiData*1.0e-05);
    samples.push_back("zmumu_3500");	lumi.push_back(lumiData*3.0e-06);
    samples.push_back("zmumu_4000");	lumi.push_back(lumiData*1.04e-06);
    samples.push_back("zmumu_4500");	lumi.push_back(lumiData*2.8071e-07);
    if(config.DataYears!="2015+2016")//missing in mc16a
      samples.push_back("zmumu_5000");	lumi.push_back(lumiData*1.2649e-07);
  }
  return;
}

void SampleCaller::ztautau(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("ztautau");    lumi.push_back(lumiData*1901.1);
  if(config.OnlyInclusive=="False"){
    samples.push_back("ztautau_120");  lumi.push_back(lumiData*17.476);    
    samples.push_back("ztautau_180");  lumi.push_back(lumiData*2.9213);	  
    samples.push_back("ztautau_250");  lumi.push_back(lumiData*1.0819);	  
    samples.push_back("ztautau_400");  lumi.push_back(lumiData*0.1955);	  
    samples.push_back("ztautau_600");  lumi.push_back(lumiData*0.037402);  
    samples.push_back("ztautau_800");  lumi.push_back(lumiData*0.010607);  
    samples.push_back("ztautau_1000"); lumi.push_back(lumiData*0.004259);  
    samples.push_back("ztautau_1250"); lumi.push_back(lumiData*0.001422);  
    samples.push_back("ztautau_1500"); lumi.push_back(lumiData*0.000545);
    samples.push_back("ztautau_1750"); lumi.push_back(lumiData*0.00023);
    samples.push_back("ztautau_2000"); lumi.push_back(lumiData*0.000104);
    samples.push_back("ztautau_2250"); lumi.push_back(lumiData*4.9e-05);
    samples.push_back("ztautau_2500"); lumi.push_back(lumiData*2.4e-05);
    samples.push_back("ztautau_2750"); lumi.push_back(lumiData*1.2e-05); 
    samples.push_back("ztautau_3000"); lumi.push_back(lumiData*1.0e-05);
    samples.push_back("ztautau_3500"); lumi.push_back(lumiData*3.0e-06);
    samples.push_back("ztautau_4000"); lumi.push_back(lumiData*1.0e-06);
    samples.push_back("ztautau_4500"); lumi.push_back(lumiData*2.8071e-07);
    samples.push_back("ztautau_5000");	lumi.push_back(lumiData*1.2648e-07);
  }
  return;
}

void SampleCaller::top(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("ttbar");              lumi.push_back(lumiData*452.3111);
  samples.push_back("single_top_s");	   lumi.push_back(lumiData*2.061459); 	
  //samples.push_back("single_antitop_s");   lumi.push_back(lumiData*1.288566)	
  samples.push_back("wt_top");	           lumi.push_back(lumiData*35.84952); 
  samples.push_back("wt_antitop");         lumi.push_back(lumiData*35.84676);	
  // samples.push_back("wt_top_dilep");	   lumi.push_back(lumiData*3.777565); 	
  // samples.push_back("wt_antitop_dilep");   lumi.push_back(lumiData*3.777431); 	
  // samples.push_back("single_top_t");	   lumi.push_back(lumiData*44.15473);
  // samples.push_back("single_antitop_t");   lumi.push_back(lumiData*26.27516);
  
  return;
}

void SampleCaller::diboson(Config config, double lumiData, std::vector<std::string>& samples, std::vector<double>& lumi){

  samples.push_back("zzqqll");        lumi.push_back(lumiData*2.192389);
  samples.push_back("wzqqll");        lumi.push_back(lumiData*3.4332);
  samples.push_back("wwpqqmlnu");     lumi.push_back(lumiData*24.71);
  samples.push_back("wwplnumqq");     lumi.push_back(lumiData*24.728);
  samples.push_back("wzlnuqq");       lumi.push_back(lumiData*11.418);
  samples.push_back("zzllll");	      lumi.push_back(lumiData*1.2519);
  samples.push_back("wzlnull");	      lumi.push_back(lumiData*4.5786);
  samples.push_back("zzllnunu");      lumi.push_back(lumiData*12.5);
  samples.push_back("wzlnununu");     lumi.push_back(lumiData*3.235);
  return;
}

void SampleCaller::systematics(Config config, std::vector<std::string>& systematics){

  systematics.push_back("nominal");
  if(config.Systematics=="True"){
    systematics.push_back("MET_SoftTrk_ResoPara");
    systematics.push_back("MET_SoftTrk_ResoPerp");
    systematics.push_back("MET_SoftTrk_ScaleDown");
    systematics.push_back("MET_SoftTrk_ScaleUp");
    systematics.push_back("MUON_ID__1down");
    systematics.push_back("MUON_ID__1up");
    systematics.push_back("MUON_MS__1down");
    systematics.push_back("MUON_MS__1up");
    systematics.push_back("MUON_SAGITTA_RESBIAS__1down");
    systematics.push_back("MUON_SAGITTA_RESBIAS__1up");
    systematics.push_back("MUON_SAGITTA_RHO__1down");
    systematics.push_back("MUON_SAGITTA_RHO__1up");
    systematics.push_back("MUON_SCALE__1down");
    systematics.push_back("MUON_SCALE__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_EtaIntercalibration_NonClosure__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_EtaIntercalibration_NonClosure__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_1__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_1__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_2__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_2__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_3__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_GroupedNP_3__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_DataVsMC__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_DataVsMC__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_1__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_1__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_2__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_2__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_3__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_3__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_4__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_4__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_5__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_5__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_6__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_6__1up");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_7restTerm__1down");
    systematics.push_back("JET_StrongReduction_Scenario1_JET_JER_EffectiveNP_7restTerm__1up");
        
  } 
  return;
}

#endif

    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_1__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_1__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_3__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_3__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_2__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_2__1up");

    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Detector1__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical1__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical4__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed2__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical3__1down");
    // systematics.push_back("JET_CategoryReduction_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling1__1down");
    // systematics.push_back("JET_CategoryReduction_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Flavor_Response__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical6__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling3__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling4__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed3__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling2__1down");
    // systematics.push_back("JET_CategoryReduction_JET_BJES_Response__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical2__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed1__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical5__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_negEta__1down");
    
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling1__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical5__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical4__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical3__1up");
    // systematics.push_back("JET_CategoryReduction_JET_BJES_Response__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical2__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Flavor_Response__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed2__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_CategoryReduction_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("ET_JER_SINGLE_NP__1up");
    // systematics.push_back("JET_CategoryReduction_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling4__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical1__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Statistical6__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_CategoryReduction_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Detector1__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed1__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling3__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Modelling2__1up");
    // systematics.push_back("JET_CategoryReduction_JET_EffectiveNP_Mixed3__1up");

    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_7__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_5__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Flavor_Response__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_4__1down");
    // systematics.push_back("JET_GlobalReduction_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_GlobalReduction_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_negEta__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_2__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_1__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_6__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_3__1down");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_8restTerm__1down");
    // systematics.push_back("JET_GlobalReduction_JET_BJES_Response__1down");
  
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_GlobalReduction_JET_BJES_Response__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Flavor_Response__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_6__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_2__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_8restTerm__1up");
    // systematics.push_back("JET_GlobalReduction_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_3__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_4__1up");
    // systematics.push_back("JET_GlobalReduction_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_7__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_GlobalReduction_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_5__1up");
    // systematics.push_back("JET_GlobalReduction_JET_EffectiveNP_1__1up");

    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat11__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_MC__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Jvt__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat10__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical1__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat9__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Veto__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat8__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat7__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed2__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat4__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat1__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat12__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat11__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat3__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat10__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_BJES_Response__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling2__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Generator__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_1__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat6__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_GamESZee__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_negEta__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical5__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_BJES_Response__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Jvt__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_3__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_dPhi__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_3__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_2__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_5__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat8__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat13__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat6__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Veto__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat12__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat5__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat9__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat1__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_Flavor_Response__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_OOC__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat3__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat2__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_dPhi__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat7__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_2__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling1__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical4__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_1__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Detector1__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed1__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling3__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat13__1down");
    // systematics.push_back("JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical6__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_7__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical3__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical2__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_BJES_Response__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_8restTerm__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_negEta__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Purity__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Flavor_Response__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed3__1down");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling4__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_ElecESZee__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Flavor_Response__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat4__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat14__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat2__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_negEta__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_ElecEsmear__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_6__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_KTerm__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_GamEsmear__1down");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_4__1down");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat5__1down");

    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_MC__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Jvt__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_dPhi__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Veto__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_GamEsmear__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat7__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_KTerm__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat3__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Flavor_Response__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat1__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat13__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat3__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat10__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed1__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat14__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat5__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_dPhi__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat8__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling3__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat11__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_3__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat11__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical3__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_3__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_4__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_5__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_8restTerm__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_1__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Veto__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat9__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Flavor_Response__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat4__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_EtaIntercalibration_NonClosure__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat2__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_1__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat6__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat10__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat4__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat13__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Detector1__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical6__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_2__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat2__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_ElecEsmear__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_7__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EffectiveNP_6__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat12__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Flavor_Response__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical4__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling1__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical2__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_SR_Scenario1_JET_GroupedNP_2__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_BJES_Response__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling2__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed2__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_GamESZee__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat5__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat6__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical1__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_BJES_Response__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Mixed3__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Modelling4__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EffectiveNP_Statistical5__1up");
    // //JET_JER_SINGLE_NP__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_ElecESZee__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_BJES_Response__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat12__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_OOC__1up");
    // systematics.push_back("JET_SR_Scenario2_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_SR_Scenario3_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat7__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat8__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Stat9__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Zjet_Stat1__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Generator__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Jvt__1up");
    // systematics.push_back("JET_SR_Scenario4_JET_Gjet_Purity__1up");
    
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat3__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat1__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_KTerm__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Jvt__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_ElecESZee__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_SingleParticle_HighPt__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_PunchThrough_MC16__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_highE__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat7__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat7__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_dPhi__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_GamEsmear__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_OffsetNPV__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_negEta__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Generator__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Flavor_Composition__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Purity__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_PtTerm__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat13__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat4__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_RhoTopology__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat2__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat12__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_OffsetMu__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat11__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat5__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat9__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_ElecEsmear__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat1__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_MC__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat10__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Veto__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat6__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_Modelling__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat13__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_OOC__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat11__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat8__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_BJES_Response__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat2__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat9__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_posEta__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Veto__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat10__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat3__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat6__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat12__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_TotalStat__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat5__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Flavor_Response__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat14__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat4__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_GamESZee__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat8__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_dPhi__1down");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Jvt__1down");
    
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat1__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_TotalStat__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat11__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Flavor_Response__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat7__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_KTerm__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat10__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_ElecEsmear__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Flavor_Composition__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_ElecESZee__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat11__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Jvt__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_RhoTopology__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_PtTerm__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat13__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat14__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_OffsetNPV__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Pileup_OffsetMu__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Generator__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat12__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat1__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_PunchThrough_MC16__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat3__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat2__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat7__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat6__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat13__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_Modelling__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Purity__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat3__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat2__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat5__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat4__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat8__1up");
    // //JET_JER_SINGLE_NP__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_dPhi__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_posEta__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat4__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat6__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat10__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_dPhi__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_SingleParticle_HighPt__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_BJES_Response__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat8__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_negEta__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_GamEsmear__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Jvt__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Veto__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_EtaIntercalibration_NonClosure_highE__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_OOC__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_GamESZee__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat12__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat5__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_MC__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Zjet_Stat9__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Stat9__1up");
    // systematics.push_back("JET_AllNuisanceParameters_JET_Gjet_Veto__1up");
