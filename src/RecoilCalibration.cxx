//RecoilCalibration includes
#ifndef RECOILCALIBRATION_CXX
#define RECOILCALIBRATION_CXX

#include "RecoilCalibration.h"
#include <cmath>
#include <iostream>
#include <vector>

#include "TMath.h"

RecoilCalibration::RecoilCalibration()
{
}

RecoilCalibration::~RecoilCalibration() {}

void RecoilCalibration::initialize(TString calibFile, TString channel) {

  calibfile = TFile::Open(calibFile);
  if(calibfile==NULL){
    std::cout<<"no calibration file "<< calibFile << " found ! " << std::endl;
    exit(10); 
    }
  h2_resol = (TH2D*) calibfile -> Get("resol_corr");
  p2_resp = (TProfile2D*) calibfile -> Get("resp_corr");
  p2_upar_data = (TProfile2D*) calibfile -> Get("upar_set_pt_data");
  fux = (TF1*) calibfile->Get("uxcorr");
  fuy = (TF1*) calibfile->Get("uycorr");
  
  h2_resp = (TH2D*) calibfile -> Get("h_bias_data_yz");
  //hux = (TProfile*) calibfile->Get("hd_ux_bin");
  //huy = (TProfile*) calibfile->Get("hd_uy_bin");
  
  if(channel=="z")
    h = (TH2F*)calibfile->Get("SETweightVsPtV_z");
  else if (channel == "wminus")
    h = (TH2F*)calibfile->Get("SETweightVsPtV_wminus");
  else if (channel == "wplus")
    h = (TH2F*)calibfile->Get("SETweightVsPtV_wplus");
  else {
    std::cout<<"no channel "<< channel << " found ! " << std::endl;
    std::cout<<"please choose between wminus, wplus and z" << std::endl;
    exit(10); 
    }
  return;
  
}


void RecoilCalibration::execute(TLorentzVector trueboson, double sumet, TLorentzVector& u, float& weight){
  applySumetReweighting(trueboson.Pt(), sumet, weight);
  applyUX(sumet, u);
  applyUY(sumet, u);
  //applyUXbin(sumet, u);
  //applyUYbin(sumet, u);
  applyResponseAndResolCorr(trueboson, sumet, u);
  return;
}


void RecoilCalibration::applySumetReweighting(double ptv, double sumet, float& weight){
  //std::cout << histname << std::endl;
  if(h==NULL){
    std::cout<<"no histogram SETweightVsPtV found ! " << std::endl;
    exit(10); 
    }
  int binx=h->GetXaxis()->FindBin(sumet);
  if(binx > h->GetNbinsX() ){
    //std::cout<< " Warning : value of SumET off limits : " << sumet << std::endl;
    //std::cout<< " Using last bin" << std::endl;
    binx=h->GetNbinsX();
    }
  int biny=h->GetYaxis()->FindBin(ptv);
  if(biny > h->GetNbinsY() ){
    //std::cout<< " Warning : value of SumET off limits : " << sumet << std::endl;
    //std::cout<< " Using last bin" << std::endl;
    biny=h->GetNbinsY();
    }
  weight=h->GetBinContent(binx, biny);
  return;
}


void RecoilCalibration::applyUX(double sumet, TLorentzVector& u){
  if(fux==NULL){
    std::cout<<"no function uxcorr found ! " << std::endl;
    exit(10); 
    }
  double corr=fux->Eval(sumet);
  //double corr=0.0155947 ;
  u.SetPx( u.Px() + corr*1000. );
  return;
}


void RecoilCalibration::applyUY(double sumet, TLorentzVector& u){
  if(fuy==NULL){
    std::cout<<"no function uycorr found ! " << std::endl;
    exit(10); 
    }
  double corr=fuy->Eval(sumet);
  //double corr=-0.194105 ;
  u.SetPy( u.Py() + corr*1000. );
  return;
}


void RecoilCalibration::applyUXbin(double sumet, TLorentzVector& u){
  if(hux==NULL){
    std::cout<<"no function uxcorr found ! " << std::endl;
    exit(10); 
    }
  double corr=hux->GetBinContent(hux->GetXaxis()->FindBin(sumet));
  //double corr=0.0155947 ;
  u.SetPx( u.Px() + corr );
  return;
}


void RecoilCalibration::applyUYbin(double sumet, TLorentzVector& u){
  if(huy==NULL){
    std::cout<<"no function uycorr found ! " << std::endl;
    exit(10); 
    }
  double corr=huy->GetBinContent(huy->GetXaxis()->FindBin(sumet));
  //double corr=-0.194105 ;
  u.SetPy( u.Py() + corr );
  return;
}

void RecoilCalibration::applyResponseAndResolCorr(TLorentzVector trueboson, double sumet, TLorentzVector& u){
  // if(p2_resp==NULL){
  //   std::cout<<"no function resp_corr found ! " << std::endl;
  //   exit(10); 
  //   }
  if(h2_resol==NULL){
    std::cout<<"no function resol_corr found ! " << std::endl;
    exit(10); 
    }
  //get bias
  /*
  int binx=p2_resp->GetXaxis()->FindBin(trueboson.Pt());
  if(binx > p2_resp->GetNbinsX() ){
    binx=p2_resp->GetNbinsX();
    }
  int biny=p2_resp->GetYaxis()->FindBin(sumet);
  if(biny > p2_resp->GetNbinsY() ){
    biny=p2_resp->GetNbinsY();
    }
    double bias=p2_resp->GetBinContent(binx, biny);*/
  int binx=h2_resp->GetXaxis()->FindBin(trueboson.Pt()/1000.);
  if(binx > h2_resp->GetNbinsX() ){
    binx=h2_resp->GetNbinsX();
    }
  int biny=h2_resp->GetYaxis()->FindBin(sumet);
  if(biny > h2_resp->GetNbinsY() ){
    biny=h2_resp->GetNbinsY();
    }
  double bias=h2_resp->GetBinContent(binx, biny);
  bias=bias*1000.;
  //get upar_data
  binx=p2_upar_data->GetXaxis()->FindBin(sumet);
  if(binx > p2_upar_data->GetNbinsX() ){
    binx=p2_upar_data->GetNbinsX();
    }
  biny=p2_upar_data->GetYaxis()->FindBin(trueboson.Pt()/1000.);
  if(biny > p2_upar_data->GetNbinsY() ){
    biny=p2_upar_data->GetNbinsY();
    }
  double upar_data=p2_upar_data->GetBinContent(binx, biny);
  upar_data=upar_data*1000.;

  //get resol
  binx=h2_resol->GetXaxis()->FindBin(trueboson.Pt()/1000.);
  if(binx > h2_resol->GetNbinsX() ){
    binx=h2_resol->GetNbinsX();
    }
  biny=h2_resol->GetYaxis()->FindBin(sumet);
  if(biny > h2_resol->GetNbinsY() ){
    biny=h2_resol->GetNbinsY();
    }
  double resol=h2_resol->GetBinContent(binx, biny);
  if(resol < 1.e-4)resol =1.;
  
  TVector2 upfoT = TVector2( u.Px(), u.Py() );
  TVector2 bosonT = TVector2( trueboson.Px(), trueboson.Py() );
  TVector2 proj = upfoT.Proj(bosonT);
 
  //upfoT = upfoT.Rotate(-bosonT.Phi() );
  //eq. 12 of
  //https://cds.cern.ch/record/2013274/files/ATL-COM-PHYS-2015-344.pdf
  double factor = proj.Mod() - upar_data;
  if( proj*bosonT < 0.){factor = -1.* (proj.Mod() + upar_data); /*std::cout<<"PUTO"<<std::endl;*/}
  
  
  // double factor = proj.Mod() ;
  // if( proj*bosonT < 0.)factor = -1.* proj.Mod() ;
  
    
  //factor = -1.* (proj.Mod() + upar_data);
  //if( proj*bosonT < 0. )factor = -factor;


  double newUpar = upar_data + bias + factor * resol;  
  //double newUpar = (factor * resol) + bias ;
  double newPhi = proj.Phi();


  if(newUpar * (proj*bosonT) < 0 ){
  	/*
  	std::cout<< "new upar negatove !!, bosonT*proj is = " << proj*bosonT <<std::endl;
	std::cout<< "new upar negatove resol =  = " << resol <<std::endl;
	std::cout<< "new upar negatove factor =  = " << factor <<std::endl;
	std::cout<< "new upar   = " << newUpar <<std::endl;
	*/
	newPhi = proj.Phi() - TMath::Pi();
	
	}
  proj.SetMagPhi( fabs(newUpar), newPhi );

  //std::cout<<"upar_data     "<<upar_data<<"    bias:   "<<bias<<"    factor:   "<<factor<<"   resol:   "<<resol<<"    newPhi   "<<newPhi<<std::endl;
  
  TVector2 bosonperp = TVector2( trueboson.Px(), trueboson.Py() );
  bosonperp = bosonperp.Rotate(TMath::Pi()/2.);
  bosonperp /= bosonperp.Mod();
  TVector2 projPerp = upfoT.Proj(bosonperp);
  //eq.6 of the W mass paper arXiV:1701.07240
  projPerp.SetMagPhi( projPerp.Mod()*resol, projPerp.Phi() );
   
  u.SetPx( (proj+projPerp).Px() );
  u.SetPy( (proj+projPerp).Py() );

  //std::cout<<"Final upfo   "<<u.Px()<<"   "<<u.Py()<<"    "<<u.Phi()<<std::endl; 
  return;
}

//recoilCalibTool->applySETrew( pTV, SET )// extra function


void RecoilCalibration::finalize() {
  calibfile->Close();
  std::cout<<"**** Finalizing Hadronic Recoil Calibration ****"<<std::endl;
  
}


//kinematic functions that could be useful in the future
float RecoilCalibration::deltaR(const float eta1,const float eta2,const float phi1,const float phi2) const
{
  float dPhi=this->deltaPhi(phi1,phi2);
  float dEta=std::fabs(eta1-eta2);
  float dR=std::sqrt(std::pow(dEta,2)+std::pow(dPhi,2));
  return dR;
}

float RecoilCalibration::deltaPhi(const float phi1,const float phi2) const
{
  float dPhi=std::fabs(phi1-phi2);
  if (dPhi>m_pi) dPhi=2.*m_pi-dPhi;
  return dPhi;
}


float RecoilCalibration::deltaEta(const float eta1,const float eta2) const
{
  float dEta=std::fabs(eta1-eta2);

  return dEta;
}


#endif
