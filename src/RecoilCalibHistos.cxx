// RecoilCalibHistos includes
#ifndef RECOILCALIBHISTOS_CXX
#define RECOILCALIBHISTOS_CXX

#include "RecoilCalibHistos.h"
#include <cmath>
#include <iostream>
#include <vector>


RecoilCalibHistos::RecoilCalibHistos()
{
}

RecoilCalibHistos::~RecoilCalibHistos() {}

void RecoilCalibHistos::initialize(TString dir, TString name) {

  m_newfile = new TFile(dir+"/RecoilCalibHistos_"+name+".root","RECREATE");
  TH1::SetDefaultSumw2(kTRUE);
  h_mu = new TH1D("h_mu",";#LT#mu#GT;Events",50,0,50);  
  h_Sumet=new TH1D("h_Sumet",";#sum E_{T}^{PFO} [GeV];Entries",1000,0.,2000.);
  h_ptvtrue=new TH1D("h_ptvtrue",";p_{T}^{V, true} [GeV];Entries",200,0.,200.);
  h_ptv=new TH1D("h_ptv",";p_{T}^{V, reco} [GeV];Entries",200,0.,200.);

  h_upar=new TH1D("h_upar",";u_{par}^{V, true}[GeV];Entries",400,-200.,200.);
  h_upar_ptv=new TH1D("h_upar_ptv",";u_{par}^{V, true}+p_{T}^{V, reco} [GeV];Entries",400,-200.,200.);
  h_upar_ptv_n=new TH1D("h_upar_ptv_n",";u_{par}^{V, true}-p_{T}^{V, reco} [GeV];Entries",400,-200.,200.);

  h_uperp=new TH1D("h_uperp",";u_{perp}^{V, true} [GeV];Entries",400,-200.,200.);
  h_upar_rec=new TH1D("h_upar_rec",";u_{par}^{V, rec} [GeV];Entries",400,-200.,200.);
  h_uperp_rec=new TH1D("h_uperp_rec",";u_{perp}^{V, rec}  [GeV];Entries",400,-200.,200.);
  h_ut=new TH1D("h_ut",";u_{T} [GeV];Entries",200,0.,200.);
  h_ux=new TH1D("h_ux",";u_{x} [GeV];Entries",200,-100.,100.);
  h_uy=new TH1D("h_uy",";u_{y} [GeV];Entries",200,-100.,100.);
  h_met=new TH1D("h_met",";E_{T}^{miss} [GeV];Entries",200,0.,200.);
  h_utphi=new TH1D("h_utphi",";#phi(u_{T});Entries",100,-5.,5.);
  h_metphi=new TH1D("h_metphi",";#phi(E_{T}^{miss});Entries",100,-5.,5.);
  
  h_uperp_pfo_vs_set = new TH2D("h_uperp_pfo_vs_set",";u_{perp}^{V, reco} [GeV];#sum E_{T}^{PFO} [GeV]",400,-200.,200.,1000,0.,2000.);
  h_uperp_pfo_vs_pt = new TH2D("h_uperp_pfo_vs_pt",";u_{perp}^{V, reco} [GeV];#p_{T}^{V} [GeV]",400,-200.,200.,200,0.,200.);

  h_upar_pfo_vs_set = new TH2D("h_upar_pfo_vs_set",";u_{par}^{V, reco} [GeV];#sum E_{T}^{PFO} [GeV]",400,-200.,200.,1000,0.,2000.);
  h_set_vs_ptv = new TH2D("h_set_vs_ptv",";#sum E_{T}^{PFO} [GeV];p_{T}^{V, reco} [GeV]",1000,0.,2000.,200,0.,200.);
  h_set_vs_ptvtrue = new TH2D("h_set_vs_ptvtrue",";#sum E_{T}^{PFO} [GeV];p_{T}^{V, true} [GeV]",1000,0.,2000.,200,0.,200.);

  h_bias_vs_set = new TH2D("h_bias_vs_set",";bias [GeV];p_{T}^{V, reco} [GeV]",200,-100.,100.,1000,0.,2000.);
  h_bias_vs_ptv = new TH2D("h_bias_vs_ptv",";bias [GeV];p_{T}^{V, reco} [GeV]",200,-100.,100.,200,0.,200.);

  h_bias_vs_ptvtrue = new TH2D("h_bias_vs_ptvtrue",";bias [GeV];p_{T}^{V, true} [GeV]",200,-100.,100.,200,0.,200.);
  h_set_vs_ut = new TH2D("h_set_vs_ut",";#sum E_{T}^{PFO} [GeV];u_{T} [GeV]",1000,0.,2000.,200,0.,200.);
  
  p_ux_phi = new TProfile("p_ux_phi", ";#phi^{lepton};u_{x} [GeV]",100,-5.,5.);
  p_uy_phi = new TProfile("p_uy_phi", ";#phi^{lepton};u_{y} [GeV]",100,-5.,5.);
  p_ux_set = new TProfile("p_ux_set", ";#sum E_{T}^{PFO} [GeV];u_{x} [GeV]",1000,0.,2000.);
  p_uy_set = new TProfile("p_uy_set", ";#sum E_{T}^{PFO} [GeV];u_{y} [GeV]",1000,0.,2000.);

  p_bias_set = new TProfile("p_bias_set", ";#sum E_{T}^{PFO} [GeV];bias [GeV]",1000,0.,2000.);
  p_bias_pt = new TProfile("p_bias_pt", ";P^{V}_{T} [GeV]; bias [GeV]",200,0.,200.);
  p_uperp_set = new TProfile("p_uperp_set", ";#sum E_{T}^{PFO} [GeV];u_{T} [GeV]",1000,0.,2000.);
  p_uperp_pt = new TProfile("p_uperp_pt", ";p_{T}^{V} [GeV];u_{T} [GeV]",200,0.,200.);
  
  p_uperp_vs_set_vs_ptvtrue = new TProfile2D("p_uperp_vs_set_vs_ptvtrue",
  ";#sum E_{T}^{PFO} [GeV];p_{T}^{V, true} [GeV];u_{perp}^{V, true} [GeV]", 1000,0.,2000.,200,0.,200.);
  p_upar_vs_set_vs_ptvtrue = new TProfile2D("p_upar_vs_set_vs_ptvtrue",
  ";#sum E_{T}^{PFO} [GeV];p_{T}^{V, true} [GeV];u_{par}^{V, true} [GeV]", 1000,0.,2000.,200,0.,200.);  
  p_uperp_vs_set_vs_ptv = new TProfile2D("p_uperp_vs_set_vs_ptv",
  ";#sum E_{T}^{PFO} [GeV];p_{T}^{V, reco} [GeV];u_{perp}^{V, reco} [GeV]", 1000,0.,2000.,200,0.,200.);
  p_upar_vs_set_vs_ptv = new TProfile2D("p_upar_vs_set_vs_ptv",
  ";#sum E_{T}^{PFO} [GeV];p_{T}^{V, reco} [GeV];u_{par}^{V, reco} [GeV]", 1000,0.,2000.,200,0.,200.);
  
  respcorr_vs_set_vs_ptv = new TH3D("respcorr_vs_set_vs_ptv",";bias [GeV];#sum E_{T}^{PFO} [GeV];p_{T}^{V, reco} [GeV]",200,-100.,100.,1000,0.,2000.,200,0.,200.);
  resolcorr_vs_set_vs_ptv = new TH3D("resolcorr_vs_set_vs_ptv",";uperp [GeV];#sum E_{T}^{PFO} [GeV];p_{T}^{V, reco} [GeV]",200,-100.,100.,1000,0.,2000.,200,0.,200.);
  
  std::cout<<"Initializing done"<<std::endl;    

}


void RecoilCalibHistos::execute(TLorentzVector tboson, //std::vector<int> truthCharge,
				std::vector<TLorentzVector> recoLeptons,  //std::vector<int> recoCharge,
				TLorentzVector upfo, double sumet, float mu, double recoWeight) {

  //For W selection it is expected to have 1 selected lepton (*Letptons.size()==1)
  //For Z selection it is expected to have 2 selected lepton (*Letptons.size()==2)
  //charge should be stored in the same order as the Leptons
  float GeV=1000.;

  bool Wselection = false, Zselection = false;
  if(recoLeptons.size()==1) Wselection = true;
  else if(recoLeptons.size()==2) Zselection = true;
  else {
    std::cout << "Error : you should feed the RecoilCalibHisto with lepton vectors of the right size" << std::endl; 
    exit(10);
  }
	//if(isMC && tboson_pt->size() > 0 )boson.SetPtEtaPhiE(tboson_pt->at(0), tboson_eta->at(0),tboson_phi->at(0),tboson_e->at(0));
	//else boson=lepton;
	//for W we use uperp and upar w.r.t. the lepton axis
	//Although, would be nice to have also the true W PT in case of MC...to be implemented
	//Because like this, the bias/resolution vs pTV/SET are not well calculated
   TLorentzVector boson;
   if(Wselection)
     boson = recoLeptons.at(0);
   else if(Zselection)
     boson = recoLeptons.at(0) + recoLeptons.at(1);
    
   TVector2 upfoT = TVector2( upfo.Px(), upfo.Py() );
   TVector2 bosonT = TVector2( boson.Px(), boson.Py() );
   TVector2 tbosonT = TVector2( tboson.Px(), tboson.Py() );
   
   TVector2 bosonperp = TVector2( boson.Px(), boson.Py() );
   bosonperp = bosonperp.Rotate(TMath::Pi()/2.);
   bosonperp /= bosonperp.Mod();
   TVector2 tbosonperp = TVector2( tboson.Px(), tboson.Py() );
   tbosonperp = tbosonperp.Rotate(TMath::Pi()/2.);
   tbosonperp /= tbosonperp.Mod();
     
   double upfo_par_t = (upfoT * tbosonT)/tbosonT.Mod();
   double upfo_perp_t = (upfoT * tbosonperp);
   double upfo_par = (upfoT * bosonT)/bosonT.Mod();
   double upfo_perp = (upfoT * bosonperp);
   
   double bias_pfo_t = -1. *( upfo_par_t + tbosonT.Mod());
   double bias_pfo = -1. *( upfo_par + bosonT.Mod());
   
   TLorentzVector met = upfo;
   if(Wselection)
     met += recoLeptons.at(0);
   if(Zselection)
     met += recoLeptons.at(0) + recoLeptons.at(1);
   met = -1 * met;
   h_mu->Fill(mu, recoWeight);
   h_Sumet->Fill(sumet/GeV, recoWeight);
   h_ptvtrue->Fill(tboson.Pt()/GeV, recoWeight);
   h_ptv->Fill(boson.Pt()/GeV, recoWeight);
   h_upar->Fill(upfo_par_t/GeV, recoWeight);

   h_upar_ptv->Fill(upfo_par/GeV + bosonT.Mod()/GeV, recoWeight);
   h_upar_ptv_n->Fill(upfo_par/GeV - bosonT.Mod()/GeV, recoWeight);

   h_uperp->Fill(upfo_perp_t/GeV, recoWeight);
   h_upar_rec->Fill(upfo_par/GeV, recoWeight);
   h_uperp_rec->Fill(upfo_perp/GeV, recoWeight);
   h_ut->Fill(upfo.Pt()/GeV, recoWeight);
   h_utphi->Fill(upfo.Phi(), recoWeight);

   h_ux->Fill(upfo.Px()/GeV, recoWeight);
   h_uy->Fill(upfo.Py()/GeV, recoWeight);

   h_met->Fill(met.Pt()/GeV, recoWeight);
   h_metphi->Fill(met.Phi(), recoWeight);
   
   h_uperp_pfo_vs_set->Fill(upfo_perp/GeV, sumet/GeV, recoWeight);
   h_uperp_pfo_vs_pt->Fill(upfo_perp/GeV, boson.Pt()/GeV, recoWeight);

   h_upar_pfo_vs_set->Fill(upfo_par/GeV, sumet/GeV, recoWeight);
   h_set_vs_ptvtrue->Fill(sumet/GeV, tboson.Pt()/GeV, recoWeight);
   h_set_vs_ptv->Fill(sumet/GeV, boson.Pt()/GeV, recoWeight);
   h_bias_vs_ptvtrue->Fill(bias_pfo_t/GeV, tboson.Pt()/GeV, recoWeight);

   h_bias_vs_ptv->Fill(-1.*bias_pfo/GeV, boson.Pt()/GeV, recoWeight);
   h_bias_vs_set->Fill(-1.*bias_pfo/GeV, sumet/GeV, recoWeight);

   h_set_vs_ut->Fill(sumet/GeV, upfo.Pt()/GeV, recoWeight);
   
   p_ux_phi->Fill(recoLeptons.at(0).Phi(), upfoT.Px()/GeV, recoWeight);
   p_uy_phi->Fill(recoLeptons.at(0).Phi(), upfoT.Py()/GeV, recoWeight);
   // p_ux_phi->Fill(boson.Phi(), upfoT.Px()/GeV, recoWeight);
   // p_uy_phi->Fill(boson.Phi(), upfoT.Py()/GeV, recoWeight);
   
   p_ux_set->Fill(sumet/GeV, upfoT.Px()/GeV, recoWeight);
   p_uy_set->Fill(sumet/GeV, upfoT.Py()/GeV, recoWeight);

   p_bias_set ->SetErrorOption("s");
   p_bias_pt  ->SetErrorOption("s");
   p_uperp_set->SetErrorOption("s");
   p_uperp_pt ->SetErrorOption("s");

   p_bias_set ->Fill(sumet/GeV, -1.*bias_pfo/GeV, recoWeight);
   p_bias_pt  ->Fill(boson.Pt()/GeV, -1.*bias_pfo/GeV, recoWeight);
   p_uperp_set->Fill(sumet/GeV, upfo_perp/GeV, recoWeight);
   p_uperp_pt ->Fill(boson.Pt()/GeV, upfo_perp/GeV, recoWeight);
   
   p_uperp_vs_set_vs_ptvtrue->Fill(sumet/GeV, tboson.Pt()/GeV, upfo_perp_t/GeV, recoWeight);
   p_upar_vs_set_vs_ptvtrue->Fill(sumet/GeV, tboson.Pt()/GeV, upfo_par_t/GeV, recoWeight);
   p_uperp_vs_set_vs_ptv->Fill(sumet/GeV, boson.Pt()/GeV, upfo_perp/GeV, recoWeight);
   p_upar_vs_set_vs_ptv->Fill(sumet/GeV, boson.Pt()/GeV, upfo_par/GeV, recoWeight);
   
   resolcorr_vs_set_vs_ptv->Fill( upfo_perp/GeV, sumet/GeV, boson.Pt()/GeV , recoWeight);
   respcorr_vs_set_vs_ptv->Fill( -1.*bias_pfo/GeV, sumet/GeV, boson.Pt()/GeV , recoWeight);   
   return;
}
 

void RecoilCalibHistos::scalehists(double scale){
  h_mu->Scale(1./h_mu->Integral()); h_mu->Scale(scale);
  h_Sumet->Scale(1./h_Sumet->Integral()); h_Sumet->Scale(scale);
  h_ptvtrue-> Scale(1./h_ptvtrue->Integral()); h_ptvtrue-> Scale(scale);
  h_ptv-> Scale(1./h_ptv->Integral()); h_ptv-> Scale(scale);
  h_upar-> Scale(1./h_upar->Integral()); h_upar-> Scale(scale);

  h_upar_ptv->Scale(1./h_upar_ptv->Integral()); h_upar_ptv->Scale(scale);
  h_upar_ptv_n->Scale(1./h_upar_ptv_n->Integral()); h_upar_ptv_n->Scale(scale);

  h_uperp-> Scale(1./h_uperp->Integral()); h_uperp-> Scale(scale);
  h_upar_rec-> Scale(1./h_upar_rec->Integral());  h_upar_rec-> Scale(scale);
  h_uperp_rec-> Scale(1./h_uperp_rec->Integral()); h_uperp_rec-> Scale(scale);
  h_ut-> Scale(1./h_ut->Integral()); h_ut-> Scale(scale);
  h_ux-> Scale(1./h_ux->Integral()); h_ux-> Scale(scale);
  h_uy-> Scale(1./h_uy->Integral()); h_uy-> Scale(scale);
  h_utphi-> Scale(1./h_utphi->Integral()); h_utphi-> Scale(scale);
  h_met-> Scale(1./h_met->Integral()); h_met-> Scale(scale);
  h_metphi-> Scale(1./h_metphi->Integral()); h_metphi-> Scale(scale);
  h_uperp_pfo_vs_set-> Scale(1./h_uperp_pfo_vs_set->Integral()); h_uperp_pfo_vs_set-> Scale(scale);
  h_uperp_pfo_vs_pt-> Scale(1./h_uperp_pfo_vs_pt->Integral()); h_uperp_pfo_vs_pt-> Scale(scale);
  h_upar_pfo_vs_set-> Scale(1./h_upar_pfo_vs_set->Integral()); h_upar_pfo_vs_set-> Scale(scale);
  h_set_vs_ptvtrue-> Scale(1./h_set_vs_ptvtrue->Integral()); h_set_vs_ptvtrue-> Scale(scale);
  h_set_vs_ptv-> Scale(1./h_set_vs_ptv->Integral());   h_set_vs_ptv-> Scale(scale);
  h_bias_vs_ptvtrue-> Scale(1./h_bias_vs_ptvtrue->Integral()); h_bias_vs_ptvtrue-> Scale(scale);
  h_bias_vs_ptv-> Scale(1./h_bias_vs_ptv->Integral());  h_bias_vs_ptv-> Scale(scale);
  h_bias_vs_set-> Scale(1./h_bias_vs_set->Integral());  h_bias_vs_set-> Scale(scale);
  h_set_vs_ut->Scale(1./h_set_vs_ut->Integral()); h_set_vs_ut-> Scale(scale);
  /*p_ux_phi-> Scale(scale);
  p_uy_phi-> Scale(scale);
  p_uperp_vs_set_vs_ptvtrue-> Scale(scale);
  p_upar_vs_set_vs_ptvtrue-> Scale(scale);
  p_uperp_vs_set_vs_ptv-> Scale(scale);
  p_upar_vs_set_vs_ptv-> Scale(scale);*/
  respcorr_vs_set_vs_ptv->Scale(1./respcorr_vs_set_vs_ptv->Integral()); respcorr_vs_set_vs_ptv->Scale(scale);
  resolcorr_vs_set_vs_ptv->Scale(1./resolcorr_vs_set_vs_ptv->Integral()); resolcorr_vs_set_vs_ptv->Scale(scale);

}


void RecoilCalibHistos::finalize(){
  
  std::cout<<"**** Finalizing Hadronic Recoil Calibration ****"<<std::endl;
  
  m_newfile->cd();
   
  h_mu-> Write();
  h_Sumet-> Write();
  h_ptvtrue-> Write();
  h_ptv-> Write();
  h_upar-> Write();

  h_upar_ptv->Write();
  h_upar_ptv_n->Write();

  h_uperp-> Write();
  h_upar_rec-> Write();
  h_uperp_rec-> Write();
  h_ut-> Write();
  h_ux-> Write();
  h_uy-> Write();
  h_utphi-> Write();
  h_met-> Write();
  h_metphi-> Write();
  h_uperp_pfo_vs_set-> Write();
  h_uperp_pfo_vs_pt-> Write();
  h_upar_pfo_vs_set-> Write();
  h_set_vs_ptvtrue-> Write();
  h_set_vs_ptv-> Write();
  h_bias_vs_ptvtrue-> Write();

  h_bias_vs_set-> Write();
  h_bias_vs_ptv-> Write();

  h_set_vs_ut-> Write();
  p_ux_phi-> Write();
  p_uy_phi-> Write();
  p_ux_set-> Write();
  p_uy_set-> Write();

  p_bias_set ->Write();
  p_bias_pt  ->Write();
  p_uperp_set->Write();
  p_uperp_pt ->Write();

  p_uperp_vs_set_vs_ptvtrue-> Write();
  p_upar_vs_set_vs_ptvtrue-> Write();
  p_uperp_vs_set_vs_ptv-> Write();
  p_upar_vs_set_vs_ptv-> Write();
  respcorr_vs_set_vs_ptv-> Write();
  resolcorr_vs_set_vs_ptv-> Write();
  
  delete h_mu;
  delete h_Sumet;
  delete h_ptvtrue;
  delete h_ptv;
  delete h_upar;

  delete h_upar_ptv;
  delete h_upar_ptv_n;

  delete h_uperp;
  delete h_upar_rec;
  delete h_uperp_rec;
  delete h_ut;
  delete h_ux;
  delete h_uy;
  delete h_utphi;
  delete h_met;
  delete h_metphi;
  delete h_uperp_pfo_vs_set;
  delete h_uperp_pfo_vs_pt;
  delete h_upar_pfo_vs_set;
  delete h_set_vs_ptvtrue;
  delete h_set_vs_ptv;
  delete h_bias_vs_ptvtrue;

  delete h_bias_vs_set;
  delete h_bias_vs_ptv;

  delete h_set_vs_ut;

  delete p_bias_set;
  delete p_bias_pt;  
  delete p_uperp_set;
  delete p_uperp_pt;

  delete p_ux_phi;
  delete p_uy_phi;
  delete p_ux_set;
  delete p_uy_set;
  delete p_uperp_vs_set_vs_ptvtrue;
  delete p_upar_vs_set_vs_ptvtrue;
  delete p_uperp_vs_set_vs_ptv;
  delete p_upar_vs_set_vs_ptv;
  delete respcorr_vs_set_vs_ptv;
  delete resolcorr_vs_set_vs_ptv;
  m_newfile->Close();
  
}

//kinematic functions that could be useful in the future
float RecoilCalibHistos::deltaR(const float eta1,const float eta2,const float phi1,const float phi2) const
{
  float dPhi=this->deltaPhi(phi1,phi2);
  float dEta=std::fabs(eta1-eta2);
  float dR=std::sqrt(std::pow(dEta,2)+std::pow(dPhi,2));
  return dR;
}

float RecoilCalibHistos::deltaPhi(const float phi1,const float phi2) const
{
  float dPhi=std::fabs(phi1-phi2);
  if (dPhi>TMath::Pi()) dPhi=2.*TMath::Pi()-dPhi;
  return dPhi;
}


float RecoilCalibHistos::deltaEta(const float eta1,const float eta2) const
{
  float dEta=std::fabs(eta1-eta2);

  return dEta;
}

#endif
