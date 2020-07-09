//Class to calculate the WZ C unfolding factor  
#ifndef CWCALCULATOR_CXX
#define CWCALCULATOR_CXX

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "CwCalculator.h"

CwCalculator::CwCalculator()
{
}

CwCalculator::~CwCalculator(){}

void CwCalculator::initialize(Config config){

  std::string dirYear="";
  if(config.DataYears=="2015+2016") dirYear="2015p2016/";
  if(config.DataYears=="2017") dirYear="2017/";
  if(config.DataYears=="2018") dirYear="2018/";

  std::string dirInclusive, total;
  if(config.OnlyInclusive=="True"){dirInclusive=""; total="";}
  if(config.OnlyInclusive!="True"){dirInclusive="Add/"; total="_Total";}

  std::string puname="";
  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";
  std::string systematic;
  systematic="_nominal";
 
  std::string wzchannel="", nameOfSample="";
  if(config.WZSelection=="zmumu"){
    wzchannel="z"; nameOfSample="zmumu";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus"; nameOfSample="wplusmunu";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus"; nameOfSample="wminmunu";}

  fmcreco   = new TFile((config.OutputFileDir+"Files/"+ dirYear + dirInclusive + nameOfSample +"_" + wzchannel + systematic + puname +total+".root").c_str());
  fmctruth  = new TFile((config.OutputFileDir+"Files/"+ dirYear +"Truth/" + dirInclusive + nameOfSample+"_" + wzchannel +total+".root").c_str());
  fout      = new TFile((config.OutputFileDir+"Files/"+ dirYear +"Truth/" + dirInclusive +"cw_factor_"+nameOfSample+total+".root").c_str(),"recreate");

  if(!fmcreco->IsOpen() || !fmctruth->IsOpen() || !fout->IsOpen()) {std::cout<<"Cw factor on, files not found, change config file, bye!"<<std::endl; exit(1);}

  return;
}

void CwCalculator::execute(Config config){
  
  hReco  = (TH1D*)fmcreco  ->GetObjectUnchecked("hMassCrec");
  hTruth = (TH1D*)fmctruth ->GetObjectUnchecked("hMassCgen");

  hCw_simple  = (TH1D*)fmcreco ->GetObjectUnchecked("hMassCrec");
  hCw_simple->Divide(hTruth);

  hMwtCstayTruth   = (TH1D*)fmctruth ->GetObjectUnchecked("hMassCstayG");
  hMwtCleaveTruth  = (TH1D*)fmctruth ->GetObjectUnchecked("hMassCleave");

  hMwtCstayReco  = (TH1D*)fmcreco ->GetObjectUnchecked("hMassCsta");
  hMwtCcomeReco  = (TH1D*)fmcreco ->GetObjectUnchecked("hMassCcome");

  unsigned int nx = hReco->GetXaxis()->GetNbins();
  double* xbins = new double[hReco->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hReco->GetXaxis()->GetNbins()+1; i++)  xbins[i]=hReco->GetXaxis()->GetBinLowEdge(i+1);
  xbins[hReco->GetXaxis()->GetNbins()]=hReco->GetXaxis()->GetBinUpEdge(hReco->GetXaxis()->GetNbins());

  hCwz     = new TH1D("hCwz"   ,"Cw",nx,xbins);
  hPtywz   = new TH1D("hPtywz" ,"Pu",nx,xbins);
  hStabwz  = new TH1D("hStabwz","St",nx,xbins);

  double sum_reco=0., sum_truth=0.;
  for(int b=1; b< (2+(hCwz->GetNbinsX())); b++){

    double mcreco=hReco->GetBinContent(b);
    if(mcreco==0) continue;
    double mctruth=hTruth->GetBinContent(b);
    if(mctruth==0) continue;

    double mcstayReco = hMwtCstayReco->GetBinContent(b);
    if(mcstayReco == 0) continue;
    double mcstayTruth = hMwtCstayTruth->GetBinContent(b);
    if(mcstayTruth == 0) continue;
    
    sum_reco+=mcreco;
    sum_truth+=mctruth;
    
    hCwz->SetBinContent(b, mcreco/mctruth);
    hCwz->SetBinError(b,0);
    
    hPtywz->SetBinContent(b, mcstayTruth/mcreco);
    hPtywz->SetBinError(b,0);
    
    hStabwz->SetBinContent(b, mcstayTruth/mctruth);
    hStabwz->SetBinError(b,0);
  }


  std::cout<<" Cw:   "<<sum_reco<<"    "<<sum_truth<<"     "<<sum_reco/sum_truth<<std::endl;

  for (int iMwt = 1; iMwt<= nx; iMwt++){

    double xlow = xbins[iMwt-1];
    double xhi  = xbins[iMwt];

    double C    = hCwz->GetBinContent(iMwt);
    double dC   = hCwz->GetBinError(iMwt);

    double Nrec = hMwtCstayReco->GetBinContent(iMwt)  + hMwtCcomeReco->GetBinContent(iMwt);
    double Ngen = hMwtCstayTruth->GetBinContent(iMwt) + hMwtCleaveTruth->GetBinContent(iMwt);

    // errors on number of events in stay, come , leave categories
    double dStay2  = pow(hMwtCstayTruth->GetBinError(iMwt), 2);
    double dCome2  = pow(hMwtCcomeReco->GetBinError(iMwt), 2);
    double dLeave2 = pow(hMwtCleaveTruth->GetBinError(iMwt), 2);

    //contributions to delta_C
    double dCstay2 = pow( (Ngen-Nrec)/(Ngen*Ngen), 2);
    double dCcome2 = pow( (1.0/Ngen), 2);
    double dCleave2= pow( Nrec/(Ngen*Ngen), 2);
    dCstay2  = (dCstay2 * dStay2);
    dCcome2  = (dCcome2 * dCome2);
    dCleave2 = (dCleave2* dLeave2);
    double dCnew = sqrt(dCstay2 + dCcome2 + dCleave2);
    hCwz->SetBinError(iMwt,dCnew);
    printf("i=%2i  M=[%4.0f-%4.0f]  C=%6.4f  dC=%6.4f  dCnew=%6.4f  Nrec=%12.2f   Ngen=%12.2f\n",iMwt,xlow,xhi,C,dC,dCnew,Nrec,Ngen);
  }


  ///repeat for eta...
  hRecoEta  = (TH1D*)fmcreco  ->GetObjectUnchecked("hEtaCrec");
  hTruthEta = (TH1D*)fmctruth ->GetObjectUnchecked("hEtaCgen");

  hCw_eta_simple = (TH1D*)fmcreco  ->GetObjectUnchecked("hEtaCrec");
  hCw_eta_simple->Divide(hTruthEta);

  hEtaCstayTruth   = (TH1D*)fmctruth ->GetObjectUnchecked("hEtaCstayG");
  hEtaCleaveTruth  = (TH1D*)fmctruth ->GetObjectUnchecked("hEtaCleave");

  hEtaCstayReco  = (TH1D*)fmcreco ->GetObjectUnchecked("hEtaCsta");
  hEtaCcomeReco  = (TH1D*)fmcreco ->GetObjectUnchecked("hEtaCcome");

  unsigned int nxEta = hRecoEta->GetXaxis()->GetNbins();
  double* xbinsEta = new double[hRecoEta->GetXaxis()->GetNbins()+1];
  for(unsigned int i=0; i< hRecoEta->GetXaxis()->GetNbins()+1; i++)  xbinsEta[i]=hRecoEta->GetXaxis()->GetBinLowEdge(i+1);
  xbinsEta[hRecoEta->GetXaxis()->GetNbins()]=hRecoEta->GetXaxis()->GetBinUpEdge(hRecoEta->GetXaxis()->GetNbins());

  hCwzEta    = new TH1D("hCwzEta"   ,"Cw",nxEta,xbinsEta);
  hPtywzEta  = new TH1D("hPtywzEta" ,"Pu",nxEta,xbinsEta);
  hStabwzEta = new TH1D("hStabwzEta","St",nxEta,xbinsEta);

  TH1D *hCwzErrors, *hCwzErrorsNew;
  hCwzErrors= new TH1D("hCwzErrors","old errors",nxEta,xbinsEta);
  hCwzErrorsNew= new TH1D("hCwzErrorsNew","new errors",nxEta,xbinsEta);

  double sum_recoEta=0., sum_truthEta=0.;
  for(int b=1; b< (2+(hCwzEta->GetNbinsX())); b++){
    
    double mcrecoEta=hRecoEta->GetBinContent(b);
    if(mcrecoEta==0) continue;
    double mctruthEta=hTruthEta->GetBinContent(b);
    if(mctruthEta==0) continue;   
    
    double mcstayRecoEta = hEtaCstayReco->GetBinContent(b);
    if(mcstayRecoEta == 0) continue;
    double mcstayTruthEta = hEtaCstayTruth->GetBinContent(b);
    if(mcstayTruthEta == 0) continue;
    
    sum_recoEta+=mcrecoEta;
    sum_truthEta+=mctruthEta;
    
    double errorReco = hRecoEta->GetBinError(b);
    double errorTruth = hTruthEta->GetBinError(b);

    double errorCw = std::sqrt( (1./std::pow(mctruthEta,2))*std::pow(errorReco,2) + (1./std::pow(mctruthEta,4))*std::pow(errorTruth,2)*std::pow(mcrecoEta,2) );

    hCwzEta->SetBinContent(b, mcrecoEta/mctruthEta);
    hCwzEta->SetBinError(b,errorCw);
    hCwzErrors->SetBinContent(b,errorCw);
    
    hPtywzEta->SetBinContent(b, mcstayTruthEta/mcrecoEta);
    hPtywzEta->SetBinError(b,0);
    
    hStabwzEta->SetBinContent(b, mcstayTruthEta/mctruthEta);
    hStabwzEta->SetBinError(b,0);
  }

  fout->cd();
  hCwzErrors->Write("hCwzErrors");

  std::cout<<" Cw:   "<<sum_recoEta<<"    "<<sum_truthEta<<"     "<<sum_recoEta/sum_truthEta<<std::endl;

  for (int iEta = 0; iEta< nxEta + 1; iEta++){
    
    double xlow = xbinsEta[iEta-1];
    double xhi  = xbinsEta[iEta];

    double C    = hCwzEta->GetBinContent(iEta);
    double dC   = hCwzEta->GetBinError(iEta);

    if(iEta == 0) C = hRecoEta->Integral() / hTruthEta->Integral();
    if(iEta == 0){
      double errorReco=0., errorTruth=0.;
      double mcreco  = hRecoEta->IntegralAndError(0,hRecoEta->GetNbinsX()+1,errorReco,"");
      double mctruth = hTruthEta ->IntegralAndError(0,hTruthEta->GetNbinsX()+1,errorTruth,"");
      dC  = std::sqrt( (1./std::pow(mctruth,2))*std::pow(errorReco,2) + (1./std::pow(mctruth,4))*std::pow(errorTruth,2)*std::pow(mcreco,2) );
    }

    double Nrec = hEtaCstayReco->GetBinContent(iEta)  + hEtaCcomeReco->GetBinContent(iEta);
    double Ngen = hEtaCstayTruth->GetBinContent(iEta) + hEtaCleaveTruth->GetBinContent(iEta);

    if(iEta == 0) Nrec = hEtaCstayReco->Integral()  + hEtaCcomeReco->Integral();  


    //check counting is correct!
    double NrecSimple = hRecoEta->GetBinContent(iEta);
    double NgenSimple = hTruthEta->GetBinContent(iEta);
    if(((NrecSimple-Nrec)/NrecSimple)*100. > 0.0001) std::cout<<"Inconsistency in counting come,stay,leave N"<<std::endl;
    if(((NgenSimple-Ngen)/NgenSimple)*100. > 0.0001) std::cout<<"Inconsistency in counting come,stay,leave N"<<std::endl;

    //errors on number of events in stay, come , leave categories
    double dStayR2  = std::pow(hEtaCstayReco->GetBinError(iEta), 2);
    double dCome2   = std::pow(hEtaCcomeReco->GetBinError(iEta), 2);
    double dStayT2  = std::pow(hEtaCstayTruth->GetBinError(iEta), 2);
    double dLeave2  = std::pow(hEtaCleaveTruth->GetBinError(iEta), 2);

    if(iEta==0){
      double dummy1=0., dummy2=0.,dummy3=0.;
      double int1 = hEtaCstayTruth ->IntegralAndError(0,hEtaCstayTruth ->GetNbinsX()+1, dummy1,""); dStayT2 = std::pow(dummy1, 2);
      double int2 = hEtaCcomeReco  ->IntegralAndError(0,hEtaCcomeReco  ->GetNbinsX()+1, dummy2,""); dCome2  = std::pow(dummy2, 2);
      double int3 = hEtaCleaveTruth->IntegralAndError(0,hEtaCleaveTruth->GetNbinsX()+1, dummy3,""); dLeave2 = std::pow(dummy3, 2);
    }

    //contributions to delta_C
    double dCstayR2 = std::pow( (1.0/Ngen), 2);
    double dCcome2  = std::pow( (1.0/Ngen), 2);
    double dCstayT2 = std::pow( Nrec/(Ngen*Ngen), 2);
    double dCleave2 = std::pow( Nrec/(Ngen*Ngen), 2);

    dCstayT2 = std::pow( (Ngen-Nrec)/(Ngen*Ngen), 2);//robs formula (un)comment

    dCstayR2 = (dCstayR2* dStayR2);
    dCcome2  = (dCcome2 * dCome2) ;
    dCstayT2 = (dCstayT2* dStayT2);
    dCleave2 = (dCleave2* dLeave2);

    double dCnew = std::sqrt(dCstayR2 + dCcome2 + dCstayT2 + dCleave2);
    dCnew = sqrt(dCstayT2 + dCcome2 + dCleave2);//robs formula (un)comment

    hCwzErrorsNew->SetBinContent(iEta,dCnew);
    hCwzEta->SetBinError(iEta,dCnew);
    printf("i=%2i  M=[%4.0f-%4.0f]  C=%6.4f  dC=%6.4f  dCnew=%6.4f  Nrec=%12.2f   Ngen=%12.2f\n",iEta,xlow,xhi,C,dC,dCnew,Nrec,Ngen);

    //for tables
    if(iEta!=0){
      Nrec_temp.push_back(hEtaCstayReco->GetBinContent(iEta)  + hEtaCcomeReco->GetBinContent(iEta));
      Ngen_temp.push_back(hEtaCstayTruth->GetBinContent(iEta) + hEtaCleaveTruth->GetBinContent(iEta));
      NstayRec_temp.push_back(hEtaCstayReco->GetBinContent(iEta));
      NstayGen_temp.push_back(hEtaCstayTruth->GetBinContent(iEta));
      NleaveGen_temp.push_back(hEtaCleaveTruth->GetBinContent(iEta));
      NComeRec_temp.push_back(hEtaCcomeReco->GetBinContent(iEta));
    }else if(iEta==0){
      Nrec_temp.push_back(hEtaCstayReco->Integral()  + hEtaCcomeReco->Integral());
      Ngen_temp.push_back(hEtaCstayTruth->Integral() + hEtaCleaveTruth->Integral());
      NstayRec_temp.push_back(hEtaCstayReco->Integral());
      NstayGen_temp.push_back(hEtaCstayTruth->Integral());
      NleaveGen_temp.push_back(hEtaCleaveTruth->Integral());
      NComeRec_temp.push_back(hEtaCcomeReco->Integral());
    }

    dCwOld_temp.push_back(dC);
    dCwNew_temp.push_back(dCnew);
    Cw_temp.push_back(C);
  }
  
  fout->cd();
  hCwzErrorsNew->Write("hCwzErrorsNew");
  
  //latex tables
  std::vector<std::string> name_eta_bin = LatexBinName();
  LatexTable("cw_factor",name_eta_bin,config.WZSelection);
  
  return;
}


void CwCalculator::LatexTable(std::string fileName, std::vector<std::string> eta_bin_name,
			      std::string channel){
  std::string labelboson;
  
  if(channel=="wminus")
    labelboson="$W^{-}\\rightarrow\\mu^{-}\\nu$";
  else if (channel=="wplus")
    labelboson="$W^{+}\\rightarrow\\mu^{+}\\nu$";
  
  std::vector<std::string> header; header.clear();
  header.push_back("$N_{tot}^{reco}$");
  header.push_back("$N_{stay}^{reco}$");
  header.push_back("$N_{come}^{reco}$");
  header.push_back("$N_{tot}^{truth}$");
  header.push_back("$N_{stay}^{truth}$");
  header.push_back("$N_{leave}^{truth}$");
  header.push_back("$C_{W}$");
  header.push_back("$\\delta^{stat}$");
  header.push_back("$\\delta_{new}^{stat}$");

  int nHeader = (int)header.size(); int nEtaBins = (int)Nrec_temp.size();
  std::ofstream CwLatex(("./"+fileName+"_"+channel+".tex").c_str());

  CwLatex << "\\begin{table}[tb!]" << std::endl;
  CwLatex << "\\begin{center}" << std::endl;
  CwLatex << "\\begin{tabular}{ c |";
  for(int k=0; k < nHeader; k++){
    if(k>=0 && k<=nHeader-4)CwLatex<<" C{2cm} "; 
    else CwLatex<<" C{1cm} ";}
  CwLatex<<" }\\hline \\hline " << std::endl;

  //header row
  CwLatex<< "Eta bin & ";
  for(int k=0; k < nHeader; k++){
    CwLatex<<std::setw(4)<<std::right;
    if(k< nHeader-1) CwLatex<<" "<<header.at(k)<<" &";
    if(k==nHeader-1) CwLatex<<" "<<header.at(k)<<" \\\\ \\hline"<<std::endl;
  }

  double totdC=0, totdCnew=0.;
  for(int j=1; j < nEtaBins; j++){
    totdC=totdC+std::pow(dCwOld_temp.at(j),2);
    totdCnew=totdCnew+std::pow(dCwNew_temp.at(j),2);
  }

  for(int j=0; j < nEtaBins; j++){
    double dC=0., dCnew=0.;
    if(j==0){
      dC    = std::sqrt(totdC);
      dCnew = std::sqrt(totdCnew);
    }else{
      dC = dCwOld_temp.at(j);
      dCnew = dCwNew_temp.at(j);
    }
    //CwLatex<<std::setw(6)<<std::right;
    CwLatex<<"  "<<eta_bin_name.at(j)         <<" &";
    CwLatex<<"  \\num{"<<Nrec_temp.at(j)      <<"} &";  /*std::fixed<<std::setprecision(2)<<*/
    CwLatex<<"  \\num{"<<NstayRec_temp.at(j)  <<"} &";
    CwLatex<<"  \\num{"<<NComeRec_temp.at(j)  <<"} &"; 
    CwLatex<<"  \\num{"<<Ngen_temp.at(j)      <<"} &";
    CwLatex<<"  \\num{"<<NstayGen_temp.at(j)  <<"} &";
    CwLatex<<"  \\num{"<<NleaveGen_temp.at(j) <<"} &";
    CwLatex<<"  \\num{"<<Cw_temp.at(j)        <<"} &";  
    CwLatex<<"  \\num{"<<dC*100               <<"} &";
    CwLatex<<"  \\num{"<<dCnew*100<<"}  \\\\ \\hline"<<std::endl;
  }
  
  CwLatex << "\\hline" << std::endl;
  CwLatex << "\\end{tabular}" << std::endl;
  CwLatex << "\\end{center}"  << std::endl;
  CwLatex << "\\caption{This table shows the elements needed to calculate the MC statistical uncertainty defined in \\cref{eq:cwerror} for the "<<labelboson<<" selection. The first column is the muon pseudorapidity bins; the second, third and fourth column are the total, stay and total reconstructed number of events, respectively; the fifth, sixth and seven are total, stay and leave truth number of events, respectively; the seventh column is the unfolding correction factor; the eighth colum is the statistical uncertainty calculated for considering the reconstruction and truth total number of events independent; finally, the right last column is the correctly calculated statistical uncertainties, in percentage.}" << std::endl;  
  CwLatex << "\\label{tab:Cwfac"<<channel<<"}" << std::endl;  
  CwLatex << "\\end{table}"<< std::endl;
  CwLatex.close();

  return;

}

std::vector<std::string> CwCalculator::LatexBinName(){

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

void CwCalculator::finalize(){

  std::cout<<"**** Finalizing Cw calculation ****"<<std::endl;
  fout->cd();

  hReco->Write("hMassCrec");
  hTruth->Write("hMassCgen");
  hMwtCstayTruth->Write("hMassCstayG");
  hMwtCleaveTruth->Write("hMassCleave");
  hMwtCstayReco->Write("hMassCsta");
  hMwtCcomeReco->Write("hMassCcome");
  hCwz->Write("hCwz");
  hPtywz  ->Write("hPtywz");
  hStabwz ->Write("hStabwz");

  hRecoEta->Write("hEtaCrec");
  hTruthEta->Write("hEtaCgen");
  hEtaCstayTruth->Write("hEtaCstayG");
  hEtaCleaveTruth->Write("hEtaCleave");
  hEtaCstayReco->Write("hEtaCsta");
  hEtaCcomeReco->Write("hEtaCcome");
  hCwzEta->Write("hCwzEta");
  hPtywzEta  ->Write("hPtywzEta");
  hStabwzEta ->Write("hStabwzEta");

  hCw_simple->Write("hCw_simple");
  hCw_eta_simple->Write("hCw_eta_simple");
  
  fout->Close();

  return;
}
#endif
