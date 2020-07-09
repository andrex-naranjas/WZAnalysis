//simple class to store reco info
//autor, Andres Ramirez: andres.ramirez.morales@cern.ch
#ifndef BOSONRECO_H
#define BOSONRECO_H

class BosonReco{

 public:
  BosonReco(){};
  ~BosonReco(){};

  ULong64_t eventNumberTruth;	    
  UInt_t    runNumberTruth;	    
  UInt_t    mcChannelNumber;	    
	                            
  Float_t   weight_mc;		    
  Float_t   weight_pileup;	    
  Double_t  KFactor_weight;
  
  double mu_pt;
  double mu_eta;
  double mu_phi;
  double mu_charge;

  /* std::vector<float>  *mu_pt; */
  /* std::vector<float>  *mu_eta; */
  /* std::vector<float>  *mu_phi; */
  /* std::vector<float>  *mu_charge; */

  Float_t met_met;
  Float_t met_phi;

};

#endif
