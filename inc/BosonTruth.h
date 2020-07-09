//simple class to store truth info
//autor, Andres Ramirez: andres.ramirez.morales@cern.ch
#ifndef BOSONTRUTH_H
#define BOSONTRUTH_H

class BosonTruth{

 public:
  BosonTruth(){};
  ~BosonTruth(){};
  ULong64_t eventNumberTruth;	    
  UInt_t    runNumberTruth;	    
  UInt_t    mcChannelNumber;	    
	                            
  Float_t   weight_mc;		    
  Float_t   weight_pileup;	    
  Double_t  KFactor_weight_truth;   
	                            
  Float_t   MC_WZ_m;		    
  Float_t   MC_WZ_pt;		    
  Float_t   MC_WZ_eta;		    
  Float_t   MC_WZ_phi;		    
	                            
  Float_t   MC_lepton_m;	    
  Float_t   MC_lepton_pt;	    
  Float_t   MC_lepton_eta;	    
  Float_t   MC_lepton_phi;	    
	                            
  Float_t   MC_antilepton_m;	    
  Float_t   MC_antilepton_pt;	    
  Float_t   MC_antilepton_eta;	    
  Float_t   MC_antilepton_phi;             
};

#endif
