//
//   @file    AtlasLabels.h         
//   
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasLabels.h 243182 2015-08-02 09:51:02Z vincter $


#ifndef __ATLASLABELS_H
#define __ATLASLABELS_H

#include "Rtypes.h"

void ATLASLabel(Double_t x,Double_t y,const char* text=NULL,Color_t color=1); 

void ATLASLabelOld(Double_t x,Double_t y,bool Preliminary=false,Color_t color=1); 

void ATLASVersion(const char* version=NULL,Double_t x=0.88,Double_t y=0.975,Color_t color=1); 

#endif // __ATLASLABELS_H
