//HaddClass
#ifndef FILEADD_CXX
#define FILEADD_CXX

#include "FileAdd.h"

#include <iostream>
#include <vector>
#include <string>


FileAdd::FileAdd()
{
}

FileAdd::~FileAdd(){}

void FileAdd::execute(Config config, std::string sample, bool recoil){

  //names for output files
  std::string wzchannel;
  if(config.WZSelection=="zmumu"){
    wzchannel="z";
  }else if(config.WZSelection=="wplus"){
    wzchannel="wplus";
  }else if(config.WZSelection=="wminus"){
    wzchannel="wminus";}
  std::string calib;
  if(config.SETCalibration=="True" && config.InsituCorrection=="True"){
    calib="_set_insitu";
  }else if(config.SETCalibration=="True" && config.InsituCorrection=="False"){
    calib="_set";
    }else if(config.SETCalibration=="False" && config.InsituCorrection=="True"){
    calib="_insitu";
  }else if(config.SETCalibration=="False" && config.InsituCorrection=="False"){
    calib="";
  }
  std::string puname="";
  if(config.OnTheFlyPileUp=="False") puname="_nonewPU";

  std::string  dir, prefix;
  if(!recoil){dir="Files/"; prefix="";}
  if(recoil) {dir="RecoilHistos/"; prefix="RecoilCalibHistos_";}

  //ouput directory according year
  std::string dirYear="";
  if(config.DataYears=="2015+2016") dirYear="2015p2016/";
  if(config.DataYears=="2017") dirYear="2017/";
  if(config.DataYears=="2018") dirYear="2018/";

  std::vector<std::string> systematic; systematic.clear();
  systematic.push_back("_nominal");

  if(config.Systematics=="All"){
    systematic.push_back("_MUONS_ID__1down");
    systematic.push_back("_MUONS_ID__1up");
    systematic.push_back("_MUONS_MS__1down");
    systematic.push_back("_MUONS_MS__1up");
    systematic.push_back("_MUONS_SCALE__1down");
    systematic.push_back("_MUONS_SCALE__1up");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_1__1down");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_1__1up");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_2__1down");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_2__1up");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_3__1down");
    systematic.push_back("_JET_NPScenario1_JET_GroupedNP_3__1up");
    systematic.push_back("_JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1down");
    systematic.push_back("_JET_NPScenario1_JET_EtaIntercalibration_NonClosure__1up");
    systematic.push_back("_MET_SoftTrk_ScaleDown");  
    systematic.push_back("_MET_SoftTrk_ScaleUp");
    systematic.push_back("_JET_JER_SINGLE_NP__1up");
    systematic.push_back("_MET_SoftTrk_ResoPerp");
    systematic.push_back("_MET_SoftTrk_ResoPara");
  }
  if(config.SFVariations=="True"){
    systematic.push_back("_IDStatD");
    systematic.push_back("_IDSysD");
    systematic.push_back("_IsoStatD");
    systematic.push_back("_IsoSysD");
    systematic.push_back("_TrigStatD");
    systematic.push_back("_TrigSysD");
    systematic.push_back("_IDStatU");
    systematic.push_back("_IDSysU");
    systematic.push_back("_IsoStatU");
    systematic.push_back("_IsoSysU");
    systematic.push_back("_TrigStatU");
    systematic.push_back("_TrigSysU");
  }

  for(int iSys=0; iSys<(int)systematic.size(); iSys++){
    
    if(recoil) systematic[iSys]="";

    Target = TFile::Open((config.OutputFileDir + dir + dirYear +"Add/"+ prefix + sample + "_"+ wzchannel + systematic[iSys] + calib + puname +"_Total.root").c_str(),"RECREATE");
    
    std::vector<std::string> massrange; massrange.clear();
    massrange.push_back("");
    massrange.push_back("_120");
    massrange.push_back("_180");
    massrange.push_back("_250");
    massrange.push_back("_400");
    massrange.push_back("_600");
    massrange.push_back("_800");
    massrange.push_back("_1000");
    massrange.push_back("_1250");
    massrange.push_back("_1500");
    massrange.push_back("_1750");
    massrange.push_back("_2000");
    massrange.push_back("_2250");
    massrange.push_back("_2500");
    massrange.push_back("_2750");
    massrange.push_back("_3000");
    massrange.push_back("_3500");
    massrange.push_back("_4000");
    massrange.push_back("_4500");
    massrange.push_back("_5000");

    FileList = new TList();
    
    for(int iRangemass = 0; iRangemass < (int)massrange.size(); ++iRangemass){
      FileList->Add( TFile::Open((config.OutputFileDir + dir + dirYear + prefix+ sample + massrange[iRangemass] + "_"+ wzchannel + systematic[iSys] + calib + puname +".root").c_str()) );
      std::cout<<(config.OutputFileDir + dir + dirYear + prefix+ sample + massrange[iRangemass] + "_"+ wzchannel + systematic[iSys] + calib + puname +".root").c_str()<<"   Adding this file..."<<std::endl;
    }
    
    MergeRootfile( Target, FileList );
    delete FileList;
    
    if(recoil) break;

  }//systematic loop
  
  return;
}


void FileAdd::MergeRootfile(TDirectory *target, TList *sourcelist){

   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );

   TFile *first_source = (TFile*)sourcelist->First();
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();

      if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         TH1 *h1 = (TH1*)obj;

         // loop over all source files and add the content of the
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         while ( nextsource ) {

            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
            if (key2) {
               TH1 *h2 = (TH1*)key2->ReadObj();
               h1->Add( h2 );
               delete h2;
            }

            nextsource = (TFile*)sourcelist->After( nextsource );
         }
      }
      else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

         // loop over all source files create a chain of Trees "globChain"
         const char* obj_name= obj->GetName();

         globChain = new TChain(obj_name);
         globChain->Add(first_source->GetName());
         TFile *nextsource = (TFile*)sourcelist->After( first_source );

         while ( nextsource ) {

            globChain->Add(nextsource->GetName());
            nextsource = (TFile*)sourcelist->After( nextsource );
         }

      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory

	std::cout << "Found subdirectory " << obj->GetName() << std::endl;

         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

         MergeRootfile( newdir, sourcelist );

      } else {

         // object is of no type that we know or can handle
	std::cout << "Unknown object type, name: "
		  << obj->GetName() << " title: " << obj->GetTitle() << std::endl;
      }

      if ( obj ) {
         target->cd();

         //!!if the object is a tree, it is stored in globChain...
         if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
         else
            obj->Write( key->GetName() );
      }

   } // while ( ( TKey *key = (TKey*)nextkey() ) )

   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);

   return;
}

#endif
