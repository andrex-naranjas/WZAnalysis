//File add
#ifndef FILEADD_H
#define FILEADD_H

#include "ConfigSettings.h"

#include <TROOT.h>
#include <TChain.h>
#include "TIterator.h"
#include "TKey.h"
#include "TObject.h"
#include "TCollection.h"
#include "TString.h"

#include "TList.h"
#include "TDirectory.h"

#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"


class FileAdd{

 public:
  FileAdd();
  virtual ~FileAdd();
  virtual void execute(Config config, std::string samples, bool recoil);

 private:
  void MergeRootfile( TDirectory *target, TList *sourcelist );
  TList *FileList;
  TFile *Target;

};

#endif //> !FILEADD_H

