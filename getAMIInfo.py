#!/usr/bin/python

import os, glob, sys
from time import strftime
import subprocess
from shutil import copyfile

from ConfigParser import ConfigParser
from ConfigParser import Error as ConfigError

class AMI:
    """
    Class to get AMI cross sections.  
    Run-II version
    """
    
    def __init__(self, **kwargs):
        "Setup options"

        # write cross section or number of events        
        self.mode = 'crossSec' # 'Nevts'
        self.outFile='dslist_crossSec.txt'

        self.configFiles    = []
        
        # Overwrite attributes explicitly set
        for k, v in kwargs.items():    setattr(self, k, v)       
        
        print  self.configFiles 

        pass
    
    
    def query(self):
        "Loop over inputs and print cross sections, efficiencies or number of events"

        # output file
        outFile = open(self.outFile, 'w')

        print 'mode', self.mode


        for filelist in self.configFiles:            
            print 'filelist',filelist
            file = open(filelist)            
            for inDS  in iter(file):
                # look out for commented or short filelists
                if '#' in inDS:
                    continue
                if len(inDS)<10:
                    continue
                #print 'inDS',inDS
                #
                pycmd='ami show dataset info '
                pycmd += ' %s' % inDS
                pycmd=pycmd.strip()                
                dsid = inDS.split('.')[1]
                #print 'dataset ',dsid
                name=inDS.split('.')[2]
                isGen=False
                if 'EVNT' in inDS:
                    isGen=True
                    pass
                isData=False
                if 'data15_13TeV' in inDS or 'data16_13TeV' in inDS or 'data17_13TeV' in inDS :
                    isData=True
                    pass
                # choose evgen/daod strings
                crossSec_str="crossSection_mean"
                filtEff_str="GenFiltEff_mean"
                totalEvents_str="totalEvents"
                DAOD=False
                if not isData and "DAOD" in inDS:
                    DAOD=True
                    crossSec_str="crossSection"    
                    filtEff_str="genFiltEff"
                    pass
                #
                #print 'pycmd', pycmd
                pycmd = subprocess.Popen(pycmd, shell=True, stdout=subprocess.PIPE)
                #
                successXS=False
                successEff=False 
                successNEvt=False                       
                for line in pycmd.stdout:
                    #print 'line',line
                    #
                    if not isData and crossSec_str in line:  
                        if DAOD:
                            if line.find("approx") < 0 : 
                                successXS=True
                                crossSec=line 
                        else:
                            successXS=True
                            crossSec=line
                    # 
                    if not isData and filtEff_str in line:  
                        successEff=True
                        filtEff=line
                        pass
                    #
                    if totalEvents_str in line:
                        successNEvt=True
                        nEvt=line                                     
                        pass
                        #
                    pass
                # end of loop over results
                #print 'success',successXS
                if not isData and not successXS:
                    print 'error for %s'  % inDS
                    crossSec='error'
                    pass
                if not isData and not successEff:
                    print 'Gen Eff error for %s'  % inDS
                    filtEff='error'
                    pass                        
                if not successNEvt:
                    print 'totalEvents error for %s'  % inDS
                    nEvt='0'
                    pass
                #
                if not isData:
                    crossSec=crossSec.replace("crossSection_mean    : ","")
                    crossSec=crossSec.replace("crossSection         : ","")
                    crossSec=crossSec.strip()
                    filtEff=filtEff.replace("GenFiltEff_mean","")
                    filtEff=filtEff.replace("genFiltEff","")
                    filtEff=filtEff.replace(":","")
                    filtEff=filtEff.strip()        
                #        
                nEvt=nEvt.replace("totalEvents","")
                nEvt=nEvt.replace(":","")
                nEvt=nEvt.strip()
                #
                #  print 'crossSec', crossSec
                #  print 'filtEff', filtEff
                #  print 'nEvt', nEvt
                if self.mode=='crossSec' and not isData:
                    #                    outstring='%s %s\n'% (dsid, crossSec)
                    outstring='%s %s %s\n'% (dsid, crossSec, filtEff)
                    print outstring
                else:
                    # outstring='%s %s\n'% (dsid, nEvt)
                    outstring='%s %s %s\n'% (dsid, nEvt, name)
                    print outstring
                    pass
                outFile.write(outstring)
                pass
            pass
        outFile.close()
    

    def treatDuplicatedDSids(self):
        "Special treatment of any duplicated dataset ids"

        # CxAOD30 fix
        # add an mc16c to all MCs and mc16c_extension to the second occurance 
        # of any duplicates

        # copy outfile before overwriting
        tmpFile=self.outFile+'_tmp'
        copyfile(self.outFile, tmpFile)
        inFile = open(tmpFile, 'r')
        outFile = open(self.outFile, 'w')

        linesA=open(tmpFile, 'r').readlines()
        linesB=open(tmpFile, 'r').readlines()
        for lineA in linesA:
            dsidA = lineA.split(' ')[0]
            nevtA = lineA.split(' ')[1]
            nameA = lineA.split(' ')[2]
            nameA = nameA.strip()
            firstDup=False
            secondDup=False
            nevtSaved='-1'
            isData=False
            if 'physics_Main' in nameA:
                isData=True
            if not isData:    
                for lineB in linesB:
                    dsidB = lineB.split(' ')[0]
                    nevtB = lineB.split(' ')[1]
                    if dsidA == dsidB: 
                        if firstDup:
                            secondDup=True
                        if not firstDup:
                            firstDup=True;
                            if nevtSaved == '-1': # compare events in first occurance
                                nevtSaved=nevtB
                    #print 'duplicate',dsidA,nevtA,nevtB,firstDup,secondDup,nevtSaved
            extra_tag="mc16c"
            if firstDup and secondDup:
                if nevtA != nevtSaved:
                    extra_tag="mc16c_extension"
            if isData:
                extra_tag=""
            outstring='%s %s %s %s'% (dsidA, nevtA, nameA, extra_tag)    
            print outstring
            outFile.write(outstring+"\n")


if __name__ == '__main__':
    total = len(sys.argv)
    # number of arguments plus 1     
    if total!=3:
        print "You need some arguments, will ABORT!"
        print "Usage: ",sys.argv[0]," myfile folder_out"
        print "Usage: ",sys.argv[0]," ./FrameworkSub/In/CxAOD28/list_sample_grid.mc15c_13TeV_25ns_signal_nominal_HIGG2D4.txt ForCxAODDownload"
        assert(False)
    # done if
        
    # inputs
    myfile=sys.argv[1]
    folder_out=sys.argv[2]

    # run
    filelist=[myfile]
    print 'filelist',filelist
    #  p=AMI(configFiles=[ datafile ]) # one file at a time
    p=AMI(configFiles=filelist)
    p.mode='Nevts' #default is to print cross section, else print number of events
    if p.mode=='Nevts':
        p.outFile=folder_out+'/dslist_NevtDxAOD.txt'
    else:
        p.outFile=folder_out+'/dslist_sigmaEffDxAOD.txt'
    p.query()
# in CxAOD30 we have mc16c datasets+extensions with the same dsid 
# this function takes the first as mc16 and the second as mc16c_extension (as done in CxAOD submission ) to be compatible with the extra_tags defined in count_Nentry_SumOfWeight.py
#    p.treatDuplicatedDSids()


