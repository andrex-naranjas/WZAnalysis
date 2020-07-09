This is a guide of how to use the WZOffLineAnalysis Framework<br />
Author:<br />
Andres Ramirez-Morales (andres.ramirez.morales@cern.ch)<br />

1.) Setup the local enviroment:<br />
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase <br />
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh' <br />
export ALRB_rootVersion=6.14.04-x86_64-slc6-gcc62-opt  <br />
setupATLAS <br />
lsetup git <br />
lsetup root <br />

2.) Setup a afs/krb5 token <br />
export CERN_USER=insert_your_lxplus_login_here<br />
kinit -5 $CERN_USER@CERN.CH<br />
aklog -c cern.ch -k CERN.CH<br />
<br />
Note: if this does not work for you, please use your ssh keys<br />
and clone the repository with ssh <br />

3.) Initial package setup<br />
Pick a working directory somewhere you like<br />
git clone https://:@gitlab.cern.ch:8443/anramire/WZOffLineAnalysis.git
cd WZOffLineAnalysis<br />
mkdir obj CheckPlots<br />
make<br />
<br />
Note: to clean and recompile from scratch please do:<br />
make clean<br />
make obj<br />
make<br />

4.) Create the directories where your files,tables and plots will be stored <br />
Pick/create a working directory somewhere you like <br />
cp createDirectories.sh yourAwesomeDirectory <br />
cd yourAwesomeDirectory <br />
source createDirectories.sh <br />


5.) Run out of the box<br />
When the code has successfully compiled,please setup the config file (below you will find an explanation of it), and do:<br />
./WZAnalysis r21config.txt

After a succesfully run (after typing make), you can check the the output of your studies organised in different directories

6.) Input and configuration<br />
This framework is intended to be used with the outuput ntuples from AnalysisTop framework, i.e. flat ntuples with trees: nominal, systematics and truth

This Framework is highly configurable and you can perform several studies, starting as simple as running event by event and produce plots for data, to calculate the xSection unfolding factors.

The input should contain all the information to perform the several studies, please check the configuration file which enables or disables functionalaties as required. <br />

Configuration file (r21config.txt or create your own!):
This configuration file should contain the relevant settings to perform the analysis, here is the explanation of each line:

     NormPlots = True   # Value: True or False. Usage: set True to display of normalised plots for the hadronic recoil variables
     DoPlots = True     # Value: True or False. Usage: set True to obtain the kinematic control plots, make sure you have run the analysis first
     DoOnlyPlots = True # Value: True or False. Usage: set False to run the looping of events. i.e. apply cuts, weights, produce histograms

     DataYears = 2017   # Value: 2015+2016 or 2017. Usage: Select the years which are to be analised

     InputFileDir = /data/morales/atlas/Nov2018/mc16d/   # Value: string indicating the path of where your samples are stored
     OutputFileDir = /data/morales/atlas/r21ControlPlots # Value: string indicating where your results will be stored, see 3.)

     NumberOfEvents = -1 # Value: -1 or any integer number. Usage: use this option to pick the number of events you would like to run (-1 == to all events)

     DoCwFactor = True     # Value: True or False. Usage: set True to calculate unfolding factors and its MC statistical errors.
     DoCrossSection = True # Value: True or False. Usage: set True to calculate cross sections and lepton charge asymmetries. Pro
     TruthAnalysis = False # Value: True or False. Usage: set True to run the truth event loop, this is needed to calculate the xSec and unfolding factors

     DoMultijet = False       # Value: True or False. Usage: set True to obtain the histograms for the multijet background in the event loop
     MultiFitRegion = Signal  # Value: Signal, FitR1 or FitR2. Usage: select the desired multijet region.
     MultijetVariation = Nominal # Value: Nominal...
     OnlyMultijetResults = True  # Value: True or False. Usage: set False to calculate the multijet normalisation and then produce results; set True to only produce the plot and tables.

     Dod0Cut = True # Value: True or False. Usage: set True to apply the d0sig cut in case you do not have (TTVA cut).
     Applyd0Weight = False # Value: True or False. Usage: set True to apply the d0sig reweight (for multijet calculation)
     Nod0Shift = False #Value: True or False. Usage: set True for not applying the d0sig MC mean gaussian mean shift (for multijet calculation)

     xBinsCw = 40000 45000 50000 55000 60000 65000 70000 200000 #Value: space-separated numbers. Usage: set the mass binning for measurements

     OnlyMC = True # Value: True or False. Usage: set True to run only Monte Carlo (no data)
     OnlyData = False # Value: True or False. Usage: set True to run only Data (no MC)
     OnlyInclusive = False # Value: True or False. Usage: set False to include the mass slices.

     WZSelection = wminus # Value: wminus, wplus, zmumu, combined. Usage: set the signature selection; combined option only works when doing only plots.

     WChannel = True  # Value: True or False. Usage: set True to include W boson Monte Carlo samples.
     ZChannel  = True # Value: True or False. Usage: set True to include Z boson Monte Carlo samples.
     TopChannel = True # Value: True or False. Usage: set True to include top quark Monte Carlo samples.
     DibosonChannel = True # Value: True or False. Usage: set True to include diboson Monte Carlo samples.

     MuFlavor = True # Value: True or False. Usage: set True to include muon Monte Carlo samples.
     TauFlavor = True # Value: True or False. Usage: set True to include tau Monte Carlo samples.

     CutFlowName = # Value: string indicating the path of where your cutflows are stored

     HasRecoilInfo = False # Value: True or False. Usage: set True to indicate that your input ntuples have the recoil variables stored.

     OnTheFlyPileUp = False # Value: True or False. Usage: set True to perform the on-the-fly pileup reweighting.

     SETCalibration = False # Value: True or False. Usage: set True to perform the SumET calibration (hadronic recoil).
     RecoilCalibFileDir = /data/morales/atlas/HadRecoilCalib/FullCalibration_15p16.root # Value: string of the directory where the calibration files are stored.
     PRWFile = /users/morales/ATLAS/ControlPlots2018/Plots/PRWFiles/prwTree_zmumu_opt.root # Value: string of the directory where the pileup weighting files are stored.

     InsituCorrection = False # Value: True or False. Usage: set True to perform the in-situ corrections (hadronic recoil).
     ResolResponse = False # Value: True or False. Usage: set True to perform the resolution and response corrections (hadronic recoil)

     TruthMatching = True # Value: True or False. Usage: set True to perform the truth matching when running the reconstrunction event loop.
     RecoMatching = True  # Value: True or False. Usage: set True to perform the reconstruction matching when running the truth event loop.

     Systematics = False # Value: True or False. Usage: set to True to enable the calculation of systematic variations.
     SFVariations = False # Value: True or False. Usage: set to True to enable the calculation of the scale factor variations.

     RecoilBG = False  # Value: True or False. Usage: set to True to enable the hadronic recoil backgrounds.

     SysInplots = 1     # Value: 1 or 0 (integers). Usage: set to 1 to include the systematic bands in the control plots.
     PlotMulti = 1      # Value: 1 or 0 (integers). Usage: set to 1 to include the multijet background background in the control plots.
     PlotData = 1       # Value: 1 or 0 (integers). Usage: set to 1 to include the data in the control plots.
     PlotZmumu = 1      # Value: 1 or 0 (integers). Usage: set to 1 to include the Zmumu boson background in the control plots.
     PlotZtautau = 1    # Value: 1 or 0 (integers). Usage: set to 1 to include the Ztautau in the control plots.
     PlotWplusTaunu = 1 # Value: 1 or 0 (integers). Usage: set to 1 to include the Ztautau in the control plots.
     PlotWplusMunu = 1
     PlotWminTaunu = 1
     PlotWminMunu = 1
     PlotTtbar = 1
     PlotWt_top = 1
     PlotWt_antitop = 1
     PlotSingle_top = 1
     PlotSingle_antitop = 0
     PlotZzqqll = 1
     PlotWwqqll = 1
     PlotWwpqqmlnu = 1
     PlotWwplnumqq = 1
     PlotWzlnuqq = 1
     PlotZzllll = 1
     PlotWzlnull = 1
     PlotWzlnununu = 1
     PlotZzllnunu = 1

The structure of the package is organised using several c++ classes to perform different tasks, perhaps the most important one, where the cuts and the calling of the files is MyWZAnalysis, here the cuts and correction scale factors are applied, this class should be modified to take into account major changes in the analysis. Please contact the author for further information.

