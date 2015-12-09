# A search for dark matter produced in association with a Higgs boson

### Introduction
This package implements an analysis of ATLAS Experiment data designed to look
for Higgs bosons produced in association with dark matter particles. The Higgs
decay is identified by a diphoton resonance, while the dark matter particle
would manifest as missing transverse energy in the detector.

The code has been structured so that all of the general analysis settings are 
stored in the configuration file in the data/ directory. *A typical user should 
only need to adjust the config file. Changes to cutflows are temporary 
exceptions, and should be added to the DMEvtSelect class. The code is designed 
to be as automatic as possible. Each class looks to see if all of the necessary 
inputs have been produced previously before generating them from scratch.

### General analysis strategy:
1)  mini-MxAODs from xAODs using the tools provided by Hgamma WG.
2)  mass points for data, backgrounds, and signal.
3)  parameterization of the SM and DM signals.
4)  background modeling.
5)  workspace to store models and PDFs.
6)  pseudoexperiment ensemble generation and analysis 
7)  CLs and p0 calculators

### Package contents:

##### settingsHDM_sys.cfg
  Luminosity, higgs mass, m_yy range, file names, script locations, production 
  mode information should all go here. Default job options also are included. 
  The sys tag refers to the fact that MxAODs with systematic variations can now
  be specified. Be careful: there are differences between the nominal and 
  systematic MxAOD path specifications. Specifically, one must specify the 
  entire file location for systematics MxAODs (root://eosatlas//eos...), since 
  the files are spread across many locations.

##### DMAnalysis
  This namespace should store all general analysis methods. The idea is to
  avoid duplication of methods to minimize the pain of changing things. 

##### DMMaster
  This is the master 'wrapper' class for the analysis. Using this class, all the
  analysis tools can be run. The file organization is automated, using a 
  directory structure based on the 'masterInput' and 'masterOutput' strings in 
  settingsHDM, as well as the masterJobName.

##### DMMassPoints
 This program uses a TTree of data events to produce a series of mass points
 that can be used as inputs for the signal parameterization or workspace 
 creation. The cutflow is implemented using the DMEvtSelect class.
  
##### DMSigParam
 This program uses signal MC to fit the resonance shape for the SM Higgs and
 the DM signal and saves the parameters for use in workspace generation. The
 fit is performed on masspoints generated with DMMassPoints. Signal cross-
 sections are provided by the BRXSReader tool.

##### BkgModel
 This program implements all of the possible background models, and can return 
 either a RooAbsPdf object, a CombinedPdf, or add a PDF directly to the analysis
 workspace. It is designed to be generic enough for use in all H->yy analyses. 

##### DMWorkspace
 This program produces the statistical model for the DM analysis. It includes SM
 signal, a single DM signal, and background PDFs, as well as associated 
 systematic uncertainties. The parameter of interest is "mu_DM", the signal 
 strength for the dark matter production process. The signal strength of the 
 Standard Model Higgs Boson "mu_SM" is set to 1. The background normalization 
 comes from data. There is the option of fitting the SM signal strengths 
 individually. 

##### DMTestStat
 This program calculates the 95% CL, CLs, and p0 values for a given DM signal.
 It is the default fitting program. 

##### DMPseudoExp
 This program implements the pseudo-dataset generation and fitting. It is 
 designed to either run locally or on a cluster, and outputs a TTree file.

##### DMToyAnalysis
 This program has tools for analyzing toy MC data. 

### Supporting Classes:

##### BRXSReader
 Reads tables of SM Higgs cross sections and branching ratios and provides an 
 easy-to-use interface.

##### DMCheckJobs
 Checks to see whether the output files for a given program are available. It
 also assembles a list of failed jobs that can be prepared for resubmission.

##### DMEvtSelect
 This class implements the cutflow and counters for the analysis. It can be 
 initialized using a pointer to the DMTree. 

##### DMTree
 This class is automatically generated based on the MxAOD structure. It provides
 a useful interface for code to the TTrees. 

##### PESReader
 This is a simple class for loading and accessing energy scale systematics based
 on an input file. 

##### PERReader
 This is a simple class for loading and accessing energy resolution systematics
 based on an input file. 

### Setting up the package. 

##### User modifications
New users will need to modify data/settingsHDM.cfg and possibly the makefile in 
order to run the code. All input and output file locations for the programs
are built around the locations provided in the settings (masterInput, 
masterOutput, packageLocation, clusterFileLocation, and fileName*). Sub-
directories will be created as necessary by the program.

##### Input files.
To be uploaded shortly. For the moment, they are located in the following lxplus
directory: ~ahard/public/GlobalInputs. The masterInput directory should point to
these files.

### Running the code
First compile the master program, which executes the analysis code:
     > make bin/DMMaster

Then to run,

     ./bin/DMMaster <Program> <SettingsFile>

The program can be any of the following options: 
  - Cleanup (clean old files from previous analysis runs)
  - MassPoints (make mass files as inputs for and model)
  - GetSystematics (make cutflows and categorizations for all syst. variations)
  - PlotVariables (plot interesting kinematic variables)
  - SigParam (build the signal PDF from MC)
  - BkgModel (build the background model)
  - Workspace (build the statistical model)
  - ResubmitWorkspace (submit failed Workspace jobs again)
  - TossPseudoExp (toss pseudo experiment ensemble)
  - PlotPseudoExp (plot the results of pseudo experiments)
  - TestStat (calculate p0 and CLs)
  - ResubmitTestStat (submit failed TestStat jobs again)
  - MuLimit (get the 95% CL limit on the parameter of interest)
  - Optimizer (optimize the analysis selection with meta job)
  - OptAnalysis (analyze the results of Optimizer)

The code will automatically run any required upstream programs in order to 
ensure that it has all required inputs. For instance, if you want to create a 
workspace from scratch, just run,

     ./bin/DMMaster Workspace data/settingsHDM.cfg

Make sure that you are running in a directory from which EOS is accesssible. 

At the moment, the bash scripts for grid or lxbatch submission in scripts/ have
lots of hard-coded user-specific information. These will be udpated in the 
future. In any case, they do not influence local jobs in any way.