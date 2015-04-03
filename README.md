# HggDM
Higgs to diphoton + dark matter search 

This package implements an analysis of ATLAS Experiment data designed to look
for Higgs bosons produced in association with dark matter particles. The Higgs
decay is identified by a diphoton resonance, while the dark matter particle
would manifest as missing transverse energy.

General analysis strategy:
     1) mini-MxAODs from xAODs using the tools provided by Hgamma WG.
     2) mass points for data, backgrounds, and signal.
     3) parameterization of the SM and DM signals.
     4) background modeling.
     5) workspace to store models and PDFs.
     6) pseudoexperiment ensemble generation and analysis 
     7) CLs and p0 calculators
     
NOTES 03/04/2015:

      The "FromFile" option should be modified. First, the program looks for the
      desired file. If it exists, open. Otherwise, call the 'createNew' method.


NOTES 24/03/2015:

      Need to think more about handling communication between analysis classes.
      For example, should DMSigParam have its own instance of DMMassPoints, or
      should DMMaster pass a reference during initialization? Same issue will 
      arise for workspace generation.
      
      Another issue: have instance of DMMassPoints for each sample? Yes. So then
      what about DMSigParam? I am inclined to do all of the signals in one shot,
      just because it is less cumbersome. 

      
OPEN ISSUES 22/03/2015: 

     * DMMassPoints only compiles on lxplus (not OSX Yosemite...). Why?
               
     * DMSigParam - need to choose RooRealVar ranges for fits. Also, need to
       consider using modified mass variable to prevent errors in the matrix	
       used for the fit.


Analysis Header:

     * DMHeader.h
          This header file should store all important analysis information. The
	  idea is to avoid hard-coding anything in the supporting classes. File
	  names, analysis luminosity, category names, production modes, mass
	  ranges are just a few of the quantities that should be defined here.

Main Classes:
     
     * DMMaster.cxx
          This is the master class for the analysis. Using this class, all the
     	  analysis tools can be run. The file organization is automated, using
     	  a directory structure based on the 'master_input' and 'master_output'
     	  strings in DMHeader.h, as well as the MasterJobName.

     * DMMassPoints.cxx
          This program uses a TTree of data events to produce a series of mass
	  points that can be used as inputs for the workspace. The cutflow is
	  implemented using the DMEvtSelect class.

     * DMSigParam.cxx
          This program uses MC to fit the resonance shape for the SM Higgs and 
	  the DM signal and saves the parameters for use in workspace 
	  generation. The cutflow is implemented using the DMEvtSelect class.

Supporting Classes:

     * BRXSReader.cxx
          Reads tables of SM Higgs cross sections and branching ratios and 
	  provides an easy-to-use interface.

     * DMEvtSelect.cxx
          This class implements the cutflow and counters for the analysis. It
	  can be initialized using a pointer to the DMTree. This class is 
	  instantiated in the DMMassPoints and DMSigParam classes. 

     * DMTree.cxx
          This class is automatically generated based on the MxAOD structure.
	  It provides a useful interface for code to the TTrees. 
