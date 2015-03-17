# HggDM
Higgs to diphoton + dark matter search 

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
