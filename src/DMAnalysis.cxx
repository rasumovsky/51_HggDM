////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMAnalysis.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
//  This namespace stores all of the global information for the H->gg + DM    //
//  search with 13 TeV data in 2015. It also has all of the includes that are //
//  necessary for the analysis.                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMAnalysis.h"

/** 
    Convert the sample name to the corresponding file list.
    @param config - The config file for the analysis settings.
    @param name - the name of the sample.
    @returns - the file list location.
*/
TString DMAnalysis::nameToFileList(Config *config, TString name) {
  TString result = Form("%s/FileLists/%s/list_H2yyMETAnalysis_%s.txt",
			(config->getStr("masterInput")).Data(),
			(config->getStr("fileListDir")).Data(), name.Data());
  // NOTE: use ggH for bbH, weight sigma_bbH/sigma_ggH:
  if (name.EqualTo("bbH")) result = nameToFileList(config, "ggH");
  return result;
}

/**
   Converts the sample name to the corresponding cutflow histogram.
   @param config - The config file for the analysis settings.
   @param name - the name of the sample.
   @returns - the file list location.
*/
TString DMAnalysis::nameToxAODCutFile(Config *config, TString name) {
  TString result = Form("%s/xAODCutFlows/%s/hist_H2yyMETAnalysis_%s.root",
			(config->getStr("masterInput")).Data(), 
			(config->getStr("fileListDir")).Data(), name.Data());
  // NOTE: use ggH for bbH, weight sigma_bbH/sigma_ggH:
  if (name.EqualTo("bbH")) result = nameToxAODCutFile(config, "ggH");
  return result;
}

/**
   Convert the DM sample name into a process name, a subset of the full name.
   @param modeName - the name of the DM production mode.
   @param processName - the name of the process.
*/
TString DMAnalysis::getMediatorName(TString modeName) {
  if (modeName.Contains("shxx_gg")) return "shxx_gg";
  else if (modeName.Contains("zphxx_gg")) return "zphxx_gg";
  else {
    std::cout << "Analysis Error: no matching mediator name: " 
	      << modeName << std::endl;
    exit(0);
  }
}

/**
   Convert the DM sample name into a mediator mass.
   @param modeName - the name of the DM production mode.
   @returns - the mass of the mediator particle.
*/
int DMAnalysis::getMediatorMass(TString modeName) {
  for (int currMass = 100; currMass < 1000; currMass += 100) {
    if (modeName.Contains(Form("ms%d",currMass)) || 
	modeName.Contains(Form("mzp%d",currMass))) {
      return currMass;
    }
  }
  std::cout << "Analysis Error: no matching mediator mass"
	    << modeName << std::endl;
  exit(0);
}

/**
   Convert the DM sample name into a dark matter particle mass.
   @param modeName - the name of the DM production mode.
   @returns - the mass of the DM particle.
*/
int DMAnalysis::getDarkMatterMass(TString modeName) {
  for (int currMass = 100; currMass < 1000; currMass += 100) {
    if (modeName.Contains(Form("mx%d",currMass))) return currMass;
  }
  std::cout << "Analysis Error: no matching DM mass" 
	    << modeName << std::endl;
  exit(0);
}

/**
   Check if the sample is among those listed as a SM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a SM signal sample.
*/
bool DMAnalysis::isSMSample(Config *config, TString sampleName) {
  std::vector<TString> sigSMModes = config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    if (sampleName.EqualTo(sigSMModes[i_SM])) return true;
  }
  return false;
}

/**
   Check if the sample is among those listed as a DM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a DM signal sample.
*/
bool DMAnalysis::isDMSample(Config *config, TString sampleName) {
  std::vector<TString> sigDMModes = config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    if (sampleName.EqualTo(sigDMModes[i_DM])) return true;
  }
  return false;
}

/**
   Check if the sample is among those listed as a DM or SM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample is a DM or SM signal sample.
*/
bool DMAnalysis::isSignalSample(Config *config, TString sampleName) {
  return (isSMSample(config, sampleName) || isDMSample(config, sampleName));
}

/**
   Check whether a sample should be weighted.
   @param config - The config file for the analysis settings.
   @param sampleName - the name of the sample being used.
   @returns - true iff the sample has associated event weights.
*/
bool DMAnalysis::isWeightedSample(Config *config, TString sampleName) {
  // First check if it is a SM or DM signal process:
  if (isSignalSample(config, sampleName)) return true;
  
  // Finally, check if it is one of the other MC processes:
  std::vector<TString> MCProcesses = config->getStrV("MCProcesses");
  for (int i_MC = 0; i_MC < (int)MCProcesses.size(); i_MC++) {
    if (sampleName.EqualTo(MCProcesses[i_MC])) return true;
  }
  return false;
}
