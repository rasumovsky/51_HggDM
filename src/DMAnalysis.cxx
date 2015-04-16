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
    @param name - the name of the sample.
    @returns - the file list location.
*/
TString DMAnalysis::nameToFileList(TString name) {
  TString result = Form("%s/FileLists/",masterInput.Data());
  if (name.EqualTo("ggH")) {
    result += "list_H2yyMETAnalysis_ggH.txt";
  }
  else if (name.EqualTo("VBF")) {
    result += "list_H2yyMETAnalysis_VBF.txt";
  }
  else if (name.EqualTo("WH")) {
    result += "list_H2yyMETAnalysis_WH.txt";
  }
  else if (name.EqualTo("ZH")) {
    result += "list_H2yyMETAnalysis_ZH.txt";
  }
  else if (name.EqualTo("ttH")) {
    result += "list_H2yyMETAnalysis_ttH.txt";
  }
  else if (name.EqualTo("gg_gjet")) {
    result += "list_H2yyMETAnalysis_gg_gjet.txt";
  }
  else if (name.EqualTo("shxx_gg_ms100_mx100")) {
    result += "list_H2yyMETAnalysis_shxx_gg_ms100_mx100.txt";
  }
  else if (name.EqualTo("shxx_gg_ms100_mx500")) {
    result += "list_H2yyMETAnalysis_shxx_gg_ms100_mx500.txt";
  }
  else if (name.EqualTo("zphxx_gg_mzp100_mx100")) {
    result += "list_H2yyMETAnalysis_zphxx_gg_mzp100_mx100.txt";
  }
  else {
    std::cout << "nameToRootFile: Error! No corresponding file" << std::endl;
  }
  return result;
}

/** 
    Determine the background PDF based on the category.
    @param category - the category name. 
    @returns - the name of the background PDF.
*/
TString DMAnalysis::cateToBkgFunc(TString category) {
  TString result = "";
  result = "Exppol01";
  //Possibilities are "BernO1",... "BernO6", "ExppolO1",... "ExppolO6"
  return result;
}

/**
   Convert the DM sample name into a process name, a subset of the full name.
   @param modeName - the name of the DM production mode.
   @param processName - the name of the process.
*/
TString DMAnalysis::getIntermediaryName(TString modeName) {
  if (modeName.Contains("shxx_gg")) {
    return "shxx_gg";
  }
  else if (modeName.Contains("zphxx_gg")) {
    return "zphxx_gg";
  }
  else {
    std::cout << "Analysis Error: no matching intermediary name" << std::endl;
    return "";
  }
}

/** 
    Convert the DM sample name into an intermediary mass.
    @param modeName - the name of the DM production mode.
    @returns - the mass of the mediator particle.
*/
int DMAnalysis::getIntermediaryMass(TString modeName) {
  for (int currMass = 100; currMass < 1000; currMass += 100) {
    if (modeName.Contains(Form("ms%d",currMass))) {
      return currMass;
    }
  }
  std::cout << "Analysis Error: no matching intermediary mass" << std::endl;
  return 0;
}

/** 
    Convert the DM sample name into a dark matter particle mass.
    @param modeName - the name of the DM production mode.
    @returns - the mass of the DM particle.
*/
int DMAnalysis::getDarkMatterMass(TString modeName) {
  for (int currMass = 100; currMass < 1000; currMass += 100) {
    if (modeName.Contains(Form("mx%d",currMass))) {
      return currMass;
    }
  }
  std::cout << "Analysis Error: no matching DM mass" << std::endl;
  return 0;
}

/** 
    Check if the sample is among those listed as a SM or DM signal. 
    @param sampleName - the name of the sample being used.
    @returns - true iff the sample is a signal sample.
*/
bool DMAnalysis::isSMSample(TString sampleName) {
  for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
    if (sampleName.EqualTo(sigSMModes[i_SM])) {
      return true;
    }
  }
  return false;
}

bool DMAnalysis::isDMSample(TString sampleName) {
  for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
    if (sampleName.EqualTo(sigDMModes[i_DM])) {
      return true;
    }
  }
  return false;
}

bool DMAnalysis::isSignalSample(TString sampleName) {
  if (isSMSample(sampleName) || isDMSample(sampleName)) {
    return true;
  }
  else {
    return false;
  }
}

/** 
    Check whether a sample should be weighted.
    @param sampleName - the name of the sample being used.
    @returns - true iff the sample has associated event weights.
*/
bool DMAnalysis::isWeightedSample(TString sampleName) {
  // First check if it is a SM or DM signal process:
  if (isSignalSample(sampleName)) {
    return true;
  }
  // Finally, check if it is one of the other MC processes:
  for (int i_MC = 0; i_MC < nMCProcesses; i_MC++) {
    if (sampleName.EqualTo(MCProcesses[i_MC])) {
      return true;
    }
  }
  return false;
}



#endif