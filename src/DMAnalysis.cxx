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
   -----------------------------------------------------------------------------
   Convert the DM sample name into a process name, a subset of the full name.
   @param modeName - the name of the DM production mode.
   @param processName - the name of the process.
*/
TString DMAnalysis::getMediatorName(TString modeName) {
  if (modeName.Contains("shxx_gg")) return "shxx_gg";
  else if (modeName.Contains("zphxx_gg")) return "zphxx_gg";
  else if (modeName.Contains("zp2hdmxx_gg")) return "zp2hdmxx_gg";
  else {
    std::cout << "Analysis Error: no matching mediator name: " 
	      << modeName << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Convert the DM sample name into a process name, a subset of the full name.
   @param modeName - the name of the DM production mode.
   @param processName - the name of the process.
*/
TString DMAnalysis::getPrintMediatorName(TString modeName) {
  if (modeName.Contains("shxx_gg")) return "Scalar";
  else if (modeName.Contains("zphxx_gg")) return "Z'";
  else if (modeName.Contains("zp2hdmxx_gg")) return "Z_{2HDM}'";
  else {
    std::cout << "Analysis Error: no matching mediator name: " 
	      << modeName << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a pretty LaTex formatted name for a variable.
   @param sampleName - The name of the sample.
   @returns - A LaTex formatted TString.
*/
TString DMAnalysis::getPrintSampleName(Config *config, TString sampleName) {
  if (sampleName.EqualTo("gg")) return "#it{#gamma#gamma}";
  else if (sampleName.EqualTo("gjet")) return "#it{#gamma}+jet";
  else if (sampleName.EqualTo("SMHiggs")) {
    return "SM #it{H#rightarrow#gamma#gamma}";
  }
  else if (sampleName.EqualTo("yy")) return "#it{#gamma#gamma}";
  else if (sampleName.EqualTo("yjet")) return "#it{#gamma}+jet";
  else if (sampleName.EqualTo("jj")) return "di-jet";
  else if (sampleName.EqualTo("Zyy")) return "#it{Z#gamma#gamma}";
  else if (sampleName.EqualTo("Wyy")) return "#it{W#gamma#gamma}";
  else if (sampleName.EqualTo("Zy")) return "#it{Z#gamma}";
  else if (sampleName.EqualTo("Wy")) return "#it{W#gamma}";
  else if (sampleName.EqualTo("ttbar")) return "#it{t#bar{t}}";
  else if (sampleName.EqualTo("diboson")) return "#it{WW}+#it{WZ}+#it{ZZ}";
  else if (sampleName.EqualTo("ggH")) {
    return "#it{ggH(#rightarrow#gamma#gamma)}";
  }
  else if (sampleName.EqualTo("VBF")) {
    return "#it{VBF (H#rightarrow#gamma#gamma)}";
  }
  else if (sampleName.EqualTo("WH")) {
    return "#it{WH(#rightarrow#gamma#gamma)}";
  }
  else if (sampleName.EqualTo("ZH")) {
    return "#it{ZH(#rightarrow#gamma#gamma)}";
  }
  else if (sampleName.EqualTo("bbH")) {
    return "#it{b#bar{b}H(#rightarrow#gamma#gamma)}";
  }
  else if (sampleName.EqualTo("ttH")) {
    return "#it{t#bar{t}H(#rightarrow#gamma#gamma)}";
  }
  else if (DMAnalysis::isDMSample(config, sampleName)) {
    TString mediatorName = DMAnalysis::getMediatorName(sampleName);
    int mediatorMass = DMAnalysis::getMediatorMass(config, sampleName);
    TString medMassForm;
    if (mediatorName.EqualTo("shxx_gg")) {
      medMassForm = Form("#it{m_{S}}=%dGeV", mediatorMass);
    }
    else if (mediatorName.EqualTo("zphxx_gg")) {
      medMassForm = Form("#it{m_{Z'}}=%dGeV", mediatorMass);
    }
    else if (mediatorName.EqualTo("zp2hdmxx_gg")) {
      medMassForm = Form("#it{m_{Z'}}=%dGeV", mediatorMass);
    }
    
    int darkMatterMass = DMAnalysis::getDarkMatterMass(config, sampleName);
    TString darkMassForm;
    if (mediatorName.EqualTo("shxx_gg")) {
      darkMassForm = Form("#it{m_{#chi}}=%dGeV", darkMatterMass);
    }
    else if (mediatorName.EqualTo("zphxx_gg")) {
      darkMassForm = Form("#it{m_{#chi}}=%dGeV", darkMatterMass);
    }
    else if (mediatorName.EqualTo("zp2hdmxx_gg")) {
      darkMassForm = Form("#it{m_{A}}=%dGeV", darkMatterMass);
    }
    //(DMAnalysis::getPrintMediatorName(sampleName)).Data(),
    TString newName = Form("%s %s", medMassForm.Data(), darkMassForm.Data());
    return newName;
  }
  else return sampleName;
}

/**
   -----------------------------------------------------------------------------
   Get a pretty LaTex formatted name for a variable.
   @param varName - The name of the variable in the file.
   @returns - A LaTex formatted TString.
*/
TString DMAnalysis::getPrintVarName(TString varName) {
  varName.ToLower();
  varName = varName.ReplaceAll("anacut", "");
  if (varName.EqualTo("ptyy") || varName.EqualTo("diphotonpt")) {
    return "p_{T}^{#gamma#gamma} [GeV]";
  }
  else if (varName.EqualTo("etmiss")) return "#slash{E}_{T} [GeV]";
  else if (varName.EqualTo("ratioetmissptyy")) {
    return "#slash{E}_{T}/p_{T}^{#gamma#gamma}";
  }
  else if (varName.EqualTo("myy")) return "M_{#gamma#gamma} [GeV]";
  else if (varName.Contains("atanratio")) {
    return "tan^{-1}(#slash{E}_{T}/p_{T}^{#gamma#gamma})";
  }
  else if (varName.EqualTo("sumsqrtetmissptyy")) {
    return "#sqrt{(p_{T}^{#gamma#gamma})^{2} + (#slash{E}_{T})^{2}} [GeV]";
  }
  else if (varName.Contains("dphiyyetmiss")) {
    return "#Delta#phi(#gamma#gamma-#slash{E}_{T})";
  }
  else if (varName.Contains("njets")) return "N_{Jets}";
  else if (varName.Contains("nleptons")) return "N_{Leptons}";
  else if (varName.Contains("p0")) {
    if (varName.Contains("exp")) return "p_{0} exp.";
    else if (varName.Contains("obs")) return "p_{0} obs.";
    else return "p_{0}";
  }
  else if (varName.Contains("cl")) {
    if (varName.Contains("exp")) return "#it{CL} exp.";
    else if (varName.Contains("obs")) return "#it{CL} obs.";
    else return "#it{CL}";
  }
  else return varName;
}

/**
   -----------------------------------------------------------------------------
   Convert the DM sample name into a mediator mass.
   @param modeName - the name of the DM production mode.
   @returns - the mass of the mediator particle.
*/
int DMAnalysis::getMediatorMass(Config *config, TString modeName) {
  std::vector<double> massList = config->getNumV("DMMediatorMasses");
  for (int i_m = (int)massList.size()-1; i_m >= 0; i_m--) {
    if (modeName.Contains(Form("ms%d",((int)massList[i_m]))) || 
	modeName.Contains(Form("mzp%d",((int)massList[i_m])))) {
      return ((int)massList[i_m]);
    }
  }
  std::cout << "Analysis Error: no matching mediator mass "
	    << modeName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Convert the DM sample name into a dark matter particle mass.
   @param modeName - the name of the DM production mode.
   @returns - the mass of the DM particle.
*/
int DMAnalysis::getDarkMatterMass(Config *config, TString modeName) {
  std::vector<double> massList = config->getNumV("DMParticleMasses");
  for (int i_m = (int)massList.size()-1; i_m >= 0; i_m--) {
    if (modeName.Contains(Form("mx%d",((int)massList[i_m]))) || 
	modeName.Contains(Form("mA%d",((int)massList[i_m])))) {
      return ((int)massList[i_m]);
    }
  }
  std::cout << "Analysis Error: no matching DM mass" 
	    << modeName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a background sample.
   @param config - The config file for the analysis settings.
   @param sampleName - The name of the sample being used.
   @returns - True iff the sample is a background sample.
*/
bool DMAnalysis::isBkgSample(Config *config, TString sampleName) {
  std::vector<TString> bkgProcesses = config->getStrV("BkgProcesses");
  for (int i_bkg = 0; i_bkg < (int)bkgProcesses.size(); i_bkg++) {
    if (sampleName.EqualTo(bkgProcesses[i_bkg])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a DM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - The name of the sample being used.
   @returns - True iff the sample is a DM signal sample.
*/
bool DMAnalysis::isDMSample(Config *config, TString sampleName) {
  std::vector<TString> sigDMModes = config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    if (sampleName.EqualTo(sigDMModes[i_DM])) return true;
  }
  return false;
}

/** 
   -----------------------------------------------------------------------------
   Check if a sample is skimmed.
   @param config - The config file for the analysis settings.
   @param fileName - The name of the file.
   @returns - True iff the sample is skimmed.
*/
bool DMAnalysis::isSkimmed(Config *config, TString fileName) {
  std::vector<TString> skimmedTags = config->getStrV("TagSkimmed");
  std::vector<TString> unskimmedTags = config->getStrV("TagUnskimmed");
  
  for (int i_s = 0; i_s < (int)skimmedTags.size(); i_s++) {
    if (fileName.Contains(skimmedTags[i_s])) {
      std::cout << "DMAnalysis: Sample is skimmed." << std::endl;
      return true;
    }
  }
  
  for (int i_u = 0; i_u < (int)unskimmedTags.size(); i_u++) {
    if (fileName.Contains(unskimmedTags[i_u])) {
      std::cout << "DMAnalysis: Sample is unskimmed." << std::endl;
      return false;
    }
  }
  
  // Exit with error message if tag has not been identified:
  std::cout << "DMAnalysis: ERROR! Unidentified tag: " << fileName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a SM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - The name of the sample being used.
   @returns - True iff the sample is a SM signal sample.
*/
bool DMAnalysis::isSMSample(Config *config, TString sampleName) {
  std::vector<TString> sigSMModes = config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    if (sampleName.EqualTo(sigSMModes[i_SM])) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if the sample is among those listed as a DM or SM signal. 
   @param config - The config file for the analysis settings.
   @param sampleName - The name of the sample being used.
   @returns - True iff the sample is a DM or SM signal sample.
*/
bool DMAnalysis::isSignalSample(Config *config, TString sampleName) {
  return (isSMSample(config, sampleName) || isDMSample(config, sampleName));
}

/**
   -----------------------------------------------------------------------------
   Check whether a sample should be weighted.
   @param config - The config file for the analysis settings.
   @param sampleName - The name of the sample being used.
   @returns - True iff the sample has associated event weights.
*/
bool DMAnalysis::isWeightedSample(Config *config, TString sampleName) {
  // First check if it is a SM or DM signal process:
  if (isDMSample(config, sampleName) || isSMSample(config, sampleName)) {
    return true;
  }
  
  // Finally, check if it is one of the other MC processes:
  std::vector<TString> BkgProcesses = config->getStrV("BkgProcesses");
  for (int i_b = 0; i_b < (int)BkgProcesses.size(); i_b++) {
    if (sampleName.EqualTo(BkgProcesses[i_b])) return true;
  }
  return false;
}

/** 
   -----------------------------------------------------------------------------
    Convert the sample name to the corresponding file list.
    @param config - The config file for the analysis settings.
    @param name - The name of the sample.
    @returns - The file list location.
*/
TString DMAnalysis::nameToFileList(Config *config, TString name) {
  // Create a new file list:
  TString newListName = Form("%s/FileLists/list_%s.txt", 
			     (config->getStr("masterInput")).Data(),
			     name.Data());
  ofstream newFileList(newListName);
  
  // Use a format specifier for DM samples:
  if (isDMSample(config, name)) {
    TString MxAODName = config->getStr("MxAODForm_DM");
    MxAODName = MxAODName.ReplaceAll("NAME",name);
    newFileList << config->getStr("MxAODDirectoryMC") << "/" << MxAODName 
		<< std::endl;
  }
  // Use MxAOD list in the settings file for SM Higgs, backgrounds, and data:
  else {
    std::vector<TString> listMxAODs
      = config->getStrV(Form("MxAODList_%s", name.Data()));
    // loop over lines in the file, each of which is one contributing MxAOD:
    for (int i_m = 0; i_m < (int)listMxAODs.size(); i_m++) {
      TString directory = listMxAODs[i_m].Contains("data") ? 
	config->getStr("MxAODDirectoryData") : 
	config->getStr("MxAODDirectoryMC");
      newFileList << directory << "/" << listMxAODs[i_m] << std::endl;
    }
  }
  newFileList.close();
  return newListName;
}
