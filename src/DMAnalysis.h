////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMAnalysis.h                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 14/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _DMAnalysis_h_
#define _DMAnalysis_h_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include "TString.h"

namespace DMAnalysis {
  
  ////////////////////////////////////////
  //         GLOBAL PARAMETERS          //
  ////////////////////////////////////////
  
  // Set True for final analysis on data:
  const bool doBlind = true;
  
  // Luminosity in pb-1:
  const double analysisLuminosity = 5000;
  
  const double higgsMass = 125.09;// GeV

  const double DMMyyRangeLo = 105.0;// GeV
  const double DMMyyRangeHi = 160.0;// GeV

  const int nSMModes = 6;
  const TString sigSMModes[nSMModes] = {"ggH","VBF","WH","ZH","bbH","ttH"};
  
  const int nDMModes = 4;
  const TString sigDMModes[nDMModes] = {"shxx_gg_ms1000_mx500",
					"shxx_gg_ms100_mx100",
					"zphxx_gg_mzp1000_mx500",
					"zphxx_gg_mzp100_mx100"};
  
  const int nMCProcesses = 1;
  const TString MCProcesses[nMCProcesses] = {"gg_gjet"};
  
  ////////////////////////////////////////
  //    INPUT AND OUTPUT DIRECTORIES    //
  ////////////////////////////////////////
  
  // Location of global input files:
  const TString masterInput = "/afs/cern.ch/work/a/ahard/files_HggDM/GlobalInputs";
  // Location of output directory:
  const TString masterOutput = "/afs/cern.ch/work/a/ahard/files_HggDM/FullAnalysis";
  
  // Location of this software package:
  const TString packageLocation = "/afs/cern.ch/user/a/ahard/analysis/51_HggDM";
  
  // Holding location of cluster job files:
  const TString clusterFileLocation = "/afs/cern.ch/work/a/ahard/jobfiles";
  
  // Sub-directory for file lists:
  const TString fileListDir = "Apr21_15";
  
  // Locations of systematic uncertainty files:
  const TString fileNamePESValues = "/afs/cern.ch/work/a/ahard/files_HggDM/GlobalInputs/Systematics/PES/table_PES.txt";
  const TString fileNamePERValues = "/afs/cern.ch/work/a/ahard/files_HggDM/GlobalInputs/Systematics/PER/table_PER.txt";
  
  ////////////////////////////////////////
  //        FOR JOB SUBMISSION          //
  ////////////////////////////////////////
  
  const TString exeWorkspace = "DMWorkspaceWrapper";
  const TString jobScriptWorkspace = "scripts/jobFileWorkspace.sh";
  
  const TString exeTestStat = "DMTestStatWrapper";
  const TString jobScriptTestStat = "scripts/jobFileTestStat.sh";
  
  const TString exeMuLimit = "DMMuLimit";
  const TString jobScriptMuLimit = "scripts/jobFileMuLimit.sh";
  
  ////////////////////////////////////////
  //           MEMBER FUNCTIONS         //
  ////////////////////////////////////////
  
  int getNumCategories(TString cateScheme);
  TString nameToFileList(TString name);
  TString cateToBkgFunc(TString category);
  TString getMediatorName(TString modeName);
  int getMediatorMass(TString modeName);
  int getDarkMatterMass(TString modeName);
  bool isSMSample(TString sampleName);
  bool isDMSample(TString sampleName);
  bool isSignalSample(TString sampleName);
  bool isWeightedSample(TString sampleName);

};

#endif
