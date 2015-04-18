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

  const double DMMyyRangeLo = 105.0;
  const double DMMyyRangeHi = 160.0;
  
  //int const nSMModes = 6;
  const int nSMModes = 6;
  const TString sigSMModes[nSMModes] = {"ggH","VBF","WH","ZH","bbH","ttH"};
  
  //int const nDMModes = 8;
  const int nDMModes = 8;
  const TString sigDMModes[nDMModes] = {"shxx_gg_ms100_mx100",
				  "shxx_gg_ms100_mx500",
				  "zphxx_gg_mzp100_mx100"};  
  
  //int const nMCProcesses = 1;
  const int nMCProcesses = 1;
  const TString MCProcesses[nMCProcesses] = {"gg_gjet"};
  
  ////////////////////////////////////////
  //    INPUT AND OUTPUT DIRECTORIES    //
  ////////////////////////////////////////
  
  // Location of global input files:
  const TString masterInput = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs";
  // Location of output directory:
  const TString masterOutput = "/afs/cern.ch/work/a/ahard/files_HDM/FullAnalysis";
  
  const TString fileNamePESValues = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/Systematics/PES/table_PES.txt";
  const TString fileNamePERValues = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/Systematics/PER/table_PER.txt";
  
  ////////////////////////////////////////
  //          SCRIPT LOCATIONS          //
  ////////////////////////////////////////
  
  //TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
  //TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";
  
  ////////////////////////////////////////
  //           MEMBER FUNCTIONS         //
  ////////////////////////////////////////
  
  TString nameToFileList(TString name);
  TString cateToBkgFunc(TString category);
  TString getIntermediaryName(TString modeName);
  int getIntermediaryMass(TString modeName);
  int getDarkMatterMass(TString modeName);
  bool isSMSample(TString sampleName);
  bool isDMSample(TString sampleName);
  bool isSignalSample(TString sampleName);
  bool isWeightedSample(TString sampleName);
};

#endif
