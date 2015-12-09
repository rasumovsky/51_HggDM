////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMAnalysis.h                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 03/08/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef _DMAnalysis_h_
#define _DMAnalysis_h_

// C++ libraries:
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "Config.h"

namespace DMAnalysis {
  
  // Member functions:
  TString getMediatorName(TString modeName);
  TString getPrintMediatorName(TString modeName);
  TString getPrintSampleName(Config *config, TString sampleName);
  TString getPrintVarName(TString varName);
  int getMediatorMass(Config *config, TString modeName);
  int getDarkMatterMass(Config *config, TString modeName);
  bool isBkgSample(Config *config, TString sampleName);
  bool isDMSample(Config *config, TString sampleName);
  bool isSkimmed(Config *config, TString fileName);
  bool isSMSample(Config *config, TString sampleName);
  bool isSignalSample(Config *config, TString sampleName);
  bool isWeightedSample(Config *config, TString sampleName);
  TString nameToFileList(Config *config, TString name, bool useSys);  
  
};

#endif
