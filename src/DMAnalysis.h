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
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "Config.h"

namespace DMAnalysis {
  
  // Member functions:
  TString nameToFileList(Config *config, TString name);
  TString nameToxAODCutFile(Config *config, TString name);
  TString getMediatorName(TString modeName);
  int getMediatorMass(TString modeName);
  int getDarkMatterMass(TString modeName);
  bool isSMSample(Config *config, TString sampleName);
  bool isDMSample(Config *config, TString sampleName);
  bool isSignalSample(Config *config, TString sampleName);
  bool isWeightedSample(Config *config, TString sampleName);

};

#endif
