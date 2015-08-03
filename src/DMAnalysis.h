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
  
  // NOTE: All global variables have been moved to the .cfg files in data/
  // Member functions:
  int getNumCategories(TString cateScheme);
  TString nameToFileList(TString name);
  TString nameToxAODCutFile(TString name);
  TString cateToBkgFunc(TString cateScheme, int cateIndex);
  TString getMediatorName(TString modeName);
  int getMediatorMass(TString modeName);
  int getDarkMatterMass(TString modeName);
  bool isSMSample(TString sampleName);
  bool isDMSample(TString sampleName);
  bool isSignalSample(TString sampleName);
  bool isWeightedSample(TString sampleName);

};

#endif
