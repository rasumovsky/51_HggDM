////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParamInterface.h                                                 //
//  Class: SigParamInterface.cxx                                              //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SigParamInterface_h
#define SigParamInterface_h

// C++ includes:
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

// ROOT includes:
#include "TString.h"

// Package includes:
#include "CommonFunc.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "SigParam.h"

class SigParamInterface 
{
  
 public:
  
  // Constructor / destructor:
  SigParamInterface(TString newJobName, TString newCateScheme,
		    TString newOptions);
  virtual ~SigParamInterface() {};
  
  // Accessors:
  bool allSignalsReady();
  bool createNew(TString signalType);
  RooDataSet *getData(TString signalType, int cateIndex);
  SigParam* getSigParam(TString signalType);
  bool loadFile(TString signalType);

 private:
  
  // Member variables:
  TString jobName;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  TString failedSigParam;
  bool signalsOK;
  
  std::map<TString,SigParam*> sigMap;
};

#endif
