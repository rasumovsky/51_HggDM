////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMassPoints.h                                                      //
//  Class: DMMassPoints.cxx                                                   //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMMassPoints_h
#define DMMassPoints_h

// C++ libraries:
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// ROOT libraries:
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"

// Package libraries:
#include "Config.h"
#include "CommonFunc.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "BRXSReader.h"
#include "DMAnalysis.h"
#include "DMEvtSelect.h"
#include "DMTree.h"
#include "DMxAODCutflow.h"

class DMMassPoints {
  
 public:
  
  DMMassPoints(TString newConfigFile, TString newSampleName, TString newOptions,
	       RooRealVar *newObservable);
  virtual ~DMMassPoints() {};
  
  // Accessors:
  RooDataSet* getCateDataSet(int cateIndex);
  RooRealVar* getMassObservable();
  TString getMassPointsFileName(int cateIndex);
  
  // Mutators:
  void setMassObservable(RooRealVar *newObservable);
   
 private:
  
  // Member methods:
  void createNewMassPoints();
  void loadMassPointsFromFile();
  
  // Member variables:
  TString m_sampleName;
  TString m_outputDir;
  bool m_isWeighted;
  
  Config *m_config;
  
  RooDataSet *m_cateData[20];
  RooRealVar *m_yy;

};

#endif
