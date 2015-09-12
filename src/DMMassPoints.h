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
  TString getMassPointsFileName(int cateIndex, TString sampleName);
  
  // Mutators:
  TString createLocalFilesAndList(TString originListName);
  void mergeMassPoints(TString newSampleName, DMMassPoints *inputMassPoints,
		       bool saveMassPoints);
  void setMassObservable(RooRealVar *newObservable);
  void removeLocalFilesAndList(TString listName);

 private:
  
  // Member methods:
  void createNewMassPoints();
  void fillHist1D(TString varName, bool allEvents, double xVal, double xWeight,
		  int cateIndex);
  void loadMassPointsFromFile();
  void newHist1D(TString varName, int nBins, double xMin, double xMax);
  void saveHists();
  
  // Member variables:
  TString m_sampleName;
  TString m_outputDir;
  bool m_isWeighted;
  TString m_options;
  
  Config *m_config;
  TString m_configFileName;
  
  RooDataSet *m_cateData[20];
  RooRealVar *m_yy;
  
  // Plotting information:
  std::map<TString,TH1F*> m_hists;
  
};

#endif
