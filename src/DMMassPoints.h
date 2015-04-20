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

// C++ includes:
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// ROOT includes:
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"

// Package includes:
#include "CommonFunc.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "DMEvtSelect.h"
#include "DMAnalysis.h"
#include "DMTree.h"
#include "BRXSReader.h"

using namespace DMAnalysis;

class DMMassPoints {
  
 public:
  
  DMMassPoints(TString newJobName, TString newSampleName, TString newCateScheme,
	       TString newOptions);
  DMMassPoints(TString newJobName, TString newSampleName, TString newCateScheme,
	       TString newOptions, RooRealVar *newObservable);
  DMMassPoints(TString newJobName, TString newSampleName, TString newCateScheme,
	       TString newOptions, RooRealVar *newObservable,
	       RooCategory *newCategories);
  virtual ~DMMassPoints() {};
  
  // Accessors:
  RooDataSet* getCateDataSet(int cateIndex);
  RooDataSet* getCombDataSet();
  RooRealVar* getMassObservable();
  TString getMassPointsFileName(int cateIndex);
  RooCategory* getRooCategory();

  // Mutators:
  void setMassObservable(RooRealVar *newObservable);
  void setRooCategory(RooCategory *newCategories);
  
 private:
  
  // Member methods:
  void createNewMassPoints();
  void loadMassPointsFromFile();
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  
  TString outputDir;
  bool isWeighted;
  
  RooDataSet *cateData[20];
  RooDataSet *combData;
  RooRealVar *m_yy;
  RooCategory *categories;
};

#endif
