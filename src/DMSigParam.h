////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMSigParam.h                                                        //
//  Class: DMSigParam.cxx                                                     //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMSigParam_h
#define DMSigParam_h

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
#include "CommonHead.h"
#include "RooFitHead.h"
#include "DMEvtSelect.h"
#include "DMHeader.h"
#include "DMTree.h"

class DMSigParam 
{
  
 public:
  
  DMSigParam(TString newJobName, TString newSampleName, TString newCateScheme,
	     TString newOptions);
  virtual ~DMSigParam() {};
  
  // Accessors:
  RooAddPdf* getCateSigPDF(int cateIndex, TString process);
  double getCateSigYield(int cateIndex, TString process);
  //RooSimultaneous* getCombSigPDF(TString process); // not sure of class...
  double getCombSigYield(TString process);
  TString getSigParamFileName(int cateIndex, TString production,
			      TString fileType);
  
 private:
  
  // Member methods:
  void createNewSigParam();
  void loadSigParamFromFile();//simple method that just counts from text file.
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  RooAddPdf* sigPDF[20];
  
};

#endif
