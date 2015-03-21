////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMBkgModel.h                                                        //
//  Class: DMBkgModel.cxx                                                     //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMBkgModel_h
#define DMBkgModel_h

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
#include "DMHeader.h"

class DMBkgModel 
{
  
 public:
  
  DMBkgModel(TString newJobName, TString newSampleName, TString newCateScheme,
	       TString newOptions);
  virtual ~DMBkgModel() {};
  
  // Accessors:
  TString getMassPointsFileName(int cateIndex);
  RooDataSet* getCateDataSet(int cateIndex);
  RooDataSet* getCombDataSet();
      
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
  
};

#endif
