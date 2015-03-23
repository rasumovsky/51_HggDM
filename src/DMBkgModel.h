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
#include "DMEvtSelect.h"
#include "RooBernsteinM.h"

class DMBkgModel 
{
  
 public:
  
  DMBkgModel(TString newJobName, TString newSampleName, TString newCateScheme,
	     TString newOptions);
  DMBkgModel(TString newJobName, TString newSampleName, TString newCateScheme,
	     TString newOptions, RooRealVar *newObservable);
  DMBkgModel(TString newJobName, TString newSampleName, TString newCateScheme,
	     TString newOptions, RooRealVar *newObservable,
	     RooCategory *newCategories);
  
  virtual ~DMBkgModel() {};
  
  // Accessors:
  //void fitCateBkgPDF(int cateIndex);
  //void fitCombBkgPDF();
  RooAbsPdf* getCateBkgPDF(int cateIndex);
  RooSimultaneous* getCombBkgPDF();
  RooAbsPdf* getBkgPDFByName(TString name, TString fitFunc);
  RooRealVar* getMassObservable();
  RooCategory* getRooCategory();
  
  // Mutators:
  void setRooCategory(RooCategory *newCategories);
  void setMassObservable(RooRealVar *newObservable);
      
 private:
  
  // Member methods:
  int getOrderFromFunc(TString fitFunc);
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  
  int nCategories;
  RooRealVar *m_yy;
  RooCategory *categories;
  DMEvtSelect *selector;
  
};

#endif
