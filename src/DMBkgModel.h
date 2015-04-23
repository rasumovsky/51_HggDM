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
#include "DMAnalysis.h"
#include "RooBernsteinM.h"

using namespace DMAnalysis;

class DMBkgModel 
{
  
 public:
  
  DMBkgModel(TString newJobName, TString newCateScheme, TString newOptions);
  DMBkgModel(TString newJobName, TString newCateScheme, TString newOptions,
	     RooRealVar *newObservable);
  DMBkgModel(TString newJobName, TString newCateScheme, TString newOptions,
	     RooRealVar *newObservable, RooCategory *newCategories);
  
  virtual ~DMBkgModel() {};
  
  // Accessors:
  //void fitCateBkgPDF(int cateIndex);
  //void fitCombBkgPDF();
  void addBkgToCateWS(RooWorkspace *&workspace, RooArgSet *&nuisParams,
		      int cateIndex);
  RooAbsPdf* getCateBkgPDF(int cateIndex);
  RooAbsPdf* getBkgPDFByName(TString name, TString fitFunc);
  RooRealVar* getMassObservable();
  
  // Mutators:
  void setMassObservable(RooRealVar *newObservable);
      
 private:
  
  // Member methods:
  int getOrderFromFunc(TString fitFunc);
  
  // Member variables:
  TString jobName;
  TString cateScheme;
  TString options;
  
  RooRealVar *m_yy;
  
};

#endif
