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
#include "statistics.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"

using namespace DMAnalysis;

class DMSigParam 
{
  
 public:
  
  DMSigParam(TString newJobName, TString newCateScheme, TString newOptions, 
	     RooRealVar *newObservable);
  virtual ~DMSigParam() {};
  
  // Accessors:
  void addSigToCateWS(RooWorkspace *&workspace, std::vector<TString> namesESS,
		      std::vector<TString> namesRes, TString process,
		      int cateIndex);
  RooCBShape* getCateCrystalBall(int cateIndex, TString process);
  RooGaussian* getCateGaussian(int cateIndex, TString process);
  RooAddPdf* getCateSigPDF(int cateIndex, TString process);
  double getCateSigYield(int cateIndex, TString process);
  double getCombSigYield(TString process);
  RooRealVar* getMassObservable(); 
  double getSigParam(TString process, TString param, int cateIndex);
  TString getSigParamFileName(TString process, TString fileType, int cateIndex);

  // Mutators:
  void setMassObservable(RooRealVar *newObservable);
  
 private:
  
  // Member methods:
  void createSigParam(TString process, bool makeNew);
  
  // Member variables:
  TString jobName;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  RooRealVar *m_yy;
  
  std::map<TString,std::vector<RooCBShape*> > sigCB;
  std::map<TString,std::vector<RooGaussian*> > sigGA;
  std::map<TString,std::vector<RooAddPdf*> > sigPDF;
  std::map<TString,std::vector<double> > sigYield;
  
};

#endif
