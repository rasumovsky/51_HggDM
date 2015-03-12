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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"

class DMMassPoints 
{
  
 public:
  
  DMMassPoints(std::string directory, 
	       std::string inFileName,
	       std::string outFileName);
  virtual ~DMMassPoints() {};
  
  // Member functions:
  int getNEvents(TString cateType, int cate);
  
  RooDataSet* createDataSet(TString cateType, int cateNum);
  TTree* createTTree(TString cateType, int cateNum);
  
 private:
  
  // count passing selection at each stage:
  int nPassing[20];
  double nPassing_weighted[20];
  TString cutNames[20];
  
  // count numbers in each category:
  int nEvents[20];
  double nEvents_weighted[20];
  
  TTree *tree;
  
  TString jobname;
  TString option;
    
};

#endif
