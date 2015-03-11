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

class DMMassPoints 
{
  
 public:
  
  DMMassPoints(std::string directory, 
	       std::string inFileName,
	       std::string outFileName);
  virtual ~DMMassPoints() {};
  
  // Member functions:
  int getTotalEvents();
  std::vector<int> getEvtPerCategory(std::string cateType);
  RooDataSet* createDataSet(std::string cateType, int cateNum);
  TTree* createTTree(std::string cateType, int cateNum);
  
 private:
  
  int nEvents;
  TTree *tree;
    
};

#endif
