////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: ESSReader.h                                                         //
//  Class: ESSReader.cxx                                                      //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ESSReader_h
#define ESSReader_h

#include <cmath>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"

class ESSReader {
 
 public:
  
  ESSReader(TString inputFileName, int nCategories);
  ~ESSReader();
  
  int sourceNameToIndex(TString name);
  double getValue(TString name, int cateIndex);
  int getSign(TString name, int cateIndex);
  int getNumberOfSources();
  TString getNameOfSource(int sourceIndex);
  
 private:
  
  int nESSParams; 
  
  // store the ESS values [#systematics][#categories]
  double valuesESS[100][20];
  
  // store the ESS names:
  std::vector<TString> nameListESS;
  
};

#endif
