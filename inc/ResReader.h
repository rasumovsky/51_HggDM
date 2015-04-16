////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: ResReader.h                                                         //
//  Class: ResReader.cxx                                                      //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ResReader_h
#define ResReader_h

#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

class ResReader {

 public:
  
  ResReader(TString inputFileName, int nCategories);
  ~ResReader();
  
  int sourceNameToIndex(TString name); 
  double getValue(TString name, int cateIndex);
  int getSign(TString name, int cateIndex);
  int getNumberOfSources();
  TString getNameOfSource(int resIndex);
  
 private:
  
  int nResParams;
  
  // store the RES values [#systematics][#categories]
  double valuesRes[100][20];
  
  // store the RES names:
  vector<TString> nameListRes;
  
};

#endif
