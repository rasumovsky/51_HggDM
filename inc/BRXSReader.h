////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BRXSReader.h                                                        //
//  Class: BRXSReader.cxx                                                     //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BRXSReader_h
#define BRXSReader_h

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <map>

class BRXSReader 
{
  
 public:
  
  BRXSReader(TString inputDirectory);
  ~BRXSReader();
  
  // Accessors:
  float getBR(int mass, TString decay, TString value);
  float getXS(int mass, TString production, TString value);
  void printBR(int mass, TString decay, TString value);
  void printXS(int mass, TString production, TString value);
  
 private:
  
  // Private accessors & mutators:
  TString getMapKey(double mass, TString type, TString value);
  bool hasKey(TString key, TString mapType);
  void loadBR(TString decayClass);
  void loadXS(TString production);
  
  // Member objects:
  TString directory;
  std::map<TString,float> valuesXS;
  std::map<TString,float> valuesBR;
  
};

#endif
