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
#include <iostream>
#include <fstream>

#include "TString.h"

class BRXSReader 
{
  
 public:
  
  BRXSReader(TString inputDirectory);
  virtual ~BRXSReader() {};

  // Accessors:
  float getSMBR(double mass, TString decay, TString value);
  float getDMXSBR(int massIntermediate, int massFermion, TString type,
		  TString value);
  float getSMXS(double mass, TString production, TString value);
  void printSMBR(double mass, TString decay);
  void printSMXS(double mass, TString production);
  
 private:
  
  // Private accessors & mutators:
  TString getDMMapKey(int massIntermediate, int massFermion, TString type,
		      TString value);
  TString getSMMapKey(double mass, TString type, TString value);
  bool hasKey(TString key, TString mapType);
  void loadSMBR(TString decayClass);
  void loadSMXS(TString production);
  
  // Member objects:
  TString directory;
  std::map<TString,float> valuesXS;
  std::map<TString,float> valuesBR;
  
};

#endif
