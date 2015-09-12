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

// C++ libraries:
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// ROOT includes:
#include "TString.h"

// Package includes:
#include "Config.h"

class BRXSReader 
{
  
 public:
  
  BRXSReader(TString configFile);
  virtual ~BRXSReader() {};

  // Accessors:
  float getSMBR(double mass, TString decay, TString value);
  float getDMXSBR(int massMediator, int massFermion, TString type,
		  TString value);
  float getMCXS(TString sampleName, TString value);
  float getSMXS(double mass, TString production, TString value);
  void printSMBR(double mass, TString decay);
  void printSMXS(double mass, TString production);
  
 private:
  
  // Private accessors & mutators:
  TString getDMMapKey(int massMediator, int massFermion, TString type,
		      TString value);
  TString getMCMapKey(TString sampleName, TString value);
  TString getSMMapKey(double mass, TString type, TString value);
  bool hasKey(TString key, TString mapType);
  void loadSMBR(TString decayClass);
  void loadSMXS(TString production);
  void loadDMXS();
  void loadMCXS();
  std::pair<double,double> getNearbySMMasses(double testMass, TString mapType);
  float getInterpolatedSMValue(double mass, TString mapType, TString process,
			       TString value);
  
  // Member objects:
  Config *m_config;
  TString directory;
  std::map<TString,float> valuesXS;
  std::map<TString,float> valuesBR;
  std::vector<double> massesHiggsXS;
  std::vector<double> massesHiggsBR;
  
};

#endif
