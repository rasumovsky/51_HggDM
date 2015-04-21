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
#include <cmath>

#include "TString.h"

class BRXSReader 
{
  
 public:
  
  BRXSReader(TString inputDirectory);
  virtual ~BRXSReader() {};

  // Accessors:
  float getSMBR(double mass, TString decay, TString value);
  float getDMXSBR(int massMediator, int massFermion, TString type,
		  TString value);
  float getSMXS(double mass, TString production, TString value);
  void printSMBR(double mass, TString decay);
  void printSMXS(double mass, TString production);
  
 private:
  
  // Private accessors & mutators:
  TString getDMMapKey(int massMediator, int massFermion, TString type,
		      TString value);
  TString getSMMapKey(double mass, TString type, TString value);
  bool hasKey(TString key, TString mapType);
  void loadSMBR(TString decayClass);
  void loadSMXS(TString production);
  void loadDMXS();
  std::pair<double,double> getNearbySMMasses(double testMass, TString mapType);
  float getInterpolatedSMValue(double mass, TString mapType, TString process,
			       TString value);
  
  // Member objects:
  TString directory;
  std::map<TString,float> valuesXS;
  std::map<TString,float> valuesBR;
  
  std::vector<double> *massesHiggsXS;
  std::vector<double> *massesHiggsBR;
  
};

#endif
