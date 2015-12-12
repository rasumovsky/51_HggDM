////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SystematicsTool.h                                                   //
//  Class: SystematicsTool.cxx                                                //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/12/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SystematicsTool_h
#define SystematicsTool_h

// C++ includes:
/*
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
*/

// ROOT includes:
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"

// Package includes:
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooBernsteinM.h"
#include "RooLandau.h"

class SystematicsTool {
  
 public:
  
  SystematicsTool(TString newConfigFile);
  virtual ~SystematicsTool() {};
  
  // Accessors:
  void getNormSys(TString sysName, TString sampleName);
  void getMigrSys(TString sysName, TString sampleName, int cateIndex);
  std::vector<TString> listAllSys();
  
  // Mutators:
  void loadAllSys(TString sampleName); 
  void loadSingleSys(TString sysName, TString sampleName);
  std::vector<TString> rankSysForSample(TString sampleName);
  
 private:
  
  // Member methods:
  TString sysKey(TString sysName, TString sampleName, int cateIndex);
  TString sysKey(TString sysName, TString sampleName);
  
  // Member variables:
  TString m_configFile;
  TString m_outputDir;
  
  // Access to analysis settings:
  Config *m_config;
  
  // A map to store systematics values:
  std::map<TString,double> m_yieldStorage;
  std::map<TString,double> m_sysStorage;
  
};

#endif
