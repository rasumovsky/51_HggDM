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
#include "Config.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooBernsteinM.h"
#include "RooLandau.h"

class SystematicsTool {
  
 public:
  
  SystematicsTool(TString newConfigFile);
  virtual ~SystematicsTool() {};
  
  // Accessors:
  double calculateMigrSys(TString sysName, TString sampleName, int cateIndex);
  double calculateNormSys(TString sysName, TString sampleName);
  double getNormSys(TString sysName, TString sampleName);
  double getMigrSys(TString sysName, TString sampleName, int cateIndex);
  double getYield(TString sysName, TString sampleName);
  double getYield(TString sysName, TString sampleName, int cateIndex);
  std::vector<TString> listAllSys();
  
  // Mutators:
  void loadAllSys(TString sampleName); 
  void loadSingleSys(TString sysName, TString sampleName);
  std::vector<TString> rankMigrSysForSample(TString sampleName);
  std::vector<TString> rankNormSysForSample(TString sampleName);
  void saveRankedMigrSys(TString sampleName);
  void saveRankedNormSys(TString sampleName);
  void setNormSys(TString sysName, TString sampleName, double sysValue);
  void setMigrSys(TString sysName, TString sampleName, int cateIndex,
		  double sysValue);
  void setYield(TString sysName, TString sampleName, double yieldValue);
  void setYield(TString sysName, TString sampleName, int cateIndex,
		double yieldValue);
  
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
