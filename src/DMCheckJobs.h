////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMCheckJobs.h                                                       //
//  Class: DMCheckJobs.cxx                                                    //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMCheckJobs_h
#define DMCheckJobs_h

// C++ libraries:
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>

// ROOT libraries:
#include "TString.h"

// Package libraries:
#include "Config.h"
#include "DMAnalysis.h"

class DMCheckJobs {

 public:
  
  // Constructor and destructor:
  DMCheckJobs(TString configFileName);
  virtual ~DMCheckJobs() {};
  
  // Mutators:
  int getNumberToResubmit(TString jobType);
  std::vector<TString> getResubmitList(TString jobType);
  void updateJobStatus(TString jobType);
  
  // Accessors:
  void printResubmitList(TString jobType);
  
 private:
  
  Config *m_config;

  std::vector<TString> m_listDMWorkspace;
  std::vector<TString> m_listDMTestStat;
  std::vector<TString> m_listDMMuLimit;
  
};

#endif
