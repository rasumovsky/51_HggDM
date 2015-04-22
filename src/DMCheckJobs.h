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

#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include "TString.h"

#include "DMAnalysis.h"

class DMCheckJobs {

 public:
  
  // Constructor and destructor:
  DMCheckJobs(TString newJobName);
  virtual ~DMCheckJobs() {};
  
  // Mutators:
  int getNumberToResubmit(TString jobType);
  std::vector<TString> getResubmitList(TString jobType);
  void updateJobStatus(TString jobType);
  
  // Accessors:
  void printResubmitList(TString jobType);
  
 private:
    
  TString jobName;
  std::vector<TString> listDMWorkspace;
  std::vector<TString> listDMTestStat;
  std::vector<TString> listDMMuLimit;
  
};

#endif
