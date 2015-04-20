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

#include "DMAnalysis.h"

class DMCheckJobs {

 public:
  
  // Constructor and destructor:
  DMCheckJobs(TString newJobName);
  virtual ~DMBkgModel() {};
  
  // Mutators:
  int getNumberToResubmit(TString jobType);
  std::vector<TString> getResubmitList(TString jobType);
  updateJobStatus(TString jobType);
  
  // Accessors:
  printResubmitList(TString jobType);
  
 private:
    
  TString jobName;
  std::vector<TString> listDMWorkspace;
  std::vector<TString> listDMTestStat;
  
};

#endif
