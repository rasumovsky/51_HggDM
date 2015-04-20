/////////////////////////////////////////
//                                     //
//  NPPV1_checkjobs.hh                 //
//                                     //
//  Author: Andrew Hard                //
//  Date: 26/11/2013                   //
//                                     //
/////////////////////////////////////////

#ifndef DMCheckJobs_h
#define DMCheckJobs_h

#include "DMAnalysis.h"

class DMCheckJobs {

 public:
  
  // Constructor and destructor:
  DMCheckJobs(TString newJobName);
  virtual ~DMBkgModel() {};
  
  // Accessors:
  getResubmitList(TString jobType);
  
 private:
  
  TString jobName;
  
  
};

#endif
