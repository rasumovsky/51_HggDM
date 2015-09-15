////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  AnaCollection.h                                                           //
//  Class: AnaCollection.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AnaCollection_h
#define AnaCollection_h

// Package libraries:
#include "AnaInfo.h"

class AnaCollection {

 public:
  
  // Class constructor and destructor:
  AnaCollection();
  virtual ~AnaCollection() {};
  
  // Public accessors:
  AnaInfo* getOptimalAnalysis(TString signal, TString statistic);
  int nAnalyses();
  int nBadAnalyses();
  std::vector<AnaInfo*>::iterator begin();
  std::vector<AnaInfo*>::iterator end();
  
  // Public mutators:
  void addAnalysis(AnaInfo *newAnalysis);
    
 private:
  
  // Private member variables:
  std::vector<AnaInfo*> m_analyses;
  
};

#endif

