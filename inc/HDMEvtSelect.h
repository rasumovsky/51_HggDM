////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: HDMEvtSelect.h                                                      //
//  Class: HDMEvtSelect.cxx                                                   //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 16/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef HDMEvtSelect_h
#define HDMEvtSelect_h

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <map>
#include "HDMTree.h"



class HDMEvtSelect 
{
  
 public:
  
  HDMEvtSelect(HDMTree *newTree);
  ~HDMEvtSelect();
  
  // Accessors:
  int getPassingEvents(TString cutname);
  double getPassingEventsWt(TString cutname);
  int getTotalEvents(TString cutname);
  double getTotalEventsWt(TString cutname);
  void printCuts();
  
  // Mutators
  void clearCounters();
  bool passesCut(TString cutname);
  bool passesCut(TString cutname, double weight);
  
 private:
  
  bool cutExists();
  
  // Member objects:
  HDMTree *evtTree;
  std::vector<TString> cutNames;
  std::map<TString,int> evtCountPass;
  std::map<TString,double> evtCountPassWt;
  std::map<TString,int> evtCountTot;
  std::map<TString,double> evtCountTotalWt;
  
};

#endif
