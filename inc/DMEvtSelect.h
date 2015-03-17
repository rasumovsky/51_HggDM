////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMEvtSelect.h                                                       //
//  Class: DMEvtSelect.cxx                                                    //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 16/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMEvtSelect_h
#define DMEvtSelect_h

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <map>
#include "DMTree.h"



class DMEvtSelect 
{
  
 public:
  
  DMEvtSelect(DMTree *newTree);
  ~DMEvtSelect();
  
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
  DMTree *evtTree;
  std::vector<TString> cutNames;
  std::map<TString,int> evtCountPass;
  std::map<TString,double> evtCountPassWt;
  std::map<TString,int> evtCountTot;
  std::map<TString,double> evtCountTotalWt;
  
};

#endif
