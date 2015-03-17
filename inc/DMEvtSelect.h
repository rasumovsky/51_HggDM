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

// C++ includes:
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

// ROOT includes:
#include "TString.h"

// Package includes:
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
  void printCutflow(bool weighted);
  void saveCutflow(TString filename, bool weighted);
  
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
  bool recursiveCall;
  
};

#endif
