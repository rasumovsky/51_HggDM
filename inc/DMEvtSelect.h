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
#include "TH1F.h"

// Package includes:
#include "Config.h"
#include "DMTree.h"

using std::vector;

class DMEvtSelect 
{
  
 public:
  
  //DMEvtSelect();
  DMEvtSelect(DMTree *newTree, TString newConfigFile);
  virtual ~DMEvtSelect() {};
  
  // Public Accessors:
  int cutIndex(TString cutName);
  int getEventsPerCate(TString cateScheme, int cate);
  double getEventsPerCateWt(TString cateScheme, int cate);
  int getPassingEvents(TString cutName);
  double getPassingEventsWt(TString cutName);
  int getTotalEvents(TString cutName);
  double getTotalEventsWt(TString cutName);
  int nCuts();
  void printCutflow(bool weighted);
  void printCategorization(bool weighted);
  TH1F* retrieveCutflowHist(bool weighted);
  void saveCutflow(TString filename, bool weighted);
  void saveCategorization(TString filename, bool weighted);

  // Public Mutators
  void clearCounters();
  int getCategoryNumber(TString cateScheme);
  int getCategoryNumber(TString cateScheme, double weight);
  bool passesCut(TString cutName);
  bool passesCut(TString cutName, double weight);
  void setTree(DMTree *newTree);

 private:
  
  // Member methods:
  bool cutExists(TString cutName);
  bool cateExists(TString cateScheme);
  
  // Member objects:
  Config *m_config;
  DMTree *evtTree;
  std::vector<TString> cutList;
  
  std::map<TString,int> evtCountPass;
  std::map<TString,double> evtCountPassWt;
  std::map<TString,int> evtCountTot;
  std::map<TString,double> evtCountTotWt;
  
  std::map<TString,int> cateSchemesAndSizes;
  std::map<TString,int> cateCount;
  std::map<TString,double> cateCountWt;
  
  // Cutflow histogram:
  TH1F *m_cutFlowHist_weighted;
  TH1F *m_cutFlowHist_unweighted;
  
};

#endif
