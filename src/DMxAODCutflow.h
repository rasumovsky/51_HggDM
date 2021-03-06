////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMxAODCutflow.h                                                     //
//  Class: DMxAODCutflow.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMxAODCutflow_h
#define DMxAODCutflow_h

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
#include "TFile.h"
#include "TH1F.h"

// Package includes:
#include "Config.h"
#include "DMAnalysis.h"

class DMxAODCutflow 
{
  
 public:
  
  //DMxAODCutflow();
  DMxAODCutflow(TString fileName, TString configFileName);
  virtual ~DMxAODCutflow() {};
  
  // Accessors:
  TString getCutNameByOrder(int order);
  int getCutOrderByName(TString cutName);
  double getEventsPassingCut(int order);
  double getEventsPassingCut(TString cutName);
  TH1F* getHist();
  double getPercentPassingCut(TString cutName);
  double getAccXEffAtCut(TString cutName);
  void printxAODCutflow();
  double nTotalEventsInFile();
  
 private:
  
  // Member methods:
  bool cutExists(TString cutName);
  
  // Member objects:
  int nCuts;
  std::vector<TString> cutList;
  std::map<TString,double> passCounter_weighted;
  std::map<TString,double> passCounter_unweighted;
  
  bool m_unskimmed;
  bool m_isWeighted;
  TFile *m_inputFile;
  TH1F *m_histCuts_weighted;
  TH1F *m_histCuts_unweighted;
  //TH1F *m_histCuts_weighted_NoDalitz;
  
};

#endif
