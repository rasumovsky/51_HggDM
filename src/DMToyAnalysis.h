////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DMToyAnalysis.h                                                           //
//  Class: DMToyAnalysis.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 24/06/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMToyAnalysis_h
#define DMToyAnalysis_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DMToyTree.h"
#include "DMTestStat.h"
#include "RooFitHead.h"
#include "statistics.h"

class DMToyAnalysis {

 public:
  
  DMToyAnalysis(TString newConfigFile, TString newDMSignal);
  virtual ~DMToyAnalysis() {};
  
  void fillToyHistograms(int muValue, DMToyTree *toyTree);
  void getAsymptoticForm(TString statistic);
  TH1F* getAsymptoticHist();
  TH1F* getGlobsHist(TString paramName, TString fitType, int toyMu);
  TH1F* getNuisHist(TString paramName, TString fitType, int toyMu);
  TH1F* getMuHist(int toyMu);
  TH1F* getStatHist(TString statistic, int toyMu);
  void plotParameter(TString paramName, TString paramType, int toyMu);
  void plotProfiledMu(); 
  void plotTestStat(TString statistic);
  void plotTestStatComparison(TString statistic);
    
 private:
  
  TString printStatName(TString statistic);
  
  // Private member variables:
  TString m_outputDir;
  
  // Classes for statistics access:
  DMTestStat *m_dmts;
  RooWorkspace *m_workspace;
  
  // Test statistic binning:
  int m_nBins;
  int m_binMin;
  int m_binMax;
  
  // Histograms:
  TH1F *m_hAsymptotic;
  TH1F *m_hMuProfiled[2];
  TH1F *m_hQ0[2];
  TH1F *m_hQMu[2];
  TH1F *m_hQMuTilde[2];
  TH1F *m_hNuisMu0[20][2];
  TH1F *m_hNuisMu1[20][2];
  TH1F *m_hNuisMuFree[20][2];
  TH1F *m_hGlobsMu0[20][2];
  TH1F *m_hGlobsMu1[20][2];
  TH1F *m_hGlobsMuFree[20][2];
  
  // Parameter data:
  std::vector<std::string> m_namesGlobs;
  std::vector<std::string> m_namesNuis;
  int m_nGlobs;
  int m_nNuis;
  
};

#endif

