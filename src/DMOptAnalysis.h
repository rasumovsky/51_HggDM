////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DMOptAnalysis.h                                                           //
//  Class: DMOptAnalysis.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 21/08/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMOptAnalysis_h
#define DMOptAnalysis_h

// Package libraries:
#include "AnaCollection.h"
#include "AnaInfo.h"
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DMAnalysis.h"


class DMOptAnalysis {

 public:
  
  // Class constructor and destructor:
  DMOptAnalysis(TString newConfigFile);
  virtual ~DMOptAnalysis() {};
  
  // Public accessors:
  std::vector<TString> listDirectoryContents(TString directory);
  void plot2DScatter(TString quantity1, TString quantity2);
  void plotCutsAndStat(TString signal, TString cutNameX, TString cutNameY,
		       TString statistic);
  void plotOptimizationPoints(TString cutName1, TString cutName2);
  
  // Public mutators:
  void loadOptimizationData(TString directory);
  
 private:
  
  // Private accessors:
  std::vector<double> checkDoubleList(std::vector<double> currList,
				      double newDouble);
  void getHistBinsAndRange(TString cutName, int &bins, double &min,
			   double &max);
  double maxEntry(std::vector<double> currList);
  double minEntry(std::vector<double> currList);
  std::vector<TString> vectorizeTString(TString originString, TString delim);
  
  // Private mutators:
  
  
  
  // Private member variables:
  TString m_outputDir;
  Config *m_config;
  AnaCollection *m_anaCollection;
  std::vector<TString> m_optimizedCutList;
  
};

#endif

