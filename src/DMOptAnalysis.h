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
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DMAnalysis.h"

class DMOptAnalysis {

 public:
  
  // Class constructor and destructor:
  DMOptAnalysis(TString newConfigFile);
  virtual ~DMOptAnalysis() {};
  
  // A data structure to store the important analysis attributes:
  struct AnalysisAttributes {
    int index; // Job index.
    bool isGood; // Data loaded successfully?
    std::map<TString,double> cutValues; // Cut names and their values.
    std::vector<TString> signals; // A list of signal names
    std::map<TString,double> statValues; // statistics for ALL signals.
  };
  
  // Public accessors:
  std::vector<TString> listDirectoryContents(TString directory);
  void plot2DScatter(TString quantity1, TString quantity2);
  void plotCutsAndStat(TString signal, TString cutNameX, TString cutNameY,
		       TString statistic, bool minimize);
  void plotOptimizationPoints(TString cutName1, TString cutName2);
  
  // Public mutators:
  void loadOptimizationData(TString directory);
  
 private:
  
  // Private accessors:
  std::vector<double> checkDoubleList(std::vector<double> currList,
				      double newDouble);
  void getHistBinsAndRange(TString cutName, int &bins, double &min,
			   double &max);
  int getOptAnaIndex(TString signal, TString statistic, bool minimize);
  TString mapKey(TString signal, TString statistic);
  double maxEntry(std::vector<double> currList);
  double minEntry(std::vector<double> currList);
  
  // Private mutators:
  
  

  // Private member variables:
  TString m_outputDir;
  Config *m_config;
  
  std::vector<AnalysisAttributes*> m_analysisList;
  
};

#endif

