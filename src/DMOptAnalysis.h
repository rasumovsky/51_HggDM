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

class DMOptAnalysis {

 public:
  
  // Class constructor and destructor:
  DMOptAnalysis(TString newConfigFile);
  virtual ~DMOptAnalysis() {};
  
  // A data structure to store the important analysis attributes:
  struct AnalysisAttributes {
    int index; // Job index.
    bool isGood; // Data loaded successfully?
    std::map<TString,double> cutNameAndVal; // Cut names and their values.
    std::vector<TString> signals; // A list of signal names
    std::map<TString,double> valuesExpCLN2; // Map of signal to exp CL (-2sigma)
    std::map<TString,double> valuesExpCLN1; // Map of signal to exp CL (-1sigma)
    std::map<TString,double> valuesExpCL;   // Map of signal to exp CL (nom)
    std::map<TString,double> valuesExpCLP1; // Map of signal to exp CL (+1sigma)
    std::map<TString,double> valuesExpCLP2; // Map of signal to exp CL (+2sigma)
    std::map<TString,double> valuesObsCL;   // Map of signal to obs CL
    std::map<TString,double> valuesExpP0;   // Map of signal to exp p0
    std::map<TString,double> valuesObsP0;   // Map of signal to obs p0
  };
  
  // Public accessors:
  std::vector<TString> listDirectoryContents(TString directory);
  
  // Public mutators:
  void plotOptimizationPoints(TString cutName1, TString cutName2);
  
 private:
  
  // Private accessors:
  TString getPrintName(TString originName);
  
  // Private mutators:
  void loadOptimizationData(TString directory);
  

  // Private member variables:
  TString m_outputDir;
  Config *m_config;
  
  std::vector<AnalysisAttributes*> analysisList;
  
  std::map<TString,std::vector<double> > cutValues;

};

#endif

