////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMTestStat.h                                                        //
//  Class: DMTestStat.cxx                                                     //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 19/04/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0 and CLs based on an input      //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMTestStat_h
#define DMTestStat_h

#include "CommonHead.h"
#include "CommonFunc.h"
#include "RooFitHead.h"
#include "statistics.h"
#include "DMWorkspace.h"

class DMTestStat {
  
 public:
    
  DMTestStat(TString newJobName, TString newDMSignal, TString newCateScheme, 
	     TString newOptions);
  DMTestStat(TString newJobName, TString newDMSignal, TString newOptions,
	     RooWorkspace *newWorkspace);
  virtual ~DMTestStat() {};
  
  double accessValue(TString testStat, bool observed, int N);
  void calculateNewCL();
  void calculateNewP0();
  void clearData();
  bool fitsAllConverged();
  double getCLsFromCL(double CL);
  double getCLsFromQMu(double qMu, TString type);
  double getCLFromQMu(double qMu, TString type);
  double getQ0FromNLL(double nllMu0, double nllMuHat, double muHat);
  double getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
		       double muTest);
  double getP0FromQ0(double q0);
  double getPMuFromQMu(double qMu);
  double getPbFromQMu(double qMu, double sigma, double mu);
  double getPbfromN(double N);
  
 private:

  double getFitNLL(TString datasetName, double muVal, bool fixMu,
		   double &profiledMu);
  void loadStatsFromFile();
  bool mapValueExists(TString mapKey);
  
  // From the initialization:
  TString jobName;
  TString DMSignal;
  TString options;
  TString outputDir;
  
  // Check whether all fits successful:
  bool allGoodFits;
  
  // The workspace for the fits:
  RooWorkspace *workspace;
  
  // Store the calculated values:
  std::map<TString,double> calculatedValues;
  
};

#endif

