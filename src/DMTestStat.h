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
	     TString newOptions, RooWorkspace *newWorkspace);
  virtual ~DMTestStat() {};
  
  double accessValue(TString testStat, bool observed, int N);
  void calculateNewCL();
  void calculateNewP0();
  void clearData();
  bool fitsAllConverged();
  double functionQMu(double x);
  double functionQMuTilde(double x, double asimovTestStat);
  double getCLFromCLs(double CLs);
  double getCLsFromCL(double CL);
  double getCLFromQMu(double qMu, bool observed, double N);
  double getCLsFromQMu(double qMu, bool observed, double N);
  double getFitNLL(TString datasetName, double muVal, bool fixMu,
		   double &profiledMu);
  std::vector<std::string> getGlobsNames();
  std::vector<double> getGlobsValues();
  std::vector<std::string> getNPNames();
  std::vector<double> getNPValues();
  double getP0FromQ0(double q0);
  double getPbFromN(double N);
  double getPbFromQMu(double qMu, double sigma, double mu);
  double getPMuFromQMu(double qMu);
  double getQ0FromNLL(double nllMu0, double nllMuHat, double muHat);
  double getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
		       double muTest);
  double getQMuTildeFromNLL(double nllMu, double nllMu0, double nllMuHat,
			    double muHat, double muTest);
  void loadStatsFromFile();
  
 private:
  
  TString getKey(TString testStat, bool observed, int N);
  bool mapValueExists(TString mapKey);
  
  // From the initialization:
  TString jobName;
  TString DMSignal;
  TString options;
  TString outputDir;
  
  TString dataForObs;
  TString dataForExp;

  // Check whether all fits successful:
  bool allGoodFits;
  
  // The workspace for the fits:
  DMWorkspace *dmw;
  RooWorkspace *workspace;
  ModelConfig *mc;

  // Store the calculated values:
  std::map<TString,double> calculatedValues;
  
  // Store fit parameters from NLL calculation:
  std::vector<std::string> namesGlobs;
  std::vector<std::string> namesNP;
  std::vector<double> valuesGlobs;
  std::vector<double> valuesNP;
  
};

#endif

