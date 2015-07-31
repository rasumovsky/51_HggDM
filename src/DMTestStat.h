////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMTestStat.h                                                        //
//  Class: DMTestStat.cxx                                                     //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 30/07/2015                                                          //
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
  void clearFitParamSettings();
  RooDataSet* createAsimovData(int valMuDH, int valMuSH);
  RooDataSet* createAsimovData(TString datasetName);
  RooDataSet* createPseudoData(int seed, int valMuDM, int valMuSM, bool fixMu);
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
  void saveSnapshots(bool doSaveSnapshot);
  void setPlotDirectory(TString directory);
  void setParams(TString paramName, double paramVal, bool doSetConstant);
  
 private:
  
  TString getKey(TString testStat, bool observed, int N);
  bool mapValueExists(TString mapKey);
  void plotFits(TString fitType, TString datasetName);
  
  // From the initialization:
  TString m_jobName;    // The name of the group of jobs (for I/O purposes).
  TString m_DMSignal;   // The specific DM signal under study.
  TString m_cateScheme; // The categorization scheme for the analysis.
  TString m_options;    // Job options.
  TString m_outputDir;  // The output directory for statistics.
  
  TString m_dataForObs; // The dataset for observed results.
  TString m_dataForExp; // The dataset for expected results.
  TString m_plotDir;    // The output directory for plots (not set by default).
  bool m_doSaveSnapshot;// Option to save snapshots of PoI and nuisance params.
  bool m_doPlot;        // Sets whether or not to plot fit results.

  // Pointer to the input file:
  TFile *m_inputFile;
  
  // Check whether all fits successful:
  bool m_allGoodFits;
  
  // The workspace for the fits:
  RooWorkspace *m_workspace;
  ModelConfig *m_mc;

  // Store the calculated values:
  std::map<TString,double> m_calculatedValues;
  
  // Store fit parameters from NLL calculation:
  std::vector<std::string> m_namesGlobs;
  std::vector<std::string> m_namesNP;
  std::vector<double> m_valuesGlobs;
  std::vector<double> m_valuesNP;
  
  // In case special parameter settings are used for a fit:
  std::vector<bool> m_setParamConsts;
  std::vector<TString> m_setParamNames;
  std::vector<double> m_setParamVals;
  
};

#endif

