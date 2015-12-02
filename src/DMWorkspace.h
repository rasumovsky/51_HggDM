////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMWorkspace.h                                                       //
//  Class: DMWorkspace.cxx                                                    //
//  Creator: Andrew Hard, Hongtao Yang, Haichen Wang                          //
//  Email: ahard@cern.ch                                                      //
//  Date: 02/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef DMWorkspace_h
#define DMWorkspace_h

// Package libraries:
#include "BkgModel.h"
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "DMTestStat.h"
#include "HggTwoSidedCBPdf.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "SigParam.h"
#include "SigParamInterface.h"
#include "statistics.h"

class DMWorkspace {

 public:
  
  DMWorkspace(TString newConfigFile, TString newDMSignal, TString newOptions);
  virtual ~DMWorkspace() {};
  
  bool fitsAllConverged();
  RooWorkspace* getCombinedWorkspace();
  ModelConfig* getModelConfig();
  
 private:
  
  void loadWSFromFile();
  void createNewWS();
  RooWorkspace* createNewCategoryWS();
  double spuriousSignal();// eliminate this as soon as possible
  void makeNP(TString varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(TString varnameNP, TString process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  void createAsimovData(RooWorkspace *cateWS, int valMuDM, int valMuSM);

  // Member variables:
  TString m_configFile;
  TString m_DMSignal;
  TString m_options;
  TString m_outputDir;

  int m_nCategories;
  int m_muNominalSM;
  TString m_dataToPlot;

  // Helper classes:
  Config *m_config;
  SigParamInterface *m_spi;

  // Updated for each call to createNewCategoryWS():
  int m_currCateIndex;
  TString m_currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *m_combinedWS;
  ModelConfig *m_modelConfig;
  
  bool m_allGoodFits;
 
};

#endif

