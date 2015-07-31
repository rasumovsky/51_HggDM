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

#include "BkgModel.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "CommonFunc.h"
#include "statistics.h"
#include "PESReader.h"
#include "PERReader.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "DMTestStat.h"
#include "SigParam.h"
#include "SigParamInterface.h"

class DMWorkspace {

 public:
  
  DMWorkspace(TString newJobName, TString newDMSignal, TString newCateScheme, 
	      TString newOptions);
  virtual ~DMWorkspace() {};
  
  bool fitsAllConverged();
  RooWorkspace* getCombinedWorkspace();
  ModelConfig* getModelConfig();
  
 private:
  
  void loadWSFromFile();
  void createNewWS();
  RooWorkspace* createNewCategoryWS();
  //void createAsimovData(RooWorkspace *cateWS, int valueMuDM, int valueMuSM);
  double spuriousSignal();// eliminate this as soon as possible
  void makeNP(TString varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(TString varnameNP, TString process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  //void plotSingleCateFit(RooWorkspace *cateWS, TString dataset);
  //void plotFinalFits(RooWorkspace *combWS, TString fitType);
  //void plotNuisParams(RooArgSet nuisParams, TString type);
  //void profileAndSnapshot(TString muDMValue, double &nllValue,
  //			  double &profiledMu);

  // Member variables:
  TString m_jobName;
  TString m_DMSignal;
  TString m_cateScheme;
  TString m_options;
  TString m_outputDir;

  int m_nCategories;
  int m_muNominalSM;
  TString m_dataToPlot;

  // Helper classes:
  PESReader *m_pes;
  PERReader* m_per;
  //DMSigParam *currSigParam;
  SigParamInterface *m_spi;

  // Updated for each call to createNewCategoryWS():
  int m_currCateIndex;
  TString m_currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *m_combinedWS;
  ModelConfig *m_mConfig;
  
  bool m_allGoodFits;
 
};

#endif

