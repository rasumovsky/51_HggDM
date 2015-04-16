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

#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "CommonFunc.h"
#include "statistics.h"
#include "ESSReader.h"
#include "ResReader.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "DMSigParam.h"
#include "DMBkgModel.h"

class DMWorkspace
{

 public:
  
  DMWorkspace(TString newJobName, TString newCateScheme, TString newOptions);
  virtual ~DMWorkspace() {};
  
 private:
  
  bool fitsAllConverged();
  RooWorkspace* DMWorkspace::getCombinedWorkspace();
  ModelConfig* DMWorkspace::getModelConfig();
  void loadWSFromFile();
  void createNewWS();
  RooWorkspace* createNewCategoryWS();
  void createAsimovData(RooWorkspace *cateWS, RooDataSet *obsData, 
			RooRealVar wt, int valueMuDM);
  // eliminate this as soon as possible:
  double spuriousSignal(TString cateName);
  void makeNP(const char* varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(const char* varnameNP, const char* process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  void plotFit(RooWorkspace *cateWS, double valMuDM);
  //void plotNuisParams();// Take this from NPP?
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  // Helper classes:
  ESSReader *ess_tool;
  ResReader* res_tool;
  DMEvtSelect *selector;
  DMSigParam *currSigParam;
  DMBkgModel *currBkgModel;
  DMMassPoints *currMassPoints;
  
  // Updated for each call to createNewCategoryWS():
  int currCateIndex;
  TString currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *combinedWS;
  ModelConfig *mconfig;
  
  bool allGoodFits;
  
};

#endif

