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
#include "PESReader.h"
#include "PERReader.h"
#include "DMAnalysis.h"
#include "DMMassPoints.h"
#include "DMSigParam.h"
#include "DMBkgModel.h"

class DMWorkspace
{

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
  void createAsimovData(RooWorkspace *cateWS, RooDataSet *obsData, 
			RooRealVar wt, int valueMuDM);
  // eliminate this as soon as possible:
  double spuriousSignal();
  void makeNP(TString varName, double setup[4], RooArgSet *&nuisParams,
	      RooArgSet *&constraints, RooArgSet *&globalObs,
	      RooArgSet *&expected);
  void makeShapeNP(TString varnameNP, TString process, double setup[4],
		   RooArgSet *&nuisParams, RooArgSet *&constraints,
		   RooArgSet *&globalObs, RooArgSet *&expected);
  void plotFit(RooWorkspace *cateWS, double valMuDM);
  //void plotNuisParams();// Take this from NPP?
  
  // Member variables:
  TString jobName;
  TString DMSignal;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  // Helper classes:
  PESReader *pes;
  PERReader* per;
  DMEvtSelect *selector;
  DMSigParam *currSigParam;
  DMBkgModel *currBkgModel;
  DMMassPoints *currMassPoints;
  
  // Updated for each call to createNewCategoryWS():
  int currCateIndex;
  TString currCateName;
  
  // The Final RooWorkspace and ModelConfig:
  RooWorkspace *combinedWS;
  ModelConfig *mConfig;
  
  bool allGoodFits;
  
};

#endif

