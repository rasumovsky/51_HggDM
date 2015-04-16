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

// Systematic Uncertainty readers (DEPRECATED):
#include "ESSReader.h"
#include "ResReader.h"
#include "DMAnalysis.h"

class DMWorkspace
{

 public:
  
  DMWorkspace(TString newJobName, TString newCateScheme, TString newOptions);
  virtual ~DMWorkspace() {};
  
 private:
  
  void backgroundPdfBuilder(RooWorkspace *&w, RooArgSet *&nuispara,
		       TString cateName);
  void createAsimovData(TString cateName, RooWorkspace *currWS, 
			RooDataSet *currData, RooRealVar currWeightVar, 
			double valuePOI);
  RooWorkspace* createNewCategoryWS();
  void createNewWS();
  void loadWSFromFile();
  
  void makeNP(const char* varname, double setup[5], RooArgSet *&nuispara,
	      RooArgSet *&constraints, RooArgSet *&globobs,
	      RooArgSet *&expected);
  void makeShapeNP(const char* varnameNP, const char* proc, double setup[5],
		   RooArgSet *&nuispara, RooArgSet *&constraints,
		   RooArgSet *&globobs, RooArgSet *&expected);
  double spuriousSignal(TString cateName);
  
  void plotFit(TString cateName, RooWorkspace *workspace);
  void plotNuisParams();
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  TString outputDir;
  
  // Helper classes:
  DMEvtSelect *selector;
  
  //ESSReader *ess_tool;
  //ResReader* res_tool;
  
  // The RooWorkspace and ModelConfig:
  RooWorkspace *combinedWS;
  ModelConfig *mconfig;
  
  // Updated inside each call to createNewCategoryWS():
  int currCateIndex;
  TString currCateName;
  RooWorkspace *currWS;
  DMSigParam *currSigParam;
  DMBkgModel *currBkgModel;
  DMMassPoints *currMassPoints;
};

#endif

