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
#include "statistics.hh"

// Systematic Uncertainty readers (DEPRECATED):
#include "ESSReader.h"
#include "ResReader.h"

#include "DMMaster.hh"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;

class DMWorkspace
{

 public:
  
  DMWorkspace(TString newJobName, TString newCateScheme, TString newOptions);
  virtual ~DMWorkspace() {};
  
 private:
  
  void backgroundPdfBuilder(RooWorkspace *&w, RooArgSet *&nuispara,
		       TString channelname);
  RooWorkspace* newCategoryWS(TString channelName);
  void createNewWS();
  void loadWSFromFile();
  
  void makeNP(const char* varname, double setup[5], RooArgSet *&nuispara,
	      RooArgSet *&constraints, RooArgSet *&globobs,
	      RooArgSet *&expected);
  void makeShapeNP(const char* varnameNP, const char* proc, double setup[5],
		   RooArgSet *&nuispara, RooArgSet *&constraints,
		   RooArgSet *&globobs, RooArgSet *&expected);
  double spuriousSignal(TString channelname);
  
  void plotFit(TString plotOptions);
  void plotNuisParams(TString plotOptions);
  
  // Member variables:
  TString jobName;
  TString sampleName;
  TString cateScheme;
  TString options;
  
  TString outputDir;
  
  // Helper classes:
  DMEvtSelect *selector;
  DMSigParam *sigParam;
  DMMassPoints *massPoints;
  DMBkgModel *bkgModel;
  //ESSReader *ess_tool;
  //ResReader* res_tool;
  
  // The RooWorkspace and ModelConfig:
  RooWorkspace* combination;
  ModelConfig *mconfig;
  
};

#endif

