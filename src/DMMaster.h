////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <time.h>

// Package libraries:
#include "BkgModel.h"
#include "Config.h"
#include "DMAnalysis.h"
#include "DMCheckJobs.h"
#include "DMMassPoints.h"
#include "DMTestStat.h"
#include "DMToyAnalysis.h"
#include "SigParamInterface.h"

bool m_isFirstJob;

Config *m_config;

void submitToOptimize(TString exeConfigOrigin, TString exeOption);

void submitWSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitMLViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitPEViaBsub(TString exeConfigFile, TString exeOption, 
		     TString exeSignal, int exeSeed, int exeToysPerJob);
