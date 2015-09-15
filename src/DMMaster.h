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
#include "DMOptAnalysis.h"
#include "DMTestStat.h"
#include "DMToyAnalysis.h"
#include "SigParamInterface.h"

// Analysis configuration file:
Config *m_config;

// Two items below are for recursive job submission:
int m_jobIndex;
ofstream m_headFile;

bool m_isFirstJob;

// Mutators:
void recursiveOptimizer(TString exeConfigOrigin, TString exeOption, 
			int cutIndex, std::vector<TString> cutN,
			std::vector<double> cutV);
void submitToOptimize(TString exeConfigOrigin, TString exeOption);

void submitWSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitMLViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitPEViaBsub(TString exeConfigFile, TString exeOption, 
		     TString exeSignal, int exeSeed, int exeToysPerJob);
