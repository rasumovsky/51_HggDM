////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMAnalysis.h"
#include "BkgModel.h"
#include "DMCheckJobs.h"
#include "DMMassPoints.h"
#include "DMSigParam.h"
#include "DMTestStat.h"
#include "DMToyAnalysis.h"

bool isFirstJob;

void makeExe(TString exeName);

void submitWSViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme);

void submitTSViaBsub(TString exeJobName, TString exeOption, TString exeSignal);

void submitMLViaBsub(TString exeJobName, TString exeOption, TString exeSignal);

void submitPEViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     int exeSeed, int exeToysPerJob);
