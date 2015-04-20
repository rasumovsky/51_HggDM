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
#include "DMMassPoints.h"
#include "DMSigParam.h"
#include "DMBkgModel.h"
#include "DMTestStat.h"

void makeExe(TString exeName);

void submitWSViaBsub(TString exeJobName, TString exeOption, TString exeSignal,
		     TString exeCateScheme);

void submitTSViaBsub(TString exeJobName, TString exeOption, TString exeSignal);
