////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMCutflowPlotter.cxx                                                //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
//  This macro plots the cutflows of several processes together.              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// ROOT libraries:
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TString.h"
#include "Rtypes.h"

// Package libraries:
#include "Config.h"
#include "CommonFunc.h"
#include "CommonHead.h"
#include "DMAnalysis.h"
#include "DMxAODCutflow.h"

int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0]
	      << " <configFile> <options>" << std::endl;
    exit(0);
  }
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("jobName");
  TString cateScheme = config->getStr("cateScheme");

  // Tool to get the total number of events at the generator level.
  DMxAODCutflow *dmx
    = new DMxAODCutflow(DMAnalysis::nameToxAODCutFile(m_config, m_sampleName));
  double nGeneratedEvt = dmx->getEventsPassingCut(1);
  
  
  return 0;
}
