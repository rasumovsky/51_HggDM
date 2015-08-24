////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: August 6, 2015                                                      //
//                                                                            //
//  This program just runs over data files from 7 and 8 TeV and produces a    //
//  smaller local file that can be used.                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
//#include "TROOT.h"
#include "TTree.h"

#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cassert>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
#include <map>
#include <cstdlib>
#include <cmath>

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void PrintProgressBar(int index, int total) {
  if (index%10000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
} 

/**
   The main method requires a file list to run.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 2) {
    std::cout << "\nUsage: " << argv[0] << " <fileList>" << std::endl;
    exit(0);
  }
  
  TString inputFileName = argv[1];
  // Declare TChain and add input files:
  TChain *chain = new TChain("Hgg_tree");
  TString currFileName;
  std::ifstream inputFile;
  inputFile.open(inputFileName);
  while (!inputFile.eof()) {
    inputFile >> currFileName;
    std::cout << "\t" << currFileName << std::endl;
    if (currFileName.Contains("data11")) {
      chain->AddFile(Form("root://eosatlas//eos/atlas/user/f/fwang/HSG1_ntuples_TwikiApr26/Data11/%s", currFileName.Data()));
    }
    else if (currFileName.Contains("data12")) {
      chain->AddFile(Form("root://eosatlas//eos/atlas/user/f/fwang/HSG1_ntuples_TwikiApr26/Data12/%s", currFileName.Data()));
    }
  }
  inputFile.close();
  Int_t nEvents = chain->GetEntries();
  std::cout << " We have " << nEvents << " events to be processed" << std::endl;
  
  // Set the relevant branch addresses:
  double mass_PV_EM; int category_MoriondMVA; bool flag_all; int Run;
  chain->SetBranchAddress("mass_PV_EM", &mass_PV_EM);
  //chain->SetBranchAddress("category_MoriondMVA", &category_MoriondMVA);
  chain->SetBranchAddress("flag_all", &flag_all);
  chain->SetBranchAddress("Run", &Run);
  
  // Create output TFile and TTree:
  TFile *outputFile = new TFile("output.root", "RECREATE");
  TTree *outputTree = new TTree("Hgg_tree", "Hgg_tree");
  outputTree->Branch("mass_PV_EM", &mass_PV_EM, "mass_PV_EM/D");
  //outputTree->Branch("category_MoriondMVA", &category_MoriondMVA, 
  //		     "category_MoriondMVA/I");
  outputTree->Branch("flag_all", &flag_all, "flag_all/O");
  outputTree->Branch("Run", &Run, "Run/I");
  
  // Loop over events
  for( int index = 0; index < nEvents; index++ ) {
    chain->GetEntry(index);
    PrintProgressBar(index, nEvents);
    if (!flag_all) continue;
    outputTree->Fill();
  }
  
  // Write everything to file:
  std::cout << "Finished event loop" << std::endl;
  outputTree->Write();
  outputFile->Close();
  
  return 0;
}
