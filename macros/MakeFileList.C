//////////////////////////////////////////
//                                      //
//  GetResubmitList.C                   //
//                                      //
//  Created: Andrew Hard 15/10/2012     //
//                                      //
//////////////////////////////////////////

#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TLatex.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

/**
   -----------------------------------------------------------------------------
   Transfers a single MxAOD and cutflow file and creates a file list.
   @param inputDir - The input directory containing the original MxAOD files.
   @param sampleName - The name of the current sample.
   @param date - The date to assign to the MxAOD file batch.
*/
void SingleTransfer(TString inputDir, TString sampleName, TString date) {
    
  // Define output locations:
  TString cutFlowDir = Form("/afs/cern.ch/user/a/ahard/work_directory/files_HggDM/GlobalInputs/xAODCutFlows/%s", date.Data());
  TString fileListDir = Form("/afs/cern.ch/user/a/ahard/work_directory/files_HggDM/GlobalInputs/FileLists/%s", date.Data());
  TString eosDir = Form("/eos/atlas/user/a/ahard/HggDM_13TeV/%s", date.Data());
  
  // Make sure output directories exist:
  system(Form("mkdir -vp %s", cutFlowDir.Data()));
  system(Form("mkdir -vp %s", fileListDir.Data()));
  system(Form("/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select mkdir %s", eosDir.Data()));
  
  // Copy cut-flow histogram to AFS:
  system(Form("cp %s/%s/hist-sample.root %s/hist_%s.root", inputDir.Data(),
	      sampleName.Data(), cutFlowDir.Data(), sampleName.Data())); 
  
  // Create a file list with the EOS location:
  ofstream outputFile(Form("%s/list_%s.txt", fileListDir.Data(), sampleName.Data()));
  outputFile << Form("root://eosatlas/%s/ntuple_%s.root", eosDir.Data(), sampleName.Data()) << std::endl;
  outputFile.close();
  
  // Copy the MxAOD to EOS:
  system(Form("xrdcp %s/%s/data-MxAOD/sample.root root://eosatlas/%s/ntuple_%s.root", inputDir.Data(), sampleName.Data(), eosDir.Data(), sampleName.Data()));
}

/**
   -----------------------------------------------------------------------------
   Main method copies files from an input directory to afs and eos and creates
   file lists for the MxAODs and a list of signal samples.
   @param inputDir - The input directory containing the original MxAOD files.
   @param date - The date to assign to the MxAOD file batch.
*/
void MakeFileList(TString inputDir, TString date) {
  system(Form("ls %s | tee sampleList.txt", inputDir.Data()));
  ofstream outputSampleList("newSampleList.txt");
  ifstream sampleFile;
  sampleFile.open("sampleList.txt");
  TString currSample;
  while (!sampleFile.eof()) {
    sampleFile >> currSample;
    if (currSample.Contains("H2yyMETAnalysis_3") || currSample.Contains("H2yyMETAnalysis_4")) {
      SingleTransfer(inputDir, currSample, date);
      outputSampleList << currSample << std::endl;
    }
  }
  sampleFile.close();
  system(Form("rm sampleList.txt"));
  outputSampleList.close();
}
