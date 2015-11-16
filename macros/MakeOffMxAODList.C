//////////////////////////////////////////
//                                      //
//  MakeOffMxAODList.C                  //
//                                      //
//  Created: Andrew Hard 12/10/2015     //
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
   Create file lists for the centrally produced MxAODs.
   @param inputDir - The input directory containing the original MxAOD files.
   @param tag - The MxAOD tag (h008, etc...).
   @param date - The date to assign to the MxAOD file batch.
*/
void MakeOffMxAODList(TString inputDir, TString inputList, TString tag,
		      TString date) {
  
  // Define and create output directory for file lists:
  TString outputDir = Form("/afs/cern.ch/user/a/ahard/work_directory/files_HggDM/GlobalInputs/FileLists/%s",date.Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Input file listing the file names:
  ifstream sampleFile;
  sampleFile.open(inputList);
  TString currSample;
  while (!sampleFile.eof()) {
    sampleFile >> currSample;
    TString newFileName = currSample;
    newFileName = newFileName.ReplaceAll(".root",".txt");
    newFileName = newFileName.ReplaceAll(".MxAOD","");
    newFileName = newFileName.ReplaceAll(Form(".%s",tag.Data()),"");
    ofstream outputSampleList(Form("%s/%s",outputDir.Data(),newFileName.Data()));
    outputSampleList << "root://eosatlas/" << inputDir << "/" << currSample 
		     << std::endl;
    outputSampleList.close();
  }
  sampleFile.close();
}
