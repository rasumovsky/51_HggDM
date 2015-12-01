////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  MakeSysList.C                                                             //
//                                                                            //
//  Created: Andrew Hard 01/12/2015                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

// ROOT libraries:
#include "TBranch.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

/**
   -----------------------------------------------------------------------------
   Main method looks at an MxAOD with systematic variations and compiles a list
   of those variations by finding all the cutFlow branches.
   @param inputFile - The input ROOT file containing the MxAOD with systematics.
*/
void MakeSysList(TString inputFileName) {
  
  // Output list:
  ofstream outputList("newSystematicsList.txt");
  
  // Load 
  TFile inputFile(inputFileName);
  TTree *tree = (TTree*)inputFile.Get("CollectionTree");
  TObjArray* branchArray = tree->GetListOfBranches();
  for (int i_b = 0; i_b < branchArray->GetEntries(); i_b++) {
    TBranch *currBranch = (TBranch*)branchArray->At(i_b);
    TString currName = currBranch->GetName();
    if (currName.Contains("cutFlow")) {
      TString printName = currName;
      printName.ReplaceAll("HGamEventInfo_", "");
      printName.ReplaceAll("HGamEventInfo.", "");
      printName.ReplaceAll("AuxDyn_cutFlow", "");
      printName.ReplaceAll("AuxDyn.cutFlow", "");
      outputList << printName << std::endl;
    }
    //delete currBranch;
  }
  
  outputList.close();
  //delete tree;
  //delete branchArray;
  //inputFile.Close();
}
