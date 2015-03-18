////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMSigParam.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/03/2015                                                          //
//                                                                            //
//  This class parameterizes the resonance shape of the SM Higgs including    //
//  the SM and DM production modes.                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMSigParam.h"

/**
   Initialize the SigParam class.
*/
DMSigParam::DMSigParam(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions) {
  std::cout << std::endl << "DMSigParam::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobname;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMSigParam/",masterOutput.Data(),jobName.Data());
  System(Form("mkdir -vp %s",outputDir.Data()));
    
  // Either load the masspoints from file or create new ones:
  if (options.Contains("FromFile")) loadSigParamFromFile();
  else createNewSigParam();
  return;
}

/**
   Get the name of the output textfile for the given category index.
*/
TString DMSigParam::getTextFileName(int cateIndex, TString process) {
  //TString name = Form("%s/%s_%d.txt",outputDir.Data(),cateScheme.Data(),i_f);
  return name;
}
   
/**
   Create new masspoints by looping over the TTree.
*/
void DMSigParam::createNewSigParam() {
  std::cout << "DMSigParam: creating new signal fit from tree." << std::endl;
}
