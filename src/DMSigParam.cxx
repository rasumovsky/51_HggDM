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
    
  outputDir = Form("%s/%s/SigParam",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/Plots",outputDir.Data())); 
  system(Form("mkdir -vp %s/all",outputDir.Data()));
  for (int i_p = 0; i_p < nProdModes; i_p++) {
    system(Form("mkdir -vp %s/%s",outputDir.Data(),(sigProdModes[i_p]).Data()));
  }
  
  // Either load the signal parameterization from file or create new ones:
  if (options.Contains("FromFile")) loadSigParamFromFile();
  else createNewSigParam();
  return;
}

/**
   Get the name of the output textfile for the given category index. fileType
   can either be "fit" or "yield".
*/
TString DMSigParam::getSigParamFileName(int cateIndex, TString production,
					TString fileType) {
  TString name = Form("%s/%s/%s_%s_%d.txt",outputDir.Data(),production.Data(),
		      fileType.Data(),cateScheme.Data(),cateIndex);
  return name;
}

/**
   Create new masspoints by looping over the TTree.
*/
void DMSigParam::createNewSigParam() {
  std::cout << "DMSigParam: creating new signal fit from tree." << std::endl;
}
