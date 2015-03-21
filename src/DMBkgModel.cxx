////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMBkgModel.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/03/2015                                                          //
//                                                                            //
//  This class implements a broad range of fit functions to use for the       //
//  background model in the H->yy analysis.                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMBkgModel.h"

/**
   Initialize the MassPoint class.
*/
DMBkgModel::DMBkgModel(TString newJobName, TString newSampleName, 
			   TString newCateScheme, TString newOptions) {
  std::cout << std::endl << "DMBkgModel::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobname;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMBkgModel",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
      
  // Either load the masspoints from file or create new ones:
  if (options.Contains("FromFile")) loadMassPointsFromFile();
  else createNewMassPoints();
  return;
}

/**
   Get the name of the output textfile for the given category index.
*/
TString DMBkgModel::getMassPointsFileName(int cateIndex) {
  TString name = Form("%s/%s_%d.txt",outputDir.Data(),cateScheme.Data(),
		      cateIndex);
  return name;
}
