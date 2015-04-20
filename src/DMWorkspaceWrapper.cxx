////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMWorkspaceWrapper.cxx                                              //
//                                                                            //
//  Creator: Andrew Hard,                                                     //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/04/2015                                                          //
//                                                                            //
//  This class builds the workspace for the dark matter analysis fits.        //
//                                                                            //
//  First: build signal and background models.                                //
//  Second: add asimov data function.                                         //
//  Third: make plots a la spin analysis or better yet NPP.                   //
//                                                                            //
//  Note: 18/4/2015. Need to specify a single DM signal process to be used in //
//  the construction of this statistical model. One workspace per DM model.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMWorkspace.h"

int main(int argc, char **argv) {

  // Check that arguments are provided.
  if (argc < 5) { 
    std::cout << "\nUsage: " << argv[0] 
	      << " <jobName> <DMSignal> <cateScheme> <options>" << std::endl;
    exit(0);
  }
  
  TString jobName = argv[1];
  TString DMSignal = argv[2];
  TString cateScheme = argv[3];
  TString options = argv[4];
  
  DMWorkspace *dmw = new DMWorkspace(jobName, DMSignal, cateScheme, options);
  if (dmw->fitsAllConverged()) {
    std::cout << "DMWorkspaceWrapper: All OK!" << std::endl;
  }
  return 0;
}
