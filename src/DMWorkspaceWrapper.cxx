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
////////////////////////////////////////////////////////////////////////////////

#include "DMWorkspace.h"

int main(int argc, char **argv) {

  // Check that arguments are provided.
  if (argc < 4) { 
    std::cout << "\nUsage: " << argv[0] 
	      << " <configFile> <DMSignal> <options>" << std::endl;
    exit(0);
  }
  TString configFile = argv[1];
  TString DMSignal = argv[2];
  TString options = argv[3];
  
  // Load the workspace tool:
  DMWorkspace *dmw = new DMWorkspace(configFile, DMSignal, options);
  if (dmw->fitsAllConverged()) {
    std::cout << "DMWorkspaceWrapper: All OK!" << std::endl;
  }
  delete dmw;
  return 0;
}
