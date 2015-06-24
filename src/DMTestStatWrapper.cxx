////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMTestStatWrapper.cxx                                               //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 19/04/2015                                                          //
//                                                                            //
//  Includes a main method for using the DMTestStat.cxx class. This is useful //
//  for grid jobs.                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMTestStat.h"

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
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/workspaces/rootfiles/workspaceDM_%s.root",
			    masterOutput.Data(), jobName.Data(), 
			    DMSignal.Data());
  TString copiedFile = Form("workspaceDM_%s.root",DMSignal.Data());
  system(Form("cp %s %s",originFile.Data(),copiedFile.Data()));
    
  // Load the RooWorkspace and ModelConfig:
  TFile inputFile(copiedFile, "read");
  RooWorkspace *workspace = (RooWorkspace*)inputFile.Get("combinedWS");
  
  DMTestStat *ts = new DMTestStat(jobName, DMSignal, cateScheme, "new",
				  workspace);
  ts->calculateNewCL();
  ts->calculateNewP0();
  if (ts->fitsAllConverged()) {
    std::cout << "DMTestStatWrapper: All OK!" << std::endl;
  }
  
  inputFile.Close();
  system(Form("rm %s",copiedFile.Data()));
  return 0;
}
