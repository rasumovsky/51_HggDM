////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This program is useful as an interface to the H->diphoton + DM analysis   //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - MassPoints                                                            //
//    - SigParam                                                              //
//    - BkgModel                                                              //
//    - Workspace                                                             //
//    - ToyMC                                                                 //
//    - CalcCLs                                                               //
//                                                                            //
//  Need to rethink the DMSigParam handling of the RooDataSet. Maybe we       //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

//   Compiles the executable required by the job options.

void MakeExe(TString exename) {
  
  // recompile all executables before running...
  system(Form("rm bin/%s",exename.Data()));
  system(Form("make bin/%s",exename.Data()));
}

void SubmitWSViaBsub(TString exeJobName, TString exeOption, TString exeSignal) {
  
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s",exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s", exeWorkspace.Data()));
  system(Form("chmod +x %s", jobScriptWorkspace.Data()));
  
  //change permissions for local input files:!!!!!!!!!!!!!
  //system(Form("chmod +x %s/%s/workspace_files/combinedWS.root",master_output.Data(),executable_jobname.Data()));
   
  system(Form("cp -f %s %s/jobFileWorkspace.sh", jobScriptWorkspace.Data(),
	      exe.Data()));
  system(Form("mv Cocoon.tar %s", exe.Data()));
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), exeJobName.Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(), 
			     exeSignal.Data());
  
  // Define the arguments for the job script:
  TString nameJobScript = Form("%s %s %s %s %s %s", jobScriptWorkspace.Data(),
			       exeJobName.Data(), inputFile.Data(),
			       exeOption.Data(), exeWorkspace.Data(),
			       exeSignal.Data());
  
  // Submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   This is the main DMMaster method:
*/
int main( int argc, char **argv ) {
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <MasterJobName> <MasterOption>\n\n",argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterJobName = argv[1];
  TString masterOption = argv[2];
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = false;
  
  // Options for each step:
  TString massPointOptions = "FromFile";
  TString sigParOptions = "FromFile";
  TString bkgModelOptions = "FromFile";
  TString workspaceOptions = "FromFile_nosys";
  TString testStatOptions = "FromFile";
  
  //--------------------------------------//
  // Step 1: Make or load mass points:
  if (masterOption.Contains("MassPoints")) {
    cout << "DMMaster: Step 1 - Make mass points." << endl;
    
    // Loop over SM, DM, MC samples:
    for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, sigSMModes[i_SM],
					  cateScheme, massPointOptions);
    }
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, sigDMModes[i_DM],
					  cateScheme, massPointOptions);
    }
    for (int i_MC = 0; i_MC < nMCProcesses; i_MC++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, MCProcesses[i_MC],
					  cateScheme, massPointOptions);
    }
  }
  
  //--------------------------------------//
  // Step 2: Make or load the signal parameterization:
  if (masterOption.Contains("SigParam")) {
    cout << "DMMaster: Step 2 - Make signal parameterization." << endl;
    DMSigparam *sp = new DMSigParam(masterJobName, cateScheme, sigParOptions);
  }
  
  //--------------------------------------//
  // Step 4: Create the background model (spurious signal calculation):
  // REPLACE WITH SPURIOUS SIGNAL CODE.
  /*
  if (masterOption.Contains("BkgModel")) {
    cout << "DMMaster: Step 4 - Making the background model." << endl;
    DMBkgModel *dmb = new DMBkgModel(masterJobName, cateScheme,
				     bkgModelOptions);
  }
  
  */
  
  //--------------------------------------//
  // Step 5.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace")) {
    std::cout << "DMMaster: Step 5.1 - Making the workspaces." << std::endl;
    
    int jobCounterWS = 0;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigmDMModes[i_DM];
      
      if (runInParallel) {
	submitWSViaBsub(masterJobName, workspaceOptions, currDMSignal);
	jobCounterWS++;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(masterJobName, currDMSignal,
					   cateScheme, workspaceOptions);
	if (dmw->fitsAllConverged()) {
	  jobCounterWS++;
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterWS << " jobs" << std::endl;
  }
  /*
  //--------------------------------------//
  // Step 5.2: Resubmit any failed workspace jobs:
  if (masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DMMaster: Step 5.2 - Resubmit failed workspace." << std::endl;
    
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(masterJobName);
    vector<TString> resubmitDMSignals = dmc->getResubmitList("DMWorkspace");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitDMSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigmDMModes[i_DM];
      
      if (runInParallel) {
	submitWSViaBsub(exeWorkspace, masterJobName, workspaceOptions, 
			currDMSignal);
	jobCounterWS++;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(masterJobName, currDMSignal,
					   cateScheme, workspaceOptions);
	if (dmw->fitsAllConverged()) {
	  jobCounterWS++;
	}
      }
    }
  }
  

  //--------------------------------------//
  // Step 6: Create pseudoexperiment ensemble (use NPP example):
  
  if (masterOption.Contains("ToyMC")) {
    cout << "DMMaster: Step 6 - Creating pseudoexperiments." << endl;
  }
  */

  //--------------------------------------//
  // Step 7: Calculate the test statistics:
  if (masterOption.Contains("TestStat")) {
    std::cout << "DMMaster: Step 7 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigmDMModes[i_DM];
      
      if (runInParallel) {
	submitTSViaBsub(exeTestStat, masterJobName, testStatOptions, 
			currDMSignal);
	jobCounterTS++;
      }
      else {
	DMTestStat *dmts = new DMTestStat(masterJobName, currDMSignal, 
					  testStatOptions, cateScheme);
	if (dmts->fitsAllConverged()) {
	  jobCounterTS++;
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterTS << " jobs" << std::endl;
  }
  
  return 0;
}
