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
//    - ResubmitWorkspace                                                     //
//    - TestStat                                                              //
//    - ResubmitTestStat                                                      //
//    - MuLimit                                                               //
//                                                                            //
//  Need to rethink the DMSigParam handling of the RooDataSet. Maybe we       //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

using namespace std;
using namespace DMAnalysis;

/**
   Compiles the executable required by the job options.
   @param exeName - the name of the executable to compile.
*/
void makeExe(TString exeName) {
  // recompile all executables before running...
  system(Form("rm %s/bin/%s", packageLocation.Data(), exeName.Data()));
  system(Form("make %s/bin/%s", packageLocation.Data(), exeName.Data()));
}

/**
   Submits the workspace jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeDMSignal - the signal to process in the executable.
   @param exeCateScheme - the categorization scheme for the executable.
*/
void submitWSViaBsub(TString exeJobName, TString exeOption, TString exeDMSignal,
		     TString exeCateScheme) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_DMWorkspace", clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s", exeWorkspace.Data()));
  system(Form("chmod +x %s/%s", packageLocation.Data(), 
	      jobScriptWorkspace.Data()));
     
  system(Form("cp -f %s/%s %s/jobFileWorkspace.sh", packageLocation.Data(), 
	      jobScriptWorkspace.Data(), exe.Data()));
  system(Form("mv Cocoon.tar %s", exe.Data()));
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), exeJobName.Data(),
			     exeDMSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(), 
			     exeDMSignal.Data());
  
  // Define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileWorkspace.sh %s %s %s %s %s %s",
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       exeOption.Data(), exeWorkspace.Data(),
			       exeDMSignal.Data(), exeCateScheme.Data());
  // Submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   Submits the test statistics jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeDMSignal - the signal to process in the executable.
*/
void submitTSViaBsub(TString exeJobName, TString exeOption, TString exeDMSignal,
		     TString exeCateScheme) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_DMTestStat", clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s", exeTestStat.Data()));
  system(Form("chmod +x %s/%s", packageLocation.Data(), 
	      jobScriptTestStat.Data()));
  system(Form("chmod +x %s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
	      masterOutput.Data(), exeJobName.Data(), exeDMSignal.Data()));
  system(Form("cp -f %s/%s %s/jobFileTestStat.sh", packageLocation.Data(), 
	      jobScriptTestStat.Data(), exe.Data()));
  system(Form("mv Cocoon.tar %s", exe.Data()));
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(),
			     exeJobName.Data(), exeDMSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(),
			     exeDMSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileTestStat.sh %s %s %s %s %s %s", 
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       exeTestStat.Data(), exeDMSignal.Data(),
			       exeCateScheme.Data(), exeOption.Data());
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(),
	      nameErrFile.Data(), nameJobScript.Data()));
}


/**
   Submits the mu limit jobs to the lxbatch server. 
   @param exeJobName - the job name.
   @param exeOption - the job options for the executable.
   @param exeDMSignal - the signal to process in the executable.
*/
void SubmitMuLimitViaBsub(TString exeJobName, TString exeOption,
			  TString exeDMSignal) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_DMMuLimit", clusterFileLocation.Data(),
		     exeJobName.Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s", exeMuLimit.Data()));
  system(Form("chmod +x %s", jobScriptMuLimit.Data()));
  system(Form("chmod +x %s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root", 
	      masterOutput.Data(), exeJobName.Data(), exeDMSignal.Data()));
  system(Form("cp -f %s/%s %s/jobFileMuLimit.sh", packageLocation.Data(), 
	      jobScriptMuLimit.Data(), exe.Data()));
  system(Form("mv Cocoon.tar %s", exe.Data()));
 
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), exeJobName.Data(),
			     exeDMSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), exeJobName.Data(),
			     exeDMSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJobScript = Form("%s/jobFileMuLimit.sh %s %s %s %s %s", 
			       exe.Data(), exeJobName.Data(), inputFile.Data(),
			       exeMuLimit.Data(), exeDMSignal.Data(),
			       exeOption.Data());
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJobScript.Data()));
}

/**
   This is the main DMMaster method:
*/
int main (int argc, char **argv) {
  // Check arguments:
  if (argc < 4) {
    printf("\nUsage: %s <MasterJobName> <MasterOption>\n\n",argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterJobName = argv[1];
  TString masterOption = argv[2];
  TString masterCateScheme = argv[3];
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = false;
  
  // Options for each step:
  TString massPointOptions = "New";//"FromFile";
  TString sigParamOptions  = "New";//"FromFile";
  TString bkgModelOptions  = "New";//"FromFile";
  TString workspaceOptions = "New_nosys";//"FromFile_nosys";
  TString testStatOptions  = "New";//"FromFile";
  TString muLimitOptions   = "null";

  //--------------------------------------//
  // Compile any wrappers that need to run remotely:
  if (runInParallel) {
    if (masterOption.Contains("Workspace")) makeExe(exeWorkspace);
    if (masterOption.Contains("TestStat")) makeExe(exeTestStat);
    if (masterOption.Contains("MuLimit")) makeExe(exeMuLimit);
  }
  
  //--------------------------------------//
  // Step 1: Make or load mass points:
  if (masterOption.Contains("MassPoints")) {
    cout << "DMMaster: Step 1 - Make mass points." << endl;
    
    // Loop over SM, DM, MC samples:
    for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, sigSMModes[i_SM],
					  masterCateScheme, massPointOptions,
					  NULL);
    }
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, sigDMModes[i_DM],
					  masterCateScheme, massPointOptions,
					  NULL);
    }
    for (int i_MC = 0; i_MC < nMCProcesses; i_MC++) {
      DMMassPoints *mp = new DMMassPoints(masterJobName, MCProcesses[i_MC],
					  masterCateScheme, massPointOptions,
					  NULL);
    }
  }
  
  //--------------------------------------//
  // Step 2: Make or load the signal parameterization:
  if (masterOption.Contains("SigParam")) {
    cout << "DMMaster: Step 2 - Make signal parameterization." << endl;
    DMSigParam *sp = new DMSigParam(masterJobName, masterCateScheme,
				    sigParamOptions, NULL);
  }
  
  //--------------------------------------//
  // Step 4: Create the background model (spurious signal calculation):
  // REPLACE WITH SPURIOUS SIGNAL CODE.
  /*
  if (masterOption.Contains("BkgModel")) {
    cout << "DMMaster: Step 4 - Making the background model." << endl;
    DMBkgModel *dmb = new DMBkgModel(masterJobName, masterCateScheme,
				     bkgModelOptions);
  }
  
  */
  
  //--------------------------------------//
  // Step 5.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace") && 
      !masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DMMaster: Step 5.1 - Making the workspaces." << std::endl;
    
    int jobCounterWS = 0;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigDMModes[i_DM];
      
      // Temporary, for testing workspaces:
      if (i_DM > 0) {
	std::cout << "DMMaster: Exiting prematurely for test reasons." 
		  << std::endl;
	exit(0);
      }
      
      if (runInParallel) {
	submitWSViaBsub(masterJobName, workspaceOptions, currDMSignal,
			masterCateScheme);
	jobCounterWS++;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(masterJobName, currDMSignal,
					   masterCateScheme, workspaceOptions);
	if (dmw->fitsAllConverged()) {
	  jobCounterWS++;
	}
	else {
	  std::cout << "DMMaster: Fits in previous workspace job failed for "
		    << currDMSignal << " and " << masterCateScheme << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterWS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Resubmit any failed workspace jobs:
  if (masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DMMaster: Step 5.2 - Resubmit failed workspace." << std::endl;
    
    int jobCounterWS = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(masterJobName);
    vector<TString> resubmitDMSignals = dmc->getResubmitList("DMWorkspace");
    dmc->printResubmitList("DMWorkspace");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitDMSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitDMSignals.size(); i_DM++) {
      TString currDMSignal = resubmitDMSignals[i_DM];
      
      if (runInParallel) {
	submitWSViaBsub(exeWorkspace, masterJobName, workspaceOptions, 
			currDMSignal);
	jobCounterWS++;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(masterJobName, currDMSignal,
					   masterCateScheme, workspaceOptions);
	if (dmw->fitsAllConverged()) {
	  jobCounterWS++;
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterWS << " jobs" << std::endl;
  }
  
  /*
  //--------------------------------------//
  // Step 6: Create pseudoexperiment ensemble (use NPP example):
  
  if (masterOption.Contains("ToyMC")) {
    cout << "DMMaster: Step 6 - Creating pseudoexperiments." << endl;
  }
  */

  //--------------------------------------//
  // Step 7.1: Calculate the test statistics:
  if (masterOption.Contains("TestStat") && 
      !masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DMMaster: Step 7.1 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigDMModes[i_DM];
      
      if (runInParallel) {
	submitTSViaBsub(masterJobName, testStatOptions, currDMSignal,
			masterCateScheme);
	jobCounterTS++;
      }
      else {
	DMTestStat *dmts = new DMTestStat(masterJobName, currDMSignal, 
					  masterCateScheme, testStatOptions,
					  NULL);
	if (dmts->fitsAllConverged()) {
	  jobCounterTS++;
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 7.2: Resubmit any failed test statistics jobs:
  if (masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DMMaster: Step 7.2 - Resubmit failed test stat." << std::endl;
    
    int jobCounterTS = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(masterJobName);
    vector<TString> resubmitDMSignals = dmc->getResubmitList("DMTestStat");
    dmc->printResubmitList("DMTestStat");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitDMSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitDMSignals.size(); i_DM++) {
      TString currDMSignal = resubmitDMSignals[i_DM];
      
      if (runInParallel) {
	submitTSViaBsub(masterJobName, testStatOptions, currDMSignal, 
			masterCateScheme);
      	jobCounterTS++;
      }
      else {
	DMTestStat *dmts = new DMTestStat(masterJobName, currDMSignal,
					  masterCateScheme, testStatOptions,
					  NULL);
	if (dmts->fitsAllConverged()) {
	  jobCounterTS++;
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 8.1: Calculate the limits on the dark matter signal strength.
  if (masterOption.Contains("MuLimit") &&
      !masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DMMaster: Step 8.1 - Calculate 95%CL mu value." << std::endl;

    int jobCounterML = 0;
    for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
      TString currDMSignal = sigDMModes[i_DM];
      
      if (runInParallel) {
	submitMLViaBsub(masterJobName, muLimitOptions, currDMSignal);
      }
      else {
	system(Form("./bin/%s %s %s %s", exeMuLimit.Data(),
		    masterJobName.Data(), currDMSignal.Data(),
		    muLimitOptions.Data()));
      }
      jobCounterML++;
    }
    std::cout << "Submitted/completed " << jobCounterML << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 8.2: Resubmit any failed test statistics jobs:
  if (masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DMMaster: Step 8.2 - Resubmit failed mu limits." << std::endl;
    
    int jobCounterML = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(masterJobName);
    vector<TString> resubmitDMSignals = dmc->getResubmitList("DMMuLimit");
    dmc->printResubmitList("DMMuLimit");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitDMSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitDMSignals.size(); i_DM++) {
      TString currDMSignal = resubmitDMSignals[i_DM];
      
      if (runInParallel) {
	submitMLViaBsub(masterJobName, muLimitOptions, currDMSignal);
      }
      else {
	system(Form("./bin/%s %s %s %s", exeMuLimit.Data(),
		    masterJobName.Data(), currDMSignal.Data(),
		    muLimitOptions.Data()));
      }
      jobCounterML++;
    }
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  
  return 0;
}
