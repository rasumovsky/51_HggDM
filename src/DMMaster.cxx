////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 03/08/2015                                                          //
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
//    - TossPseudoExp                                                         //
//    - PlotPseudoExp                                                         //
//    - TestStat                                                              //
//    - ResubmitTestStat                                                      //
//    - MuLimit                                                               //
//    - Optimizer                                                             //
//                                                                            //
//  Need to rethink the DMSigParam handling of the RooDataSet. Maybe we       //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

/**
   -----------------------------------------------------------------------------
   Submit the DMMainMethod to run remotely on lxbatch.
   @param exeConfigOrigin - The original config file for batch jobs.
   @param exeOption - The option for the executable...
*/
void submitToOptimize(TString exeConfigOrigin, TString exeOption) {
  std::cout << "DMMaster: Preparing to run myself remotely for optimization!"
	    << std::endl;
  
  int jobIndex = 0;
  
  // An output file to track job indices and cuts:
  ofstream headFile;
  headFile.open(Form("%s/%s/DMMaster/jobSummary.txt", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data()));
  
  // Loop over first cut value:
  double increment = ((TMath::Pi()/2.0) / m_config->getNum("nPointsPerVar"));
  double theta2 = increment;
  for (double theta2 = 2*increment; theta2 <= (TMath::Pi()/2.0); 
       theta2 += increment) {
    // Loop over second cut value:
    for (double theta1 = increment; theta1 < theta2; theta1 += increment) {
      
      // Convert thetas to ratio:
      double currCut2 = TMath::ATan(theta2);
      double currCut1 = TMath::ATan(theta1);
      
      // First copy the config file (special config file) and make changes:
      //  - change cut values.
      //  - change output directory.
      ifstream inputConfig;
      inputConfig.open(exeConfigOrigin);
      TString exeConfigNew = Form("exeConfig%d.cfg", jobIndex);
      ofstream outputConfig;
      outputConfig.open(exeConfigNew);
      //string key;// need to pipe whole line into key
      char key[256];// need to pipe whole line into key
      while (!inputConfig.eof()) {
	inputConfig.getline(key,256);
	TString currLine = TString(key);
	// Change the category 0/1 division:
	if (currLine.Contains("RatioCut1:") && !currLine.Contains("#")) {
	  outputConfig << "RatioCut1: \t" << currCut1 << std::endl;
	}
	// Change the category 1/2 division:
	else if (currLine.Contains("RatioCut2:") && !currLine.Contains("#")) {
	  outputConfig << "RatioCut2: \t" << currCut2 << std::endl;
	}
	// Make the output local (on lxbatch as opposed to afs):
	else if (currLine.Contains("masterOutput:") && !currLine.Contains("#")){
	  outputConfig << "masterOutput: \t" << "." << std::endl;
	}
	else if (currLine.Contains("NoSMProdModes:")&&!currLine.Contains("#")) {
	  outputConfig << "NoSMProdModes: \t" << "NO" << std::endl;
	}
	// Tell DMMassPoints to copy MxAODs before running off of afs or eos:
	//else if (currLine.Contains("massPointOptions:")) {
	//outputConfig << "massPointOptions: \t" << "NewCopyFiles" << std::endl;
	//}
	// Just copy the configuration otherwise:
	else {
	  outputConfig << currLine << std::endl;
	}
      }
      inputConfig.close();
      outputConfig.close();
            
      ////////////////////////////////////////
      // Begin job submission section:

      // Make directories for job info:
      TString dir = Form("%s/%s_DMMaster",
			 (m_config->getStr("clusterFileLocation")).Data(),
			 (m_config->getStr("jobName")).Data());
      TString out = Form("%s/out", dir.Data());
      TString err = Form("%s/err", dir.Data());
      TString exe = Form("%s/exe", dir.Data());
      system(Form("mkdir -vp %s", out.Data()));
      system(Form("mkdir -vp %s", err.Data()));
      system(Form("mkdir -vp %s", exe.Data()));
      
      // Create .tar file with everything needed to run remotely:
      TString tempDir = Form("KillMe%d", jobIndex);
      TString tarFile = Form("Cocoon%d.tar", jobIndex);
      system(Form("mkdir -vp %s", tempDir.Data()));
      system(Form("cp %s/bin/%s %s/", 
		  (m_config->getStr("packageLocation")).Data(), 
		  (m_config->getStr("exeMaster")).Data(), tempDir.Data()));
      system(Form("mv %s %s/", exeConfigNew.Data(), tempDir.Data()));
      system(Form("cp %s/%s %s/",
		  (m_config->getStr("packageLocation")).Data(),
		  (m_config->getStr("jobScriptMaster")).Data(), 
		  tempDir.Data()));
      //system(Form("tar zcf %s %s/*", tarFile.Data(), tempDir.Data()));
      //system(Form("mv %s %s/", tarFile.Data(), exe.Data()));
      system(Form("mv %s %s/", tempDir.Data(), exe.Data()));
      
      // Is this necessary? Probably...
      system(Form("cp -f %s/%s %s/jobFileWorkspace.sh", 
		  (m_config->getStr("packageLocation")).Data(), 
		  (m_config->getStr("jobScriptMaster")).Data(), exe.Data()));
      
      //TString inputFile = Form("%s/%s", exe.Data(), tarFile.Data());
      TString inputFile = Form("%s/%s", exe.Data(), tempDir.Data());
      TString nameOutFile = Form("%s/out/%s_%d.out", dir.Data(), 
				 (m_config->getStr("jobName")).Data(),jobIndex);
      TString nameErrFile = Form("%s/err/%s_%d.err", dir.Data(), 
				 (m_config->getStr("jobName")).Data(),jobIndex);
      
      // Define the arguments for the job script:
      TString nameJScript = Form("%s/jobFileWorkspace.sh %s %s %s %s %s %d",
				 exe.Data(),
				 (m_config->getStr("jobName")).Data(),
				 exeConfigNew.Data(),
				 inputFile.Data(),
				 exeOption.Data(),
				 (m_config->getStr("exeMaster")).Data(),
				 jobIndex);
      // Submit the job:
      system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
		  nameErrFile.Data(), nameJScript.Data()));
      
      // End job submission section!
      ////////////////////////////////////////
      
      // Note: job script should copy output files to a new output directory.
      system(Form("rm -rf %s", tempDir.Data()));
      
      headFile << jobIndex << " " << currCut1 << " " << currCut2 << std::endl;
      jobIndex++;
      
      if (m_config->getBool("doTestMode")) break;
    }// End of theta1 loop
    if (m_config->getBool("doTestMode")) break;
  }// End of theta2 loop

  // Close the file that records job indices and cuts:
  headFile.close();
  
  std::cout << "DMMaster: Submitted " << jobIndex << " self-optimizing jobs."
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Submits the workspace jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
 */
void submitWSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal) {
  // Make directories for job info:
  TString dir = Form("%s/%s_DMWorkspace",
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeWorkspace")).Data()));
    system(Form("chmod +x %s/%s", (m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptWorkspace")).Data()));
    system(Form("cp -f %s/%s %s/jobFileWorkspace.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptWorkspace")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("jobName")).Data(), 
			     exeSignal.Data());
  
  // Define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileWorkspace.sh %s %s %s %s %s %s",
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     exeOption.Data(),
			     (m_config->getStr("exeWorkspace")).Data(),
			     exeSignal.Data());
  // Submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the test statistics jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
 */
void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal) {  
  // Make directories for job info:
  TString dir = Form("%s/%s_DMTestStat", 
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeTestStat")).Data()));
    system(Form("chmod +x %s/%s", (m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data()));
    system(Form("chmod +x %s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFileTestStat.sh",
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptTestStat")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileTestStat.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeTestStat")).Data(),
			     exeSignal.Data(), exeOption.Data());
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(),
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
*/
void SubmitMuLimitViaBsub(TString exeConfigFile, TString exeOption,
			  TString exeSignal) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_DMMuLimit", 
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeMuLimit")).Data()));
    system(Form("chmod +x %s", (m_config->getStr("jobScriptMuLimit")).Data()));
    system(Form("chmod +x %s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root", 
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFileMuLimit.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptMuLimit")).Data(), exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s.out", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  TString nameErrFile = Form("%s/err/%s_%s.err", dir.Data(), 
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data());
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFileMuLimit.sh %s %s %s %s %s %s", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeMuLimit")).Data(),
			     exeSignal.Data(), exeOption.Data());
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSignal - the signal to process in the executable.
   @param int exeSeed - the seed for the randomized dataset generation.
   @param int exeToysPerJob - the number of toy datasets to create per job.
*/
void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal, int exeSeed, int exeToysPerJob) {
  
  // Make directories for job info:
  TString dir = Form("%s/%s_PseudoExp",
		     (m_config->getStr("clusterFileLocation")).Data(),
		     (m_config->getStr("jobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exePseudoExp")).Data()));
    system(Form("chmod +x %s", (m_config->getStr("jobScriptPseudoExp")).Data()));
    system(Form("chmod +x %s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root", 
		(m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data(), exeSignal.Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptPseudoExp")).Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%s_%d.out", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data(), exeSeed);
  TString nameErrFile = Form("%s/err/%s_%s_%d.err", dir.Data(),
			     (m_config->getStr("jobName")).Data(),
			     exeSignal.Data(), exeSeed);
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %s %d %d", 
			     exe.Data(), (m_config->getStr("jobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exePseudoExp")).Data(),
			     exeSignal.Data(), exeOption.Data(), exeSeed,
			     exeToysPerJob);
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   This is the main DMMaster method:
*/
int main (int argc, char **argv) {
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <option> <configFileName> \n\n", argv[0]);
    exit(0);
  }
  
  // Clock the program:
  clock_t time;
  time = clock();

  // The job name and options (which analysis steps to perform):
  TString masterOption = argv[1];
  TString configFileName = argv[2];
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = false;
  m_isFirstJob = true;
    
  // Load the config class and file:
  std::cout << "DMMaster: Loading the global config file." << std::endl;
  m_config = new Config(configFileName);
  m_config->printDB();
  TString fullConfigPath = Form("%s/%s",
				(m_config->getStr("packageLocation")).Data(),
				configFileName.Data());
  
  // Options for each analysis step:
  TString massPointOptions = m_config->getStr("massPointOptions");
  TString sigParamOptions  = m_config->getStr("sigParamOptions");
  TString bkgModelOptions  = m_config->getStr("bkgModelOptions");
  TString workspaceOptions = m_config->getStr("workspaceOptions");
  TString pseudoExpOptions = m_config->getStr("pseudoExpOptions");
  TString toyPlotOptions   = m_config->getStr("toyPlotOptions");
  TString testStatOptions  = m_config->getStr("testStatOptions");
  TString muLimitOptions   = m_config->getStr("muLimitOptions");
  
  //--------------------------------------//
  // Step 1: Make or load mass points:
  if (masterOption.Contains("MassPoints")) {
    cout << "DMMaster: Step 1 - Make mass points." << endl;
    
    // Loop over SM, DM, MC samples:
    std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
    for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
      DMMassPoints *mp = new DMMassPoints(configFileName, sigSMModes[i_SM],
					  massPointOptions, NULL);
    }
    std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
      DMMassPoints *mp = new DMMassPoints(configFileName, sigDMModes[i_DM],
					  massPointOptions, NULL);
    }
    std::vector<TString> MCProcesses = m_config->getStrV("MCProcesses");
    for (int i_MC = 0; i_MC < (int)MCProcesses.size(); i_MC++) {
      DMMassPoints *mp = new DMMassPoints(configFileName, MCProcesses[i_MC],
					  massPointOptions, NULL);
    }
  }
  
  //--------------------------------------//
  // Step 2: Make or load the signal parameterization:
  if (masterOption.Contains("SigParam")) {
    cout << "DMMaster: Step 2 - Make signal parameterization." << endl;
    SigParamInterface *spi = new SigParamInterface(configFileName,
						   sigParamOptions);
  }
  
  //--------------------------------------//
  // Step 3: Create the background model (spurious signal calculation):
  // REPLACE WITH SPURIOUS SIGNAL CODE.
  if (masterOption.Contains("BkgModel")) {
    cout << "DMMaster: Step 4 - Making the background model." << endl;
    exit(0);
  }
  
  
  //--------------------------------------//
  // Step 4.1: Create the workspace for fitting:
  if (masterOption.Contains("Workspace") && 
      !masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DMMaster: Step 4.1 - Making the workspaces." << std::endl;
    
    int jobCounterWS = 0;
    std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
      TString currSignal = sigDMModes[i_DM];
      if (runInParallel) {
	submitWSViaBsub(fullConfigPath, workspaceOptions, currSignal);
	jobCounterWS++;
	m_isFirstJob = false;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(configFileName, currSignal,
					   workspaceOptions);
	if (dmw->fitsAllConverged()) jobCounterWS++;
	else {
	  std::cout << "DMMaster: Problem with workspace fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterWS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 4.2: Resubmit any failed workspace jobs:
  if (masterOption.Contains("ResubmitWorkspace")) {
    std::cout << "DMMaster: Step 4.2 - Resubmit failed workspace." << std::endl;
    
    int jobCounterWS = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(configFileName);
    vector<TString> resubmitSignals = dmc->getResubmitList("DMWorkspace");
    dmc->printResubmitList("DMWorkspace");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitSignals.size(); i_DM++) {
      TString currSignal = resubmitSignals[i_DM];
      
      if (runInParallel) {
	submitWSViaBsub(fullConfigPath, workspaceOptions, currSignal);
	jobCounterWS++;
	m_isFirstJob = false;
      }
      else {
	DMWorkspace *dmw = new DMWorkspace(configFileName, currSignal,
					   workspaceOptions);
	if (dmw->fitsAllConverged()) jobCounterWS++;
	else {
	  std::cout << "DMMaster: Problem with workspace fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterWS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.1: Create pseudoexperiment ensemble:
  TString currToySignal = m_config->getStr("exampleSignal");
  if (masterOption.Contains("TossPseudoExp")) {
    cout << "DMMaster: Step 5.1 - Creating pseudoexperiments for signal "
	 << currToySignal << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(fullConfigPath, pseudoExpOptions, currToySignal, i_s,
		      nToysPerJob);
      m_isFirstJob = false;
    }
    std::cout << "DMMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Plot pseudo-experiment ensemble results:
  if (masterOption.Contains("PlotPseudoExp")) {
    std::cout << "DMMaster: Step 5.2 - Plot pseudoexperiment results for "
	      << currToySignal << std::endl;    
    DMToyAnalysis *dmta = new DMToyAnalysis(configFileName, currToySignal);
  }
  
  //--------------------------------------//
  // Step 6.1: Calculate the test statistics:
  if (masterOption.Contains("TestStat") && 
      !masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DMMaster: Step 6.1 - Calculating CL and p0." << std::endl;

    int jobCounterTS = 0;
    std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
      TString currSignal = sigDMModes[i_DM];
      
      if (runInParallel) {
	submitTSViaBsub(fullConfigPath, testStatOptions, currSignal);
	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DMTestStat *dmts = new DMTestStat(configFileName, currSignal, 
					  testStatOptions, NULL);
	dmts->calculateNewCL();
	dmts->calculateNewP0();
	if (dmts->fitsAllConverged()) jobCounterTS++;
	else {
	  std::cout << "DMMaster: Problem with test-stat fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Submitted/completed " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 6.2: Resubmit any failed test statistics jobs:
  if (masterOption.Contains("ResubmitTestStat")) {
    std::cout << "DMMaster: Step 6.2 - Resubmit failed test stat." << std::endl;
    
    int jobCounterTS = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(configFileName);
    vector<TString> resubmitSignals = dmc->getResubmitList("DMTestStat");
    dmc->printResubmitList("DMTestStat");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " test statistic jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitSignals.size(); i_DM++) {
      TString currSignal = resubmitSignals[i_DM];
      
      if (runInParallel) {
	submitTSViaBsub(fullConfigPath, testStatOptions, currSignal);
      	jobCounterTS++;
	m_isFirstJob = false;
      }
      else {
	DMTestStat *dmts = new DMTestStat(configFileName, currSignal,
					  testStatOptions, NULL);
	dmts->calculateNewCL();
	dmts->calculateNewP0();
	if (dmts->fitsAllConverged()) jobCounterTS++;
	else {
	  std::cout << "DMMaster: Problem with test-stat fit!" << std::endl;
	  exit(0);
	}
      }
    }
    std::cout << "Resubmitted " << jobCounterTS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 7.1: Calculate the limits on the dark matter signal strength.
  if (masterOption.Contains("MuLimit") &&
      !masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DMMaster: Step 7.1 - Calculate 95%CL mu value." << std::endl;

    int jobCounterML = 0;
    std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
      TString currSignal = sigDMModes[i_DM];
      
      if (runInParallel) {
	submitMLViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	TString muCommand = Form(".%s/bin/%s %s %s %s", 
				 (m_config->getStr("packageLocation")).Data(), 
				 (m_config->getStr("exeMuLimit")).Data(),
				 fullConfigPath.Data(), currSignal.Data(), 
				 muLimitOptions.Data());
	std::cout << "Executing following system command: \n\t"
		  << muCommand << std::endl;
	system(muCommand);
      }
      jobCounterML++;
    }
    std::cout << "Submitted/completed " << jobCounterML << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 7.2: Resubmit any failed mu limit jobs:
  if (masterOption.Contains("ResubmitMuLimit")) {
    std::cout << "DMMaster: Step 7.2 - Resubmit failed mu limits." << std::endl;
    
    int jobCounterML = 0;
    // Get the points to resubmit:
    DMCheckJobs *dmc = new DMCheckJobs(configFileName);
    vector<TString> resubmitSignals = dmc->getResubmitList("DMMuLimit");
    dmc->printResubmitList("DMMuLimit");
    
    // Then resubmit as necessary:
    std::cout << "Resubmitting " << (int)resubmitSignals.size()
	      << " workspace jobs." << std::endl;
    for (int i_DM = 0; i_DM < (int)resubmitSignals.size(); i_DM++) {
      TString currSignal = resubmitSignals[i_DM];
      
      if (runInParallel) {
	submitMLViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	system(Form(".%s/bin/%s %s %s %s", 
		    (m_config->getStr("packageLocation")).Data(), 
		    (m_config->getStr("exeMuLimit")).Data(),
		    fullConfigPath.Data(), currSignal.Data(),
		    muLimitOptions.Data()));
      }
      jobCounterML++;
    }
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 8: Optimize the analysis!
  if (masterOption.Contains("Optimizer")) {
    submitToOptimize(configFileName, m_config->getStr("masterJobOptions"));
  }
  
  //--------------------------------------//
  // Terminate successfully:
  std::cout << "DMMaster: All requested processes have completed." <<std::endl;
  // Clock the program:
  time = clock() - time;
  printf("\tRunning required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return 0;
}

