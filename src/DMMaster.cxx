////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 03/08/2015                                                          //
//                                                                            //
//  The main() method is an interface to the H->diphoton + DM analysis tools. //
//  It centralizes the commands for creating inputs, plots, workspaces, and   //
//  statistical results. Some of the commands will rely on accessing classes  //
//  (mass points, signal parameterization), while others will use system      //
//  commands to submit jobs to various clusters.                              //
//                                                                            //
//  To run:                                                                   //
//    ./bin/DMMaster <MasterOption> <configFileName>                          //
//                                                                            //
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - Cleanup                                                               //
//    - MassPoints                                                            //
//    - GetSystematics                                                        //
//    - RankSystematics                                                       //
//    - PlotVariables                                                         //
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
//    - OptAnalysis                                                           //
//                                                                            //
//  Need to rethink the DMSigParam handling of the RooDataSet. Maybe we       //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

/**
   -----------------------------------------------------------------------------
   Recursive method to alter selection cuts and submit jobs. The recursive case
   is when the current cut index is less than the total number of cuts to be 
   varied. The base case is reached when the cutIndex = # cuts - 1 (all cuts 
   have been varied). Upon reaching the base case, a job is submitted. 
   @param exeConfigOrigin - The original analysis config file.
   @param exeOption - The executable options for the job.
   @param cutIndex - The index of the current cut to be modified.
   @param cutName - A vector of cut names.
   @param cutVal - A vector of cut values.
*/
void recursiveOptimizer(TString exeConfigOrigin, TString exeOption, 
			int cutIndex, std::vector<TString> cutName,
			std::vector<double> cutVal) {
  // Get the number of cuts:
  int nCuts = m_config->getInt("NOptVar");
  
  // The recursive case:
  if (cutIndex < nCuts) {
    // Get current cut information:
    std::vector<double> cutPositions
      = m_config->getNumV(Form("OptVarPos%d",cutIndex));
    
    // Loop over the current cut positions:
    for (int i_c = 0; i_c < (int)cutPositions.size(); i_c++) {
      
      // Update the value of the cut identified by "cutIndex":
      std::vector<TString> currCutN = cutName;
      std::vector<double> currCutV = cutVal;
      if ((int)cutName.size() == cutIndex) {
	currCutN.push_back(m_config->getStr(Form("OptVar%d",cutIndex)));
	currCutV.push_back(cutPositions[i_c]);
      }
      else {
	currCutN[cutIndex] = m_config->getStr(Form("OptVar%d",cutIndex));
	currCutV[cutIndex] = cutPositions[i_c];
      }
      
      //Then call recursiveOptimizer.
      recursiveOptimizer(exeConfigOrigin, exeOption, cutIndex+1, currCutN,
			 currCutV);
      
      // Only reach one base case (submit 1 job) if in test mode:
      if (m_config->getBool("doTestMode")) break;
    }
  }
  
  // The base case (No more cuts to change, so just submit job:
  else {
    // First copy the config file (special config file) and make changes:
    //  - change cut values.
    //  - change output directory.
    ifstream inputConfig;
    inputConfig.open(exeConfigOrigin);
    TString exeConfigNew = Form("exeConfig%d.cfg", m_jobIndex);
    ofstream outputConfig; outputConfig.open(exeConfigNew);
    std::string key;
    // Loop over each line of the config file:
    while (!inputConfig.eof()) {
      std::getline(inputConfig, key);
      TString currLine = TString(key);
            
      // Check if the current line specifies cut information:
      bool lineSpecifiesCut = false; 
      TString specifiedName = "";
      double specifiedVal = 0.0;
      for (int i_c = 0; i_c < (int)cutName.size(); i_c++) {
	if (currLine.Contains(cutName[i_c]) && !currLine.Contains("#")) {
	  lineSpecifiesCut = true;
	  specifiedName = cutName[i_c];
	  specifiedVal = cutVal[i_c];
	}
      }
      
      // Change the cut value, if the current line specifies cut information:
      if (lineSpecifiesCut) {
	outputConfig << specifiedName << ": \t" << specifiedVal << std::endl;
      }
      // Make the output local (on lxbatch as opposed to afs):
      else if (currLine.Contains("masterOutput:") && !currLine.Contains("#")){
	outputConfig << "masterOutput: \t" << "." << std::endl;
      }
      else if (currLine.Contains("NoSMProdModes:")&&!currLine.Contains("#")) {
	outputConfig << "NoSMProdModes: \t" << "NO" << std::endl;
      }
      // Just copy the configuration otherwise:
      else {
	outputConfig << currLine << std::endl;
      }
    }
    inputConfig.close();
    outputConfig.close();
    
    // BEGIN JOB SUBMISSION SECTION!
    
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
    TString tempDir = Form("KillMe%d", m_jobIndex);
    TString tarFile = Form("Cocoon%d.tar", m_jobIndex);
    system(Form("mkdir -vp %s", tempDir.Data()));
    system(Form("cp %s/bin/%s %s/", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("exeMaster")).Data(), tempDir.Data()));
    system(Form("mv %s %s/", exeConfigNew.Data(), tempDir.Data()));
    system(Form("cp %s/%s %s/",
		(m_config->getStr("packageLocation")).Data(),
		(m_config->getStr("jobScriptMaster")).Data(), 
		tempDir.Data()));
    system(Form("mv %s %s/", tempDir.Data(), exe.Data()));
    
    // Is this necessary? Probably...
    system(Form("cp -f %s/%s %s/jobFileWorkspace.sh", 
		(m_config->getStr("packageLocation")).Data(), 
		(m_config->getStr("jobScriptMaster")).Data(), exe.Data()));
    
    //TString inputFile = Form("%s/%s", exe.Data(), tarFile.Data());
    TString inputFile = Form("%s/%s", exe.Data(), tempDir.Data());
    TString nameOutFile = Form("%s/out/%s_%d.out", dir.Data(), 
			       (m_config->getStr("jobName")).Data(),m_jobIndex);
    TString nameErrFile = Form("%s/err/%s_%d.err", dir.Data(), 
			       (m_config->getStr("jobName")).Data(),m_jobIndex);
    
    // Define the arguments for the job script:
    TString nameJScript = Form("%s/jobFileWorkspace.sh %s %s %s %s %s %d",
			       exe.Data(),
			       (m_config->getStr("jobName")).Data(),
			       exeConfigNew.Data(),
			       inputFile.Data(),
			       exeOption.Data(),
			       (m_config->getStr("exeMaster")).Data(),
			       m_jobIndex);
    // Submit the job:
    system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
		nameErrFile.Data(), nameJScript.Data()));
    
    // END JOB SUBMISSION SECTION!
    
    // Note: job script should copy output files to a new output directory.
    system(Form("rm -rf %s", tempDir.Data()));
    
    m_headFile << m_jobIndex;
    for (int i_c = 0; i_c < (int)cutName.size(); i_c++) {
      if (i_c == (int)cutName.size()-1) {
	m_headFile << " " << cutVal[i_c] << std::endl;
      }
      else {
	m_headFile << " " << cutVal[i_c];
      }
    }
    m_jobIndex++;
  }
}

/**
   -----------------------------------------------------------------------------
   Submit the DMMainMethod to run remotely on lxbatch.
   @param exeConfigOrigin - The original config file for batch jobs.
   @param exeOption - The option for the executable...
*/
void submitToOptimize(TString exeConfigOrigin, TString exeOption) {
  std::cout << "DMMaster: Preparing to run myself remotely for optimization!"
	    << std::endl;
  
  // Count the jobs as they are submitted, for bookkeeping:
  m_jobIndex = 0;
  
  // An output file to track job indices and cuts:
  system(Form("mkdir -vp %s/%s/DMMaster",
	      (m_config->getStr("masterOutput")).Data(), 
	      (m_config->getStr("jobName")).Data()));
  m_headFile.open(Form("%s/%s/DMMaster/jobSummary.txt", 
		       (m_config->getStr("masterOutput")).Data(), 
		       (m_config->getStr("jobName")).Data()));
  
  // The first line of the output file should just list the cut names:
  int nCuts = m_config->getInt("NOptVar");
  m_headFile << "Index \t";
  for (int i_c = 0; i_c < nCuts; i_c++) {
    if (i_c == nCuts - 1) {
      m_headFile << m_config->getStr(Form("OptVar%d",i_c)) << std::endl;
    }
    else {
      m_headFile << m_config->getStr(Form("OptVar%d",i_c)) << " \t";
    }
  }
  // The remaining lines are filled with values in the recursive method below.
  
  // Call the recursive function:
  std::vector<TString> cutName; cutName.clear();
  std::vector<double> cutVal; cutVal.clear();
  recursiveOptimizer(exeConfigOrigin, exeOption, 0, cutName, cutVal);
  
  // Close the file that records job indices and cuts:
  m_headFile.close();
  
  std::cout << "DMMaster: Submitted " << m_jobIndex
	    << " recursive optimization jobs." << std::endl;
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
void submitMuLimitViaBsub(TString exeConfigFile, TString exeOption,
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
    system(Form("chmod +x %s",(m_config->getStr("jobScriptPseudoExp")).Data()));
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
      
  // Load the config class and file:
  std::cout << "DMMaster: Loading the global config file." << std::endl;
  m_config = new Config(configFileName);
  m_config->printDB();
  TString fullConfigPath = Form("%s/%s",
				(m_config->getStr("packageLocation")).Data(),
				configFileName.Data());
  
  // Submit jobs to bsub or grid, etc.:
  bool runInParallel = m_config->getStr("RunInParallel");
  m_isFirstJob = true;
  
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
  // Step 0: Remove any prior outputs
  if (masterOption.Contains("Cleanup")) {
    std::cout << "DMMaster: Step 0 - Cleaning up previous runs." << std::endl;
    system(Form("rm -rf %s/%s", (m_config->getStr("masterOutput")).Data(), 
		(m_config->getStr("jobName")).Data()));
  }
  
  //--------------------------------------//
  // Step 1.1: Make or load mass points:
  if (masterOption.Contains("MassPoints")) {
    std::cout << "DMMaster: Step 1.1 - Make mass points." << std::endl;
    
    std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
    std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
    std::vector<TString> bkgProcesses = m_config->getStrV("BkgProcesses");
    std::vector<TString> allSamples; allSamples.clear();
    allSamples.push_back("Data");
    allSamples.insert(allSamples.end(),sigSMModes.begin(),sigSMModes.end());
    allSamples.insert(allSamples.end(),sigDMModes.begin(),sigDMModes.end());
    allSamples.insert(allSamples.end(),bkgProcesses.begin(),bkgProcesses.end());
    
    // Loop over all samples:
    for (int i_s = 0; i_s < (int)allSamples.size(); i_s++) {
      // Option to submit remote jobs:
      if (runInParallel) {
      }
      // Otherwise run locally:
      else {
	DMMassPoints *mp = new DMMassPoints(configFileName, allSamples[i_s],
					    massPointOptions, NULL);
	delete mp;
      }
    }
  }
  
  //--------------------------------------//
  // Step 1.2: Make cutflows with experimental systematics:
  if (masterOption.Contains("GetSystematics")) {
    std::cout << "DMMaster: Step 1.2 - Make cutflows w/ systematic variations."
	      << std::endl;
    // Load MxAODs with systematic uncertainties:
    std::vector<TString> sysSamples = m_config->getStrV("SystematicsSamples");
    for (int i_s = 0; i_s < (int)sysSamples.size(); i_s++) {
      // Mass point options will include "Syst":
      TString sysOptions = massPointOptions + "_Syst";
      DMMassPoints *mp 
	= new DMMassPoints(configFileName, sysSamples[i_s], sysOptions, NULL);
      delete mp;
    }
  }

  //--------------------------------------//
  // Step 1.3: Rank the systematic uncertainties for each sample:
  if (masterOption.Contains("RankSystematics")) {
    std::cout << "DMMaster: Step 1.3 - Rank systematic variations by impact."
	      << std::endl;
    SystematicsTool *sh = new SystematicsTool(configFileName);
    
    std::vector<TString> sysSamples = m_config->getStrV("SystematicsSamples");
    for (int i_s = 0; i_s < (int)sysSamples.size(); i_s++) {
      sh->loadAllSys(sysSamples[i_s]);
      sh->saveRankedNormSys(sysSamples[i_s]);
      sh->saveRankedMigrSys(sysSamples[i_s]);
    }
    delete sh;
  }
  
  //--------------------------------------//
  // Step 1.4: Make variable plots:
  if (masterOption.Contains("PlotVariables")) {
    std::cout << "DMMaster: Step 1.4 - Make kinematic variable plots."
	      << std::endl;
    std::vector<TString> plotVariables = m_config->getStrV("PlotVariables");
    for (int i_v = 0; i_v < (int)plotVariables.size(); i_v++) {
      system(Form("./bin/PlotVariables %s %s %s", configFileName.Data(), 
		  plotVariables[i_v].Data(), 
		  (m_config->getStr("PlotVariableOptions")).Data()));
    }
  }
  
  //--------------------------------------//
  // Step 2: Make or load the signal parameterization:
  if (masterOption.Contains("SigParam")) {
    std::cout << "DMMaster: Step 2 - Make signal parameterization." 
	      << std::endl;
    SigParamInterface *spi = new SigParamInterface(configFileName,
						   sigParamOptions);
    delete spi;
  }
  
  //--------------------------------------//
  // Step 3: Create the background model (spurious signal calculation):
  // REPLACE WITH SPURIOUS SIGNAL CODE.
  if (masterOption.Contains("BkgModel")) {
    std::cout << "DMMaster: Step 4 - Making the background model." << std::endl;
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
	delete dmw;
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
	delete dmw;
      }
    }
    delete dmc;
    std::cout << "Resubmitted " << jobCounterWS << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.1: Create pseudoexperiment ensemble:
  TString currToySignal = m_config->getStr("exampleSignal");
  if (masterOption.Contains("TossPseudoExp")) {
    std::cout << "DMMaster: Step 5.1 - Creating pseudoexperiments for signal "
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
    delete dmta;
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
	delete dmts;
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
	delete dmts;
      }
    }
    delete dmc;
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
	//submitMuLimitViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	/*
	TString muCommand = Form(".%s/bin/%s %s %s %s", 
				 (m_config->getStr("packageLocation")).Data(), 
				 (m_config->getStr("exeMuLimit")).Data(),
				 fullConfigPath.Data(), currSignal.Data(), 
				 muLimitOptions.Data());
	*/
	TString muCommand = Form("./bin/%s %s %s %s", 
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
	submitMuLimitViaBsub(fullConfigPath, muLimitOptions, currSignal);
	m_isFirstJob = false;
      }
      else {
	/*
	TString muCommand = Form(".%s/bin/%s %s %s %s", 
				 (m_config->getStr("packageLocation")).Data(), 
				 (m_config->getStr("exeMuLimit")).Data(),
				 fullConfigPath.Data(), currSignal.Data(), 
				 muLimitOptions.Data());
	*/
	TString muCommand = Form("./bin/%s %s %s %s", 
				 (m_config->getStr("exeMuLimit")).Data(),
				 fullConfigPath.Data(), currSignal.Data(), 
				 muLimitOptions.Data());
	std::cout << "Executing following system command: \n\t"
		  << muCommand << std::endl;
	system(muCommand);
      }
      jobCounterML++;
    }
    delete dmc;
    std::cout << "Resubmitted " << jobCounterML << " jobs" << std::endl;
  }
  
  //--------------------------------------//
  // Step 8: Optimize the analysis!
  if (masterOption.Contains("Optimizer")) {
    submitToOptimize(configFileName, m_config->getStr("masterJobOptions"));
  }
  
  //--------------------------------------//
  // Step 9: Plot the results of the optimization
  if (masterOption.Contains("OptAnalysis")) {
    DMOptAnalysis *dmoa = new DMOptAnalysis(configFileName);
    delete dmoa;
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
