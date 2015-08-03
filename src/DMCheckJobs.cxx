////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMCheckJobs.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/04/2015                                                          //
//                                                                            //
//  This class checks to see whether jobs of a particular type have finished. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMCheckJobs.h"

/**
   Initialize the class for checking job statuses.
   @param configFileName - the name of the config file for the analysis.
*/
DMCheckJobs::DMCheckJobs(TString configFileName) {
  m_config = new Config(configFileName);
  updateJobStatus("DMWorkspace");
  updateJobStatus("DMTestStat");
  updateJobStatus("DMMuLimit");
  return;
}

/**
   Get the number of jobs that need to be resubmitted.
   @param jobType - the type of job (DMWorkspace, DMTestStat, DMMuLimit).
   @returns - the number of jobs that failed the first attempt.
*/
int DMCheckJobs::getNumberToResubmit(TString jobType) {
  std::vector<TString> currList = getResubmitList(jobType);
  return (int)currList.size();
}

/**
   Get the list of signal points that must be resubmitted.
   @param jobType - the type of job (DMWorkspace, DMTestStat).
   @returns - a vector of the signal point names.
*/
std::vector<TString> DMCheckJobs::getResubmitList(TString jobType) {
  std::vector<TString> result; result.clear();
  if (jobType.EqualTo("DMWorkspace")) result = m_listDMWorkspace;
  else if (jobType.EqualTo("DMTestStat")) result = m_listDMTestStat;
  else if (jobType.EqualTo("DMMuLimit")) result = m_listDMMuLimit;
  return result;
}

/**
   Update the status of jobs from a particular program.
   @param jobType - the type of job (DMWorkspace, DMTestStat)
   @returns - void. Updates the lists of failed jobs.
*/
void DMCheckJobs::updateJobStatus(TString jobType) {
  
  // Save names of failed jobs:
  m_listDMWorkspace.clear();
  m_listDMTestStat.clear();
  m_listDMMuLimit.clear();
  
  // Then loop over submissions to see whether output files exist:
  std::vector<TString> sigDMModes = m_config->getStr("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    TString currDMSignal = sigDMModes[i_DM];
    
    TString fileName;
    TString fullName;
    if (jobType.EqualTo("DMWorkspace")) {
      fileName = Form("workspaceDM_%s.root", currDMSignal.Data());
      fullName = Form("%s/%s/DMWorkspace/rootfiles/%s", 
		      (m_config->getStr("masterOutput")).Data(),
		      (m_config->getStr("jobName")).Data(), fileName.Data());
    }
    if (jobType.EqualTo("DMTestStat")) {
      fileName = Form("CL_values_%s.txt", currDMSignal.Data());
      fullName = Form("%s/%s/DMTestStat/CL/%s", 
		      (m_config->getStr("masterOutput")).Data(), 
		      (m_config->getStr("jobName")).Data(), fileName.Data());
    }
    
    if (jobType.EqualTo("DMMuLimit")) {
      fileName = Form("text_CLs_%s.txt", currDMSignal.Data());
      fullName = Form("%s/%s/DMMuLimit/single_files/%s", 
		      (m_config->getStr("masterOutput")).Data(), 
		      (m_config->getStr("jobName")).Data(), fileName.Data());
    }

    // Then test the existence of the file.
    std::ifstream testFile(fullName);
    if (!testFile) {
      if (jobType.EqualTo("DMWorkspace")) {
	m_listDMWorkspace.push_back(currDMSignal);
      }
      if (jobType.EqualTo("DMTestStat")) {
	m_listDMTestStat.push_back(currDMSignal);
      }
      if (jobType.EqualTo("DMMuLimit")) {
	m_listDMMuLimit.push_back(currDMSignal);
      }
    }
  }
}

/**
   Print the list and number of failed jobs.
   @param jobType - the type of job (DMWorkspace, DMTestStat).
   @returns void. Prints to the terminal.
*/
void DMCheckJobs::printResubmitList(TString jobType) {
  // Get the relevant list of failed jobs:
  std::vector<TString> currList = getResubmitList(jobType);
  std::cout << "Failed to make the following " << jobType << " files ("
	    << currList.size() << " total)" << std::endl;
  for (int i_f = 0; i_f < (int)currList.size(); i_f++) {
    std::cout << "\t" << currList[i_f] << std::endl;
  }
}
