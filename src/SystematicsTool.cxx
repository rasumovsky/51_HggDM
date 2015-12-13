////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SystematicsTool.cxx                                                 //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 12/12/2015                                                          //
//                                                                            //
//  WARNING!!!!!! THIS PROGRAM HAS YET TO IMPLEMENT THE CALCULATION OF SYST   //
//  FROM THE LOADED YIELD VALUES. THIS IS ESSENTIAL!!!
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "SystematicsTool.h"

/**
   -----------------------------------------------------------------------------
   Initialize the SystematicsTool class and make a new RooCategory.
   @param newConfigFile - The name of the analysis config file.
*/
SystematicsTool::SystematicsTool(TString newConfigFile) {
  std::cout << "\nSystematicsTool::Initializing..." << std::endl;
  
  // Assign member variables:
  m_configFile = newConfigFile;
  
  // Load the configuration for the analysis:
  m_config = new Config(m_configFile);
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/SystematicsTool", 
		     (m_config->getStr("masterOutput")).Data(),
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Also clear the class data:
  m_yieldStorage.clear();
  m_sysStorage.clear();
  
  std::cout << "\nSystematicsTool::Initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the normalization systematic uncertainty from a given systematic
   source for a given sample.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
*/
double SystematicsTool::getNormSys(TString sysName, TString sampleName) {
  if (m_sysStorage.count(sysKey(sysName, sampleName)) > 0) {
    return m_sysStorage[sysKey(sysName, sampleName)];
  }
  else {
    std::cout << "SystematicsTool: ERROR! normalization systematic " << sysName 
	      << " not defined for sample " << sampleName << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Retieve the migration systematic uncertainty from a given systematic source 
   for a given sample in a specified category.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @param cateIndex - The index of the category.
*/
double SystematicsTool::getMigrSys(TString sysName, TString sampleName,
				 int cateIndex) {
  if (m_sysStorage.count(sysKey(sysName, sampleName, cateIndex)) > 0) {
    return m_sysStorage[sysKey(sysName, sampleName, cateIndex)];
  }
  else {
    std::cout << "SystematicsTool: ERROR! migration systematic " << sysName 
	      << " not defined for sample " << sampleName 
	      << " in category " << cateIndex << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a list of all of the defined systematic uncertainties.
   @return - A list of systematic variations.
*/
std::vector<TString> SystematicsTool::listAllSys() {
  return m_config->getStrV("SystematicsList");
}

/**
   -----------------------------------------------------------------------------
   Load all of the systematic variations for a given sample.
   @param sampleName - The name of the sample for which sys. will be loaded.
*/
void SystematicsTool::loadAllSys(TString sampleName) {
  // loop over systematics:
  std::vector<TString> sysList = listAllSys(); 
  for (int i_s = 0; i_s < (int)sysList.size(); i_s++) {
    loadSingleSys(sysList[i_s], sampleName);
  }
}

/**
   -----------------------------------------------------------------------------
   Load the values of a single systematic uncertainty for a single sample.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
*/
void SystematicsTool::loadSingleSys(TString sysName, TString sampleName) {
  TString inputDir = Form("%s/%s/DMMassPoints/Systematics", 
			  (m_config->getStr("masterOutput")).Data(),
			  (m_config->getStr("jobName")).Data());
  
  // First load cutflow to get yield:
  std::ifstream sysInput(Form("%s/cutflow_%s_%s.txt",inputDir.Data(), 
			      sysName.Data(), sampleName.Data()));
  TString str1; TString str2; double passYield; double totalYield;
  while (!sysInput.eof()) {
    sysInput >> str1 >> passYield >> str2 >> totalYield;
    if (str1.EqualTo("AllCuts")) {
      m_yieldStorage[sysKey(sysName, sampleName)] = passYield;
    }
  }
  
  // Then calculate systematic. If nominal has not already been created, create:
  if (m_yieldStorage.count(sysKey(sysName, sampleName)) == 0 &&
      !sysName.EqualTo("Nominal")) {
    loadSingleSys("Nominal", sampleName);
  }
  
  
  // Then load categorization (WARNING! improper output in DMMassPoints)
}

/**
   -----------------------------------------------------------------------------
*/
//std::vector<TString> SystematicsTool::rankMigrSysForSample(TString sampleName) {
//}

/**
   -----------------------------------------------------------------------------
   Get a list of systematic uncertainties ordered from largest effect to 
   smallest effect for the specified sample.
   @param sampleName - The name of the sample for which sys. will be loaded.
   @return - An ordered vector of systematic uncertainty names. 
*/
std::vector<TString> SystematicsTool::rankNormSysForSample(TString sampleName) {
  
  // An unordered vector of all systematics:
  std::vector<TString> unorderedSys = listAllSys();
  
  // A vector ordered from largest to smallest systematic effect:
  std::vector<TString> orderedSys; orderedSys.clear();
  
  // Loop over all the unordered systematics:
  for (int i_s = 0; i_s < (int)unorderedSys.size(); i_s++) {
    
    bool wasInserted = false;
    
    // Iterate over the ordered systematics:
    for (std::vector<TString>::iterator orderIter = orderedSys.begin();
	 orderIter != orderedSys.end(); orderIter++) {
      
      // Insert unordered val if the value of orderIter is lower:
      double currUnorderedVal = getNormSys(unorderedSys[i_s], sampleName);
      double currOrderedVal = getNormSys(*orderIter, sampleName);
      if (currUnorderedVal > currOrderedVal) {
	orderedSys.insert(orderIter, unorderedSys[i_s]);
	wasInserted = true;
      }
      // else continue iteration.
    }
    
    // if it was not inserted during iteration, push back on end:
    if (!wasInserted) orderedSys.push_back(unorderedSys[i_s]);
    
  }
  
  // Check that output list makes sense:
  if ((int)orderedSys.size() != (int)unorderedSys.size()) {
    std::cout << "SystematicsTool: ERROR! Something went wrong in ordering algo"
	      << std::endl;
    exit(0);
  }
  
  return orderedSys;
}

/**
   -----------------------------------------------------------------------------
   Prints a ranked list of the normalization systematics and their values. 
   @param sampleName - The name of the sample for which sys. will be loaded.
*/
void SystematicsTool::saveRankedNormSys(TString sampleName) {
  std::ofstream outputRanking(Form("%s/normSysRank_%s.txt", 
				   m_outputDir.Data(), sampleName.Data()));
  std::vector<TString> rankedSysNames = rankNormSysForSample(sampleName);
  for (std::vector<TString>::iterator iterSys = rankedSysNames.begin(); 
       iterSys != rankedSysNames.end(); iterSys++) {
    outputRanking << *iterSys << " \t" << getNormSys(*iterSys,sampleName)
		  << std::endl;
  }
  outputRanking.close();
}

/**
   -----------------------------------------------------------------------------
   A private method for accessing map data.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @param cateIndex - The index of the category.
*/
TString SystematicsTool::sysKey(TString sysName, TString sampleName,
				int cateIndex) {
  return Form("%s_%s_c%d", sysName.Data(), sampleName.Data(), cateIndex);
}

/**
   -----------------------------------------------------------------------------
   A private method for accessing map data.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
*/
TString SystematicsTool::sysKey(TString sysName, TString sampleName) {
  return Form("%s_%s", sysName.Data(), sampleName.Data());
}
