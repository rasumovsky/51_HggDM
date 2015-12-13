////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SystematicsTool.cxx                                                 //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/12/2015                                                          //
//                                                                            //
//  This class loads and stores information on systematic uncertainties. For  //
//  the moment, these include normalization and migration systematics that    //
//  are calculated using the output of DMMassPoints.cxx. Eventually, the      //
//  framework can be extended to include resolution and mass scale systematic //
//  uncertainties, if useful.                                                 //
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
   Calculate the size of a migration systematic uncertainty.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @param cateIndex - The index of the category.
   @return - The fractional size of a systematic effect.
*/
double SystematicsTool::calculateMigrSys(TString sysName, TString sampleName,
					 int cateIndex) {
  // Normalize by the total event yield when calculating migration effect.
  double normFactor = (getYield("Nominal", sampleName) / 
		       getYield(sysName, sampleName));
  
  double sysValue = ((normFactor * getYield(sysName, sampleName, cateIndex) -
		      getYield("Nominal", sampleName, cateIndex)) / 
		     getYield("Nominal", sampleName, cateIndex));

  // Store the result and return it:
  setMigrSys(sysName, sampleName, cateIndex, sysValue);
  return sysValue;
}

/**
   -----------------------------------------------------------------------------
   Calculate the size of a normalization systematic uncertainty.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @return - The fractional size of a systematic effect.
*/
double SystematicsTool::calculateNormSys(TString sysName, TString sampleName) {
  double sysValue
    = ((getYield(sysName,sampleName) - getYield("Nominal",sampleName)) / 
       getYield("Nominal", sampleName));
  // Store the result and return it:
  setNormSys(sysName, sampleName, sysValue);
  return sysValue;
}

/**
   -----------------------------------------------------------------------------
   Retieve the migration systematic uncertainty from a given systematic source 
   for a given sample in a specified category.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @param cateIndex - The index of the category.
   @return - The fractional size of the systematic effect.
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
   Retrieve the normalization systematic uncertainty from a given systematic
   source for a given sample.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @return - The fractional size of the systematic effect.
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
   Retieve the total yield for a given sample in a specified category with a 
   specified systematic variation.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @return - The event yield.
*/
double SystematicsTool::getYield(TString sysName, TString sampleName) {
  if (m_yieldStorage.count(sysKey(sysName, sampleName)) > 0) {
    return m_yieldStorage[sysKey(sysName, sampleName)];
  }
  else {
    std::cout << "SystematicsTool: ERROR! yield undefined for sys. " << sysName 
	      << " and sample " << sampleName << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Retieve the yield in a given category for a given sample with a specified 
   systematic variation.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample for which sys. will be loaded.
   @param cateIndex - The index of the category.
   @return - The event yield.
*/
double SystematicsTool::getYield(TString sysName, TString sampleName,
				 int cateIndex) {
  if (m_yieldStorage.count(sysKey(sysName, sampleName, cateIndex)) > 0) {
    return m_yieldStorage[sysKey(sysName, sampleName, cateIndex)];
  }
  else {
    std::cout << "SystematicsTool: ERROR! yield undefined for sys. " << sysName 
	      << " and sample " << sampleName << std::endl;
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
  
  // If nominal has not already been created, create (necessary for sys. calc):
  if (m_yieldStorage.count(sysKey("Nominal", sampleName)) == 0 &&
      !sysName.EqualTo("Nominal")) {
    loadSingleSys("Nominal", sampleName);
  }
  
  // Step 1: load cutflow to get normalization systematic:
  TString normFileName = Form("%s/cutflow_%s_%s.txt",inputDir.Data(), 
			      sysName.Data(), sampleName.Data());
  std::ifstream normSysInput(normFileName);
  if (!normSysInput.is_open()) {
    std::cout << "SystematicsTool: ERROR opening norm sys. file: "
	      << normFileName << std::endl;
    exit(0);
  }
  
  // Load the event yield after all cuts:
  TString str1; TString str2; double passYield; double totalYield;
  while (!normSysInput.eof()) {
    normSysInput >> str1 >> passYield >> str2 >> totalYield;
    if (str1.EqualTo("AllCuts")) {
      setYield(sysName, sampleName, passYield);
    }
  }
  normSysInput.close();
  
  // Calculate the systematic effect on normalization:
  calculateNormSys(sysName, sampleName);
  
  // Step 2: load categorization to get migration systematic:
  TString migrFileName = Form("%s/categorization_%s_%s.txt",inputDir.Data(), 
			      sysName.Data(), sampleName.Data());
  std::ifstream migrSysInput(migrFileName);
  if (!migrSysInput.is_open()) {
    std::cout << "SystematicsTool: ERROR opening migration sys. file: "
	      << migrFileName << std::endl;
    exit(0);
  }
  
  // Load the event categorization after all cuts:
  std::string line;
  while (!migrSysInput.eof()) {
    std::getline(migrSysInput, line);
    TString currLine = TString(line);
    TObjArray *tokenizedLine = currLine.Tokenize(" ");
    for (int i_e = 1; i_e < tokenizedLine->GetEntries(); i_e++) {
      TString *currValue = (TString*)(tokenizedLine->At(i_e));
      setYield(sysName, sampleName, i_e-1, currValue->Atof());
      delete currValue;
    }
    delete tokenizedLine;
  }
  migrSysInput.close();
  
  // Calculate the systematic effect of migration:
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    calculateMigrSys(sysName, sampleName, i_c);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a list of migration systematic uncertainties ordered from largest to 
   smallest effect for the specified sample.
   @param sampleName - The name of the sample.
   @return - An ordered vector of migration systematic uncertainty names. 
*/
std::vector<TString> SystematicsTool::rankMigrSysForSample(TString sampleName) {
  
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
      
      // Find the maximum shift:
      double maxUnorderedVal = getMigrSys(unorderedSys[i_s], sampleName, 0); 
      double maxOrderedVal = getMigrSys(*orderIter, sampleName, 0);
      for (int i_c = 1; i_c < m_config->getInt("nCategories"); i_c++) {
	double currUnorderedVal = getMigrSys(unorderedSys[i_s], sampleName,i_c);
	double currOrderedVal = getMigrSys(*orderIter, sampleName, i_c);
	if (currUnorderedVal > maxUnorderedVal) {
	  maxUnorderedVal = currUnorderedVal;
	}
	if (currOrderedVal > maxOrderedVal) {
	  maxOrderedVal = currOrderedVal;
	}
      }
      
      // Insert unordered val if the value of orderIter is lower:
      if (maxUnorderedVal > maxOrderedVal) {
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
   Get a list of systematic uncertainties ordered from largest effect to 
   smallest effect for the specified sample.
   @param sampleName - The name of the sample.
   @return - An ordered vector of normalization systematic uncertainty names. 
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
   Prints a ranked list of the migration systematics and their values. 
   @param sampleName - The name of the sample.
*/
void SystematicsTool::saveRankedMigrSys(TString sampleName) {
  std::ofstream outputRanking(Form("%s/migrSysRank_%s.txt", 
				   m_outputDir.Data(), sampleName.Data()));
  std::vector<TString> rankedSysNames = rankMigrSysForSample(sampleName);
  for (std::vector<TString>::iterator iterSys = rankedSysNames.begin(); 
       iterSys != rankedSysNames.end(); iterSys++) {
    outputRanking << *iterSys;
    for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
      outputRanking << " " << getMigrSys(*iterSys, sampleName, i_c);
    }
    outputRanking << std::endl;
  }
  outputRanking.close();
}

/**
   -----------------------------------------------------------------------------
   Prints a ranked list of the normalization systematics and their values. 
   @param sampleName - The name of the sample.
*/
void SystematicsTool::saveRankedNormSys(TString sampleName) {
  std::ofstream outputRanking(Form("%s/normSysRank_%s.txt", 
				   m_outputDir.Data(), sampleName.Data()));
  std::vector<TString> rankedSysNames = rankNormSysForSample(sampleName);
  for (std::vector<TString>::iterator iterSys = rankedSysNames.begin(); 
       iterSys != rankedSysNames.end(); iterSys++) {
    outputRanking << *iterSys << " " << getNormSys(*iterSys, sampleName)
		  << std::endl;
  }
  outputRanking.close();
}

/**
   -----------------------------------------------------------------------------
   Set the value of the systematic uncertainty on normalization.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample.
   @param sysValue - The new value of the systematic uncertainty.
*/
void SystematicsTool::setNormSys(TString sysName, TString sampleName, 
				 double sysValue) {
  m_sysStorage[sysKey(sysName, sampleName)] = sysValue;
}

/**
   -----------------------------------------------------------------------------
   Set the value of the systematic uncertainty on migration.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample.
   @param cateIndex - The index of the category.
   @param sysValue - The new value of the systematic uncertainty.
*/
void SystematicsTool::setMigrSys(TString sysName, TString sampleName,
				 int cateIndex, double sysValue) {
  m_sysStorage[sysKey(sysName, sampleName, cateIndex)] = sysValue;
}
/**
   -----------------------------------------------------------------------------
   Set the value of the yield.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample.
   @param yieldValue - The new value of the yield.
*/
void SystematicsTool::setYield(TString sysName, TString sampleName, 
				 double yieldValue) {
  m_yieldStorage[sysKey(sysName, sampleName)] = yieldValue;
}

/**
   -----------------------------------------------------------------------------
   Set the value of the yield.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample.
   @param cateIndex - The index of the category.
   @param yieldValue - The new value of the yield.
*/
void SystematicsTool::setYield(TString sysName, TString sampleName,
				 int cateIndex, double yieldValue) {
  m_yieldStorage[sysKey(sysName, sampleName, cateIndex)] = yieldValue;
}

/**
   -----------------------------------------------------------------------------
   A private method for accessing map data.
   @param sysName - The name of the systematic uncertainty. 
   @param sampleName - The name of the sample.
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
   @param sampleName - The name of the sample.
*/
TString SystematicsTool::sysKey(TString sysName, TString sampleName) {
  return Form("%s_%s", sysName.Data(), sampleName.Data());
}
