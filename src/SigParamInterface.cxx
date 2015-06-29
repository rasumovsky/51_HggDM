////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParamInterface.cxx                                               //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/06/2015                                                          //
//                                                                            //
//  This is an interface class for the SigParam class. It allows one to load  //
//  preexisting workspaces or create new ones if they cannot be loaded.       //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "SigParamInterface.h"

/**
   -----------------------------------------------------------------------------
   Initialize the SigParamInterface class with a new RooCategory.
   @param newJobName - The name of the job 
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
*/
SigParamInterface::SigParamInterface(TString newJobName, TString newCateScheme,
				     TString newOptions) {
  std::cout << "\nSigParamInterface::Initializing..."
	    << "\n\tjobName = " << newJobName
	    << "\n\tcateScheme = " << newCateScheme 
	    << "\n\toptions = " << newOptions << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  signalsOK = true;
  failedSigParam = "";
  sigMap.clear();
  
  CommonFunc::SetAtlasStyle();

  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMSigParam", DMAnalysis::masterOutput.Data(),
		   jobName.Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Load the SM signal parameterization from file or start from scratch:
  for (int i_SM = 0; i_SM < DMAnalysis::nSMModes; i_SM++) {
    if ((options.Contains("FromFile") && loadFile(DMAnalysis::sigSMModes[i_SM]))
	|| createNew(DMAnalysis::sigSMModes[i_SM])) {
      std::cout << "SigParamInterface: " << DMAnalysis::sigSMModes[i_SM] 
		<< " ready!" << std::endl;
    }
    else signalsOK = false;
  }
  
  // Load the DM signal parameterization from file or start from scratch:
  for (int i_DM = 0; i_DM < DMAnalysis::nDMModes; i_DM++) {
    if ((options.Contains("FromFile") && loadFile(DMAnalysis::sigDMModes[i_DM]))
	|| createNew(DMAnalysis::sigDMModes[i_DM])) {
      std::cout << "SigParamInterface: " << DMAnalysis::sigDMModes[i_DM]
		<< " ready!" << std::endl;
    }
    else signalsOK = false;
  }
  
  // Also load or create the total SM parameterization:
  
  if ((options.Contains("FromFile") && loadFile("SM")) || createNew("SM")) {
    std::cout << "SigParamInterface: Total SM signal ready!" << std::endl;
  }
  else signalsOK = false;
  
  if (signalsOK) {
    std::cout << "SigParamInterface: Successfully initialized!" << std::endl;
  }
  else {
    std::cout << "SigParamInterface: Problem initializing :(" << std::endl;
    std::cout << "\tFailed fits: " << failedSigParam << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Check if all of the signals were either loaded successfully or created.
*/
bool SigParamInterface::allSignalsReady() {
  if (!signalsOK) {
    std::cout << "SigParamInterface: Problems detected with following signals"
	      << failedSigParam << std::endl;
  }
  return signalsOK;
}

/**
   -----------------------------------------------------------------------------
   Create a new signal parameterization.
   @param signalType - The type of signal for parameterization.
   @returns - True iff successfully created.
*/
bool SigParamInterface::createNew(TString signalType) {
  std::cout << "SigParamInterface: Creating new signal for "
	    << signalType << std::endl;
  
  bool signalConverged = true;
  SigParam *sp = new SigParam(signalType, outputDir);
  for (int i_c = 0; i_c < DMAnalysis::getNumCategories(cateScheme); i_c++) {
    RooDataSet *currDataSet = getData(signalType, i_c);
    sp->addDataSet(DMAnalysis::higgsMass, i_c, currDataSet, "m_yy");
    
    if (sp->makeSingleResonance(DMAnalysis::higgsMass, i_c, 
				DMAnalysis::resonancePDF)) {
      sp->saveAll();
      sp->plotSingleResonance(DMAnalysis::higgsMass, i_c);
      sigMap[signalType] = sp;
    }
    else {
      signalConverged = false;
      failedSigParam += signalType + ", ";
    }
  }
  return signalConverged;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the relevant datasets for the signal and category of interest.
   @param signalType - The type of signal for parameterization.
   @param cateIndex - The category index.
   @returns - A RooDataSet with all data points necessary for fitting.
*/
RooDataSet* SigParamInterface::getData(TString signalType, int cateIndex) {
  if (signalType.EqualTo("SM")) {
    RooDataSet *currData = NULL;
    for (int i_SM = 0; i_SM < DMAnalysis::nSMModes; i_SM++) {
      DMMassPoints *mp = new DMMassPoints(jobName, DMAnalysis::sigSMModes[i_SM],
					  cateScheme, "FromFile", NULL);
      if (i_SM == 0) currData = mp->getCateDataSet(cateIndex);
      else currData->append(*(mp->getCateDataSet(cateIndex)));
    }
    return currData;
  }
  else {
    DMMassPoints *mp = new DMMassPoints(jobName, signalType, cateScheme, 
					"FromFile", NULL);
    return mp->getCateDataSet(cateIndex);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a pointer to the SigParam class for the specified signal type.
   @param signalType - The type of signal for parameterization.
*/
SigParam* SigParamInterface::getSigParam(TString signalType) {
  return sigMap[signalType];
}

/**
   -----------------------------------------------------------------------------
   Load signal parameterization from file. If that fails, create new.
   @param signalType - The type of signal for parameterization.
   @returns - True iff successfully loaded or created.
*/
bool SigParamInterface::loadFile(TString signalType) {
  SigParam *sp = new SigParam(signalType, outputDir);
  if (sp->loadParameterization(outputDir, signalType)) {
    sigMap[signalType] = sp;
    std::cout << "SigParamInterface: Successful load from file for "
	      << signalType << std::endl;
    return true;
  }
  else {
    std::cout << "SigParamInterface: Failed to load from file for "
	      << signalType << std::endl;
    return createNew(signalType);
  }
}
