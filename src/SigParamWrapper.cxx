////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParamWrapper.cxx                                                 //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/06/2015                                                          //
//                                                                            //
//  This program demonstrates how one can implement the SigParam.cxx signal   //
//  parameterization class. It relies on methods that are not provided, which //
//  the user must replace (for instance, DMMassPoints and CommonFunc.         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "SigParam.h"// the signal parameterization class
#include "DMMassPoints.h"// a specialized class for retrieving datasets
#include "CommonFunc.h"// a header with various functions.

using namespace DMAnalysis;// a header specific to the mono-H analysis.

/**
   This is the main method:
*/
int main (int argc, char **argv) {
  // Check arguments:
  if (argc < 1) {
    printf("\nUsage: %s\n\n", argv[0]);
    exit(0);
  }
  
  // Just a function to set the ATLAS style
  CommonFunc::SetAtlasStyle();
  
  // For this example, we only have one mass point and one category:
  double resonanceMass = 125.0;
  int categoryIndex = 0;
  
  //----------------------------------------//
  // EXAMPLE 1: Fit the inclusive SM shape by inputting RooDataSets.
  //----------------------------------------//
  
  // Step 1: Instantiate the signal parameterization with no special options:
  SigParam *sp_inclusive = new SigParam("");
  
  // Loop over all of the SM production modes:
  for (int i_SM = 0; i_SM < DMAnalysis::nSMModes; i_SM++) {
 
    // DMMassPoints is a tool for loading the mono-H datasets:
    DMMassPoints *mp = new DMMassPoints("TestingSigParam",
					DMAnalysis::sigSMModes[i_SM],
					"inclusive", "FromFile", NULL);
    RooDataSet *currDataSet = mp->getCateDataSet(categoryIndex);
    
    // Step 2: Import the data for each production mode into the SigParam tool:
    sp_inclusive->addDataSet(resonanceMass, categoryIndex, currDataSet, "m_yy");
  }
  
  // Step 3: Fit the resonance:
  bool converged_inclusive = sp_inclusive->makeSingleResonance(resonanceMass, categoryIndex, "DoubleCB");
  //bool converged_inclusive = sp_inclusive->makeSingleResonance(resonanceMass, categoryIndex, "CBGA");
  
  if (converged_inclusive) {
    std::cout << "SigParamWrapper: fit converged!" << std::endl;
  }
  else {
    std::cout << "SigParamWrapper: fit did not converge :(" << std::endl;
  }
  
  // Step 4: Plot the result:
  sp_inclusive->plotSingleResonance(resonanceMass, categoryIndex,
				    "plot_inclusive.eps");
  // Step 5: save the parameterization:
  sp_inclusive->saveParameterization("sigParamFile_inclusive.root");
  
  //----------------------------------------//
  // EXAMPLE 2: Fit a DM shape by inputting single events.
  //----------------------------------------//
  
  // Step 1: Instantiate the signal parameterization with no special options:
  SigParam *sp_DM = new SigParam("");
  
  // DMMassPoints is just a tool for loading the mono-H datasets:
  DMMassPoints *mp = new DMMassPoints("TestingSigParam",
				      DMAnalysis::sigDMModes[0],
				      "inclusive", "FromFile", NULL);
  RooDataSet *currDataSet = mp->getCateDataSet(categoryIndex);
  // Loop over dataset entries:
  for (int i_d = 0; i_d < currDataSet->numEntries(); i_d++) {
    // Get the RooRealVars and values in the current event:
    const RooArgSet *currArgs = (RooArgSet*)currDataSet->get(i_d);
    
    double massValue = -999.0; double weightValue = -999.0;
    
    // Iterate through the RooArgSet, find mass and weight variables:
    TIterator *iterArgs = currArgs->createIterator();
    RooRealVar* currArg;
    while ((currArg = (RooRealVar*)iterArgs->Next())) {
      if (TString(currArg->GetName()).EqualTo("m_yy")) {
	massValue = currArg->getVal();
      }
    }
    weightValue = currDataSet->weight();
    
    // Step 2: Add the mass and weight values to the dataset:
    sp_DM->addMassPoint(resonanceMass, categoryIndex, massValue, weightValue);
  }
  
  // Step 3: Fit the resonance:
  bool converged_DM = sp_DM->makeSingleResonance(resonanceMass,
						 categoryIndex,
						 "CBGA");
  if (converged_DM) {
    std::cout << "SigParamWrapper: fit converged!" << std::endl;
  }
  else {
    std::cout << "SigParamWrapper: fit did not converge :(" << std::endl;
  }
  
  // Step 4: Plot the result:
  sp_DM->plotSingleResonance(resonanceMass, categoryIndex, "plot_DM.eps");
  
  // Step 5: save the parameterization:
  sp_DM->saveParameterization("sigParamFile_DM.root");
  
  return 0;
}
