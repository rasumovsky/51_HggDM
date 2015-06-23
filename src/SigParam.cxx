////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/06/2015                                                          //
//                                                                            //
//  This class implements the resonance modeling for the ATLAS Hgamma group.  //
//
//  - "DoubleCB" or "CBGA" for the shape. Defaults to CBGA.                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/////////////////////
// DEV NOTES
//
// 
// METHODS:
//   getFileName()
//   loadFromFile()
//   printParameterValues()
//   printTable()
//   saveParameterValues()
//
/////////////////////

#include "SigParam.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the SigParam class.
   @param options - the parameterization options (see documentation).
*/
SigParam::SigParam(TString options) {
  std::cout << "\nSigParam::Initializing..."
    	    << "\n\toptions = " << options << std::endl;
  
  // Assign member variables:
  m_options = options;
  
  // Create RooWorkspace and RooCategory to store all signal models:
  m_ws = new RooWorkspace("signalWS");
  m_cat = new RooCategory("signalCates","signalCates");
  
  // Import mass and weight variables, then make pointers:
  m_ws->factory(Form("m_yy[%f,%f]", 10.0, 7000.0));
  m_ws->factory(Form("wt[%f]", 1.0));
  m_yy = m_ws->var("m_yy");
  m_wt = m_ws->var("wt");
  
  // Define the data sets at each mass in each category for resonance fit:
  m_massCatePairs.clear();

  // Set the default functions for each variable's parameterization
  m_funcList.clear();
  setVarParameterization("muCBNom", "pol2");
  setVarParameterization("sigmaCBNom", "pol1");
  setVarParameterization("alphaCB", "pol1");
  setVarParameterization("nCB", "pol0");
  setVarParameterization("sigmaGANom", "pol1");
  setVarParameterization("fracCB", "pol1");
  setVarParameterization("alphaCBLo", "pol2");
  setVarParameterization("nCBLo", "pol0");
  setVarParameterization("alphaCBHi", "pol1");
  setVarParameterization("nCBHi", "pol0");
  
  // Lists of mass resolution and mass scale systematics:
  m_listMRS = "";
  m_listMSS = "";
  
  std::cout << "SigParam: Successfully initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Adds the data from a RooDataSet object to the dataset for parameterization.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the event category index (starting at 0).
   @param dataSet - the RooDataSet object to import.
   @param observableName - the name of the observable.
*/
void SigParam::addDataSet(double resonanceMass, int cateIndex,
			  RooDataSet* dataSet, TString observableName) {
  
  // Loop over events stored in the input dataset:
  for (int i_d = 0; i_d < dataSet->numEntries(); i_d++) {
    // Get the RooRealVars and values in the current event:
    const RooArgSet *currArgs = (RooArgSet*)dataSet->get(i_d);
    
    double massValue = -999.0; double weightValue = -999.0;
    
    // Iterate through the RooArgSet, find mass and weight variables:
    TIterator *iterArgs = currArgs->createIterator();
    RooRealVar* currArg;
    while ((currArg = (RooRealVar*)iterArgs->Next())) {
      if (TString(currArg->GetName()).EqualTo(observableName)) {
	massValue = currArg->getVal();
      }
    }
    weightValue = dataSet->weight();
    // Then add the mass and weight values to the dataset:
    if (massValue > -990 && weightValue > -990) {
      addMassPoint(resonanceMass, cateIndex, massValue, weightValue);
    }
    // Exit if the variable names were incorrect:
    else {
      std::cout << "SigParam: No match for observable or weight." << std::endl;
      exit(0);
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Adds the data from a Root TTree object or an MxAOD to the dataset for
   signal parameterization. Note: the input 
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the event category index (starting at 0).
   @param dataTree - the TTree or MxAOD object to import.
   @param massBranchName - the name of the branch storing mass values.
   @param weightBranchName - the name of the branch storing event weights.
*/
void SigParam::addDataTree(double resonanceMass, int cateIndex,
			   TTree* dataTree, TString massBranchName,
			   TString weightBranchName) {
  // Get the mass and weight branches.
  double massValue;
  double weightValue;
  dataTree->SetBranchAddress(massBranchName, &massValue);
  dataTree->SetBranchAddress(weightBranchName, &weightValue);
  
  // Loop over the TTree.
  for (Long64_t i_t = 0; i_t < dataTree->GetEntries(); i_t++) {
    
    // Add current event values to the dataset used for parameterization.
    addMassPoint(resonanceMass, cateIndex, massValue, weightValue);
  }
}

/**
   -----------------------------------------------------------------------------
   Adds a mass point to the dataset for fitting the signal diphoton resonance.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the event category index (starting at 0).
   @param diphotonMass - the reconstructed diphoton invariant mass.
   @param eventWeight - the weight of the event.
*/
void SigParam::addMassPoint(double resonanceMass, int cateIndex, 
			    double diphotonMass, double eventWeight) {

  TString currKey = getKey(resonanceMass, cateIndex);
  
  // Create new dataset if the corresponding one doesn't yet exist.
  if (!dataExists(resonanceMass, cateIndex)) {
    
    RooDataSet* newData 
      = new RooDataSet(Form("data_%s",currKey.Data()),
		       Form("data_%s",currKey.Data()),
		       RooArgSet(*m_yy,*m_wt), RooFit::WeightVar(*m_wt));
    m_ws->import(*newData);
    
    // Each dataset also corresponds to a unique category in the fit:
    m_cat->defineType(currKey);
    
    // Keep track of all resonance mass - category pairs for fitting.
    std::pair<double,int> currPair;
    currPair.first = resonanceMass;
    currPair.second = cateIndex;
    m_massCatePairs.push_back(currPair);
  }
  
  // Set the observable and weight values and then fill dataset:
  m_yy->setVal(diphotonMass);
  m_wt->setVal(eventWeight);
  RooDataSet* currData =(RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
  currData->add(RooArgSet(*m_yy,*m_wt), eventWeight);
}

/**
   -----------------------------------------------------------------------------
   Add a single mass resolution systematic uncertainty to the signal shape.
   Note: the constraint terms must be defined separately. 
   @param nameMResSys - the name of the systematic nuisance parameter.
*/
void SigParam::addMResSystematic(TString nameMResSys) {

  // Name of nuisance parameter:
  TString atlasExpMRS = Form("atlas_expected_%s",nameMResSys.Data());
  
  // Check to see if the nuisance parameter is already in the workspace:
  if (!(bool)m_ws->obj(atlasExpMRS)) {
    m_ws->factory(Form("%s[1]",atlasExpMRS.Data()));
  }
  
  // Add to the list storing signal resolution nuisance parameters:
  if (!m_listMRS.Contains(nameMResSys)) {
    m_listMRS.Append(Form(",%s",atlasExpMRS.Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Add several mass resolution systematic uncertainties to the signal shape. 
   Note: the constraint terms must be defined separately. 
   @param namesMResSys - list of the names of systematic nuisance parameters.
*/
void SigParam::addMResSystematics(std::vector<TString> namesMResSys) {
  for (std::vector<TString>::iterator sysIter = namesMResSys.begin();
       sysIter != namesMResSys.end(); sysIter++) {
    SigParam::addMResSystematic(*sysIter);
  }
}

/**
   -----------------------------------------------------------------------------
   Add a single mass scale systematic uncertainty to the signal shape.
   Note: the constraint terms must be defined separately. 
   @param nameMScaleSys - the name of the systematic nuisance parameter.
*/
void SigParam::addMScaleSystematic(TString nameMScaleSys) {
  
  // Name of nuisance parameter:
  TString atlasExpMSS = Form("atlas_expected_%s", nameMScaleSys.Data());
  
  // Check to see if the nuisance parameter is already in the workspace:
  if (!(bool)m_ws->obj(atlasExpMSS)) {
    m_ws->factory(Form("%s[1]",atlasExpMSS.Data()));
  }
  
  // Add to the list storing signal mass scale nuisance parameters:
  if (!m_listMSS.Contains(nameMScaleSys)) {
    m_listMSS.Append(Form(",%s",atlasExpMSS.Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Add several mass scale systematic uncertainties to the signal shape. 
   Note: the constraint terms must be defined separately. 
   @param namesMScaleSys - list of the names of systematic nuisance parameters.
*/
void SigParam::addMScaleSystematics(std::vector<TString> namesMScaleSys) {
  for (std::vector<TString>::iterator sysIter = namesMScaleSys.begin();
       sysIter != namesMScaleSys.end(); sysIter++) {
    SigParam::addMScaleSystematic(*sysIter);
  }
}

/**
   -----------------------------------------------------------------------------
   Adds a signal PDF from this class to a pre-existing workspace. 
   @param workspace - the pre-existing workspace.
   @param nuisParams - a RooArgSet containing fit parameters (set NULL if none).
   @param cateIndex - the index of the category of the desired PDF.
*/
void SigParam::addSigToWS(RooWorkspace *&workspace, RooArgSet *&nuisParams, 
			  int cateIndex) {
  std::cout << "SigParam: Adding parameterized signal in category "
	    << cateIndex << " to pre-existing workspace." << std::endl;
  
  // Add the signal model to the workspace:
  RooAbsPdf* currSignal = m_ws->pdf(Form("sigPdf_c%d",cateIndex));
  workspace->import(*currSignal);
  
  // Then explicitly add the parameters to the workspace:
  RooArgSet *currSet = currSignal->getVariables();
  TIterator *iterParam = currSet->createIterator();
  RooRealVar* currParam = NULL;
  while ((currParam = (RooRealVar*)iterParam->Next())) {
    workspace->import(*currParam);
    if (nuisParams) {
      nuisParams->add(*currParam);
    }
  }
  std::cout << "SigParam:: Finished adding parameterized signal." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Adds a signal PDF from this class to a pre-existing workspace. 
   @param workspace - the pre-existing workspace.
   @param nuisParams - a RooArgSet containing fit parameters (set NULL if none).
   @param resonanceMass - the mass of the resonance.
   @param cateIndex - the index of the category of the desired PDF.
*/
void SigParam::addSigToWS(RooWorkspace *&workspace, RooArgSet *&nuisParams,
			  double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Adding parameterized signal in category "
	    << cateIndex << " with mass " << resonanceMass 
	    << " to pre-existing workspace." << std::endl;
  
  // Add the signal model to the workspace:
  RooAbsPdf* currSignal = m_ws->pdf(Form("sigPdf_%s",(getKey(resonanceMass, cateIndex)).Data()));
  workspace->import(*currSignal);
  
  // Then explicitly add the parameters to the workspace:
  RooArgSet *currSet = currSignal->getVariables();
  TIterator *iterParam = currSet->createIterator();
  RooRealVar* currParam = NULL;
  while ((currParam = (RooRealVar*)iterParam->Next())) {
    workspace->import(*currParam);
    if (nuisParams) {
      nuisParams->add(*currParam);
    }
  }
  std::cout << "SigParam:: Finished adding parameterized signal." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Get a list of categories corresponding to a single mass point.
   @param resonanceMass - the mass value.
   @returns - a vector of category indices.
*/
std::vector<int> SigParam::categoriesForMass(double resonanceMass) {
  // Create a list of mass points in this category.
  std::vector<int> currCategories; currCategories.clear();
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if (equalMasses((m_massCatePairs[i_p]).first, resonanceMass)) {
      currCategories.push_back((m_massCatePairs[i_p]).second);
    }
  }
  std::sort(currCategories.begin(), currCategories.end());
  return currCategories;
}

/*
RooFormulaVar* SigParam::constructFormulaVar(TString varName, int polyOrder, double resonanceMass){
  
  TString ALF[6] = {"a","b","c","d","e","f"};
  
  TString varName = Form("%s_%s", varName.Data(), currKey.Data());
  
  TString funcForm = "@0";
  for (int i_v = 0; i_v < polyOrder; i_v++) {
    
    // Define a new RooRealVar in the workspace:
    m_ws->factory(Form("%s_%s_c%d[-0.38,-1.0,1.0]", ALF[i_v].Data(), varName.Data(), cateIndex));
    argList->add(*m_ws->var(Form("a_%s_c%d", varName.Data(), cateIndex)));
    funcForm += Form("+@%d",i_v+1);
    for (int i_p = 0; i_p < i_v+1; i_p++) {
      funcForm += Form("*%f",regularizedMass(currMass));
    }
  }
  funcForm += Form("%f",regularizedMass(currMass));
  
  RooFormulaVar *currVar = new RooFormulaVar(varName, funcForm, argList);

  m_ws->factory(Form("a_muCB_c%d[-0.38,-1.0,1.0]",cateIndex));
  m_ws->factory(Form("b_muCB_c%d[-0.06,-0.1,0.1]",cateIndex));
  m_ws->factory(Form("c_muCB_c%d[-0.02,-0.1,0.1]",cateIndex));
  
  RooFormulaVar *muCBNom = new RooFormulaVar(Form("muCBNom_%s",currKey.Data()), Form("@0+@1*%f+@2*%f+%f",regularizedMass(currMassPoints[i_m]),regularizedMass(currMassPoints[i_m]),currMassPoints[i_m]), RooArgList(*m_ws->var(Form("a_muCB_c%d",cateIndex)), *m_ws->var(Form("b_muCB_c%d",cateIndex)), *m_ws->var(Form("c_muCB_c%d",cateIndex))));

}
*/
	
/**
   -----------------------------------------------------------------------------
   Check if the dataset being requested has been instantiated (exists in map).
   @param massIndex - the index of the signal mass.
   @param cateIndex - the index of the category.
   @returns - true iff the dataset has been defined.
*/
bool SigParam::dataExists(double resonanceMass, int cateIndex) {
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if (equalMasses((m_massCatePairs[i_p]).first, resonanceMass) &&
	m_massCatePairs[i_p].second == cateIndex) {
      return true;
    }
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Check if two doubles are equal.
   @param massValue1 - the first mass value to compare.
   @param massValue2 - the second mass value to compare.
   @returns - true iff the masses are equal within 0.001 GeV.
*/
bool SigParam::equalMasses(double massValue1, double massValue2) {
  return (fabs(massValue1 - massValue2) <= 0.001);// mass precision in GeV
}

/**
   -----------------------------------------------------------------------------
   Perform a single or simultaneous fit.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the index of the category.
   @returns - the RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(double resonanceMass, int cateIndex) {
  std::cout << "SigParam: Preparing to fit resonance" << std::endl;
  TString currKey = (resonanceMass < 0.0) ? 
    Form("c%d",cateIndex) : getKey(resonanceMass,cateIndex);
  RooAbsPdf *currSignal = m_ws->pdf(Form("sigPdf_%s",currKey.Data()));
  RooAbsData *currData = m_ws->data(Form("data_%s",currKey.Data()));
  
  // Free, fit, fix, then save the nuisance parameter values:
  SigParam::setParamsConstant(currSignal, false);
  RooFitResult *result = currSignal->fitTo(*currData, RooFit::PrintLevel(0),
					   RooFit::SumW2Error(kTRUE),
					   RooFit::Save(true));
  SigParam::setParamsConstant(currSignal, true);
  std::cout << "SigParam: Fit procedure concluded." << std::endl;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Perform a simultaneous fit across multiple masses.
   @param cateIndex - the index of the category.
   @returns - the RooFitResult, which gives fit status.
*/
RooFitResult* SigParam::fitResult(int cateIndex) {
  return SigParam::fitResult(-999.9, cateIndex);
}

/**
   -----------------------------------------------------------------------------
   Retrieve the map key name for the dataset map.
   @param resonanceMass - the floating value of the signal mass.
   @param cateIndex - the index of the category.
   @returns - a key string for the dataset map.
*/
TString SigParam::getKey(double resonanceMass, int cateIndex) {
  TString key = Form("m%d_c%d", massDoubleToInt(resonanceMass), cateIndex);
  return key;
}

/**
   -----------------------------------------------------------------------------
   Get the number of categories contained in the datasets for fitting. Note:
   it is possible that there are different numbers of categories defined for 
   different mass points. This is up to the user to sort out. 
   @returns - the total number of categories for the parameterization.
*/
int SigParam::getNCategories() {
  m_nCategories = 0;
  // Loop over mass-category pairs, find highest category index. 
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    // The number of categories is equal to the highest index + 1:
    if ((m_massCatePairs[i_p]).second > m_nCategories + 1) {
      m_nCategories = (m_massCatePairs[i_p]).second + 1;
    }
  }
  return m_nCategories;
}

/**
   -----------------------------------------------------------------------------
   Get the value of the fit error for a particular parameter of the signal PDF. 
   @param paramName - the name of the shape parameter of interest.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - The index of the category.
   @returns - the value of the specified signal parameter. 
*/
double SigParam::getParameterError(TString paramName, double resonanceMass,
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("fracCB")) {
    return SigParam::getParameterError(paramName, cateIndex);
  }
  else {
    RooRealVar *var = m_ws->var(Form("%s_%s",paramName.Data(),
				     (getKey(resonanceMass,cateIndex)).Data()));
    if (!var) {
      std::cout << "SigParam: requested signal parameter not found: " 
		<< paramName << std::endl;
      return 0.0;
    }
    else {
      return var->getError();
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of the fit error for a particular parameter of the signal PDF. 
   @param paramName - the name of the shape parameter of interest.
   @param cateIndex - The index of the category.
   @returns - the value of the specified signal parameter. 
*/
double SigParam::getParameterError(TString paramName, int cateIndex) {
  RooRealVar *var = m_ws->var(Form("%s_c%d",paramName.Data(),cateIndex));
  if (!var) {
    std::cout << "SigParam: requested signal parameter not found: "
	      << paramName << std::endl;
    return 0.0;
  }
  else {
    return var->getError();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular parameter of the signal PDF. 
   @param paramName - the name of the shape parameter of interest.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - The index of the category.
   @returns - the value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, double resonanceMass, 
				   int cateIndex) {
  if (paramName.Contains("nCB") || paramName.Contains("fracCB")) {
    return SigParam::getParameterValue(paramName, cateIndex);
  }
  else {
    RooRealVar *var = m_ws->var(Form("%s_%s",paramName.Data(),
				     (getKey(resonanceMass,cateIndex)).Data()));
    if (!var) {
      std::cout << "SigParam: requested signal parameter not found: "
		<< paramName << std::endl;
      return 0.0;
    }
    else {
      return var->getVal();
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Get the value of a particular parameter of the signal PDF. 
   @param paramName - the name of the shape parameter of interest.
   @param cateIndex - The index of the category.
   @returns - the value of the specified signal parameter. 
*/
double SigParam::getParameterValue(TString paramName, int cateIndex) {
  RooRealVar *var = m_ws->var(Form("%s_c%d",paramName.Data(),cateIndex));
  if (!var) {
    std::cout << "SigParam: requested signal parameter not found: "
	      << paramName << std::endl;
    return 0.0;
  }
  else {
    return var->getVal();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the resonance shape for a single category and mass.
   @param resonanceMass - the truth mass of the resonance
   @param cateIndex - The index of the category for which we want the PDF.
   @returns - A pointer to the signal PDF.
*/
RooAbsPdf* SigParam::getSingleResonance(double resonanceMass, int cateIndex) {
  RooAbsPdf* pdf = m_ws->pdf(Form("sigPdf_%s",
				  (getKey(resonanceMass,cateIndex)).Data()));
  return pdf;
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular mass in a particular category.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - The index of the category for which we want the PDF.
   @returns - The signal yield for the specified mass in the given category.
*/
double SigParam::getYieldInCategory(double resonanceMass, int cateIndex) {
  // maybe add a protection in case the dataset is not defined.
  if (dataExists(resonanceMass, cateIndex)) {
    return ((RooAbsData*)m_ws->data(Form("data_%s",(getKey(resonanceMass,cateIndex)).Data())))->sumEntries();
  }
  else {
    std::cout << "SigParam: requested dataset not found." << std::endl;
    return 0.0;
  }
}

/**
   -----------------------------------------------------------------------------
   Get the signal yield for a particular resonance mass in all categories.
   @param resonanceMass - the truth mass of the resonance.
   @returns - The signal yield in all categories for the specified mass.
*/
double SigParam::getYieldTotal(double resonanceMass) {
  // Create a list of categories for this mass point:
  std::vector<int> currCategories = SigParam::categoriesForMass(resonanceMass);
  // Loop through names of datasets, add components:
  double sum = 0.0;
  for (int i_c = 0; i_c < (int)currCategories.size(); i_c++) {
    sum += getYieldInCategory(resonanceMass, currCategories[i_c]);
  }
  return sum;
}

/**
   -----------------------------------------------------------------------------
   Load the signal parameterization from file.
   @param fileName - the name of the .root file containing input workspace.
   @returns - true iff the file is successfully loaded.
*/
bool SigParam::loadParameterization(TString fileName) {
  std::cout << "SigParam: Load parameterization file " << fileName << std::endl;
  bool parameterizationExists = false;
  TFile inputFile(fileName);
  if (inputFile.IsOpen()) {
    m_ws = (RooWorkspace*)inputFile.Get("signalWS");
    if (m_ws) parameterizationExists = true;
    m_yy = m_ws->var("m_yy");
    m_wt = m_ws->var("wt");
    std::cout << "SigParam: Successfully loaded from file!" << std::endl;
  }
  return parameterizationExists;
}

/**
   -----------------------------------------------------------------------------
   Parameterize the resonance shape in all categories.
   @param function - the functional form of the resonance.
   @returns - true iff. all fits converge.
*/
bool SigParam::makeAllParameterizations(TString function) {
  std::cout << "SigParam: Engage full signal parameterization!" << std::endl;
  
  bool result = true;
  // Define models in each category independently:
  for (int i_c = 0; i_c < getNCategories(); i_c++) {
    if (!makeCategoryParameterization(i_c, function)) result = false;
  }
  
  std::cout << "SigParam: Complete full signal parameterization!" << std::endl;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Parameterize the resonance shape as a function of mass for a single category.
   @param cateIndex - the index of the category to fit.
   @param function - the functional form of the resonance.
   @returns - true iff. all fits converge.
*/
bool SigParam::makeCategoryParameterization(int cateIndex, TString function) {
  std::cout << "SigParam: parameterizing category " << cateIndex << std::endl;
  
  // Create a list of mass points in this category.
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
    
  // Create a simultaneous PDF just for this category:
  RooSimultaneous *currSim = new RooSimultaneous(Form("sigPdf_c%d",cateIndex),
						 Form("sigPdf_c%d",cateIndex),
						 *m_cat);
  
  // If no mass points, cannot fit signal, silly!
  if (currMassPoints.size() == 0) {
    std::cout << "SigParam: No masspoints for cate. " << cateIndex << std::endl;
    return false;
  }
  // If only one mass point, no need for parameterization:
  else if (currMassPoints.size() == 1) {
    std::cout << "SigParam: 1 mass point -> no parameterization." << std::endl;
    return makeSingleResonance(currMassPoints[0], cateIndex, function);
  }
  // If more than 1 mass points, parameterize variables:
  else {
    m_ws->factory(Form("a_muCBNom_c%d[-0.38,-1.0,1.0]",cateIndex));
    m_ws->factory(Form("b_muCBNom_c%d[-0.06,-0.1,0.1]",cateIndex));
    m_ws->factory(Form("c_muCBNom_c%d[-0.02,-0.1,0.1]",cateIndex));
    m_ws->factory(Form("a_sigmaCBNom_c%d[1.54,0.5,4.0]",cateIndex));
    m_ws->factory(Form("b_sigmaCBNom_c%d[0.90,0.1,2.0]",cateIndex));

    if (function.EqualTo("CBGA")) {
      m_ws->factory(Form("a_alphaCB_c%d[2.2,0.0,4.0]",cateIndex));
      m_ws->factory(Form("b_alphaCB_c%d[0.0,-0.1,0.1]",cateIndex));
      m_ws->factory(Form("nCB_c%d[5.0,0.1,10.0]",cateIndex));
      m_ws->factory(Form("a_sigmaGANom_c%d[5.0,0.1,20.0]",cateIndex));
      m_ws->factory(Form("b_sigmaGANom_c%d[1.0,0.1,2.0]",cateIndex));
      m_ws->factory(Form("fracCB_c%d[0.9,0.0,1.0]",cateIndex));
    }
    else if (function.EqualTo("DoubleCB")) {
      m_ws->factory(Form("a_alphaCBLo_c%d[2.42,1.0,4.0]",cateIndex));
      m_ws->factory(Form("b_alphaCBLo_c%d[-483,-1000,0]",cateIndex));
      m_ws->factory(Form("c_alphaCBLo_c%d[380,100,500]",cateIndex));
      m_ws->factory(Form("nCBLo_c%d[9.0,0.1,20.0]", cateIndex));
      m_ws->factory(Form("a_alphaCBHi_c%d[2.2,0.0,4.0]",cateIndex));
      m_ws->factory(Form("b_alphaCBHi_c%d[0.0,-0.1,0.1]",cateIndex));
      m_ws->factory(Form("nCBHi_c%d[5.0,0.1,10.0]",cateIndex));
    }
    
    // Loop over mass points, define resonance model in each:
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {

      TString currKey = getKey(currMassPoints[i_m],cateIndex);
      
      RooFormulaVar *muCBNom = new RooFormulaVar(Form("muCBNom_%s",currKey.Data()), Form("@0+@1*%f+@2*%f*%f+%f",regularizedMass(currMassPoints[i_m]),regularizedMass(currMassPoints[i_m]),regularizedMass(currMassPoints[i_m]),currMassPoints[i_m]), RooArgList(*m_ws->var(Form("a_muCBNom_c%d",cateIndex)), *m_ws->var(Form("b_muCBNom_c%d",cateIndex)), *m_ws->var(Form("c_muCBNom_c%d",cateIndex))));
      RooFormulaVar *sigmaCBNom = new RooFormulaVar(Form("sigmaCBNom_%s",currKey.Data()), Form("@0+@1*%f",regularizedMass(currMassPoints[i_m])), RooArgList(*m_ws->var(Form("a_sigmaCBNom_c%d",cateIndex)), *m_ws->var(Form("b_sigmaCBNom_c%d",cateIndex))));
      
      // Define the RooFormulaVars which control mH parameterization:
      if (function.EqualTo("CBGA")) {
	RooFormulaVar *alphaCB = new RooFormulaVar(Form("alphaCB_%s",currKey.Data()), Form("@0+@1*%f",regularizedMass(currMassPoints[i_m])), RooArgList(*m_ws->var(Form("a_alphaCB_c%d",cateIndex)), *m_ws->var(Form("b_alphaCB_c%d",cateIndex))));
	RooFormulaVar *sigmaGANom = new RooFormulaVar(Form("sigmaGANom_%s",currKey.Data()), Form("@0+@1*%f",regularizedMass(currMassPoints[i_m])), RooArgList(*m_ws->var(Form("a_sigmaGANom_c%d",cateIndex)), *m_ws->var(Form("b_sigmaGANom_c%d",cateIndex))));
	
	// Import the RooFormulaVar into the workspace:
	m_ws->import(*muCBNom);
	m_ws->import(*sigmaCBNom);
	m_ws->import(*alphaCB);
	m_ws->import(*sigmaGANom);
      }
      
      // Define the RooFormulaVars which control mH parameterization:
      else if (function.EqualTo("DoubleCB")) {
	RooFormulaVar *alphaCBLo = new RooFormulaVar(Form("alphaCBLo_%s",currKey.Data()), Form("@0+@1/(%f+@2)",regularizedMass(currMassPoints[i_m])), RooArgList(*m_ws->var(Form("a_alphaCBLo_c%d",cateIndex)), *m_ws->var(Form("b_alphaCBLo_c%d",cateIndex)), *m_ws->var(Form("c_alphaCBLo_c%d",cateIndex))));
	RooFormulaVar *alphaCBHi = new RooFormulaVar(Form("alphaCBHi_%s",currKey.Data()), Form("@0+@1*%f",regularizedMass(currMassPoints[i_m])), RooArgList(*m_ws->var(Form("a_alphaCBHi_c%d",cateIndex)), *m_ws->var(Form("b_alphaCBHi_c%d",cateIndex))));
		
	// Import the RooFormulaVar into the workspace:
	m_ws->import(*muCBNom);
	m_ws->import(*sigmaCBNom);
	m_ws->import(*alphaCBLo);
	m_ws->import(*alphaCBHi);
      }
      
      // Create the individual resonance shapes for simultaneous fitting:
      resonanceCreator(currMassPoints[i_m], cateIndex, function);
      // Add single resonance from workspace to simultaneous PDF:
      currSim->addPdf(*m_ws->pdf(Form("sigPdf_%s",currKey.Data())), currKey);
    }

    // Import simultaneous PDF into workspace:
    std::cout << "SigParam: Importing simultaneous PDF" << std::endl;
    m_ws->import(*currSim);
    
    // Get the RooDataSets:
    std::cout << "SigParam: Beginning dataset combination." << std::endl;
    std::map<std::string,RooDataSet*> currDataMap; currDataMap.clear();
    for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
      // maybe add if statement here...
      TString currKey = getKey(currMassPoints[i_m], cateIndex);
      currDataMap[((std::string)currKey)]
	= (RooDataSet*)m_ws->data(Form("data_%s",currKey.Data()));
    }
    std::cout << "SigParam: Continuing dataset combination 1." << std::endl;
    RooArgSet *args = new RooArgSet();
    args->add(*(m_ws->var("m_yy")));
    args->add(*(m_ws->var("wt")));
    
    std::cout << "SigParam: Continuing dataset combination 2." << std::endl;
    RooDataSet *obsData = new RooDataSet(Form("data_c%d",cateIndex),
					 Form("data_c%d",cateIndex), *args,
					 RooFit::Index(*m_cat),
					 RooFit::Import(currDataMap),
					 RooFit::WeightVar(*m_wt));
    m_ws->import(*obsData);
    
    // Then fit simultaneous PDF to combined dataset:
    RooFitResult *result = fitResult(cateIndex);
    // Return the fit status:
    return (result->status() == 0);
  }
}

/**
   -----------------------------------------------------------------------------
   Create the resonance for a single mass point and category. 
   @param resonanceMass - the truth mass of the resonance
   @param cateIndex - the index of the category to fit.
   @param function - the functional form of the resonance.
   @returns - true iff. all fits converge.
*/
bool SigParam::makeSingleResonance(double resonanceMass, int cateIndex,
				   TString function) {
  // Before calling resonanceCreator, need to define dependent variables.
  TString currKey = getKey(resonanceMass, cateIndex);
  if (function.EqualTo("CBGA")) {
    m_ws->factory(Form("muCBNom_%s[%f,%f,%f]",currKey.Data(),
		       resonanceMass, 0.9*resonanceMass, 1.1*resonanceMass));
    m_ws->factory(Form("sigmaCBNom_%s[%f,%f,%f]",currKey.Data(),2.0,0.01,40.0));
    m_ws->factory(Form("alphaCB_%s[%f,%f,%f]",currKey.Data(),1.5,1.0,2.5));
    if (!m_ws->var(Form("nCB_c%d",cateIndex))) {
      m_ws->factory(Form("nCB_c%d[%f,%f,%f]",cateIndex,9.0,0.01,20.0));
    }
    m_ws->factory(Form("sigmaGANom_%s[%f,%f,%f]",currKey.Data(),2.0,0.01,40.0));
    m_ws->factory(Form("fracCB_c%d[%f,%f,%f]",cateIndex,0.9,0.01,1.0));
  }
  else if (function.EqualTo("DoubleCB")) {
    m_ws->factory(Form("muCBNom_%s[%f,%f,%f]",currKey.Data(),
		       resonanceMass, 0.9*resonanceMass, 1.1*resonanceMass));
    m_ws->factory(Form("sigmaCBNom_%s[%f,%f,%f]",currKey.Data(),2.0,0.01,40.0));
    m_ws->factory(Form("alphaCBLo_%s[%f,%f,%f]",currKey.Data(),1.5,1.0,2.5));
    if (!m_ws->var(Form("nCBLo_c%d",cateIndex))) {
      m_ws->factory(Form("nCBLo_c%d[%f,%f,%f]",cateIndex,17.0,0.01,30.0));
    }
    m_ws->factory(Form("alphaCBHi_%s[%f,%f,%f]",currKey.Data(),2.2,1.00,3.0));
    if (!m_ws->var(Form("nCBHi_c%d",cateIndex))) {
      m_ws->factory(Form("nCBHi_c%d[%f,%f,%f]",cateIndex,5.2,0.01,10.0));
    }
  }
  // Adds an object called "sigPdf_%s" to the workspace, %s=getKey() value.
  resonanceCreator(resonanceMass, cateIndex, function);
  
  // Then fit single PDF to single dataset:
  RooFitResult *result = fitResult(resonanceMass,cateIndex);
  // Return the fit status:
  return (result->status() == 0);
}

/**
   -----------------------------------------------------------------------------
   Convert the resonance mass integer to a floating value.
   @param massInteger - an integer value representing the mass.
   @returns - the floating value of the mass in GeV.
*/
double SigParam::massIntToDouble(int massInteger) {
  return ((double)massInteger) / 1000.0;
}

/**
   -----------------------------------------------------------------------------
   Convert the resonance mass value to an integer representation. 
   @param resonanceMass - the value of the mass.
   @returns - the integer representation of the mass.
*/
int SigParam::massDoubleToInt(double resonanceMass) {
  return (int)(resonanceMass * 1000.0);
}

/**
   -----------------------------------------------------------------------------
   Get a list of mass points corresponding to a single category.
   @param cateIndex - the index of the category.
   @returns - a vector of mass values.
*/
std::vector<double> SigParam::massPointsForCategory(int cateIndex) {
  // Create a list of mass points in this category.
  std::vector<double> currMassPoints; currMassPoints.clear();
  for (int i_p = 0; i_p < (int)m_massCatePairs.size(); i_p++) {
    if ((m_massCatePairs[i_p]).second == cateIndex) {
      currMassPoints.push_back((m_massCatePairs[i_p]).first);
    }
  }
  std::sort(currMassPoints.begin(), currMassPoints.end());
  return currMassPoints;
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for all masses defined for one category.
   @param cateIndex - the index of the category.
   @param fileName - the location of the plot file (include format, .pdf etc.)
*/
void SigParam::plotCategoryResonances(int cateIndex, TString fileName) {
  std::cout << "SigParam: Plotting resonances in category " << cateIndex
	    << std::endl;

  // Get a list of the resonance masses:
  std::vector<double> currMassPoints = massPointsForCategory(cateIndex);
  int xMin = currMassPoints[0] - 10;
  int xMax = currMassPoints[currMassPoints.size()-1] + 10;
  int xBins = 20 + xMax - xMin;
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  RooPlot* frame = m_yy->frame(RooFit::Bins(xBins), RooFit::Range(xMin, xMax));
  frame->SetYTitle("Events/0.5 GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  
  // Loop over mass points, drawing data and PDF for each.
  for (int i_m = 0; i_m < (int)currMassPoints.size(); i_m++) {
    TString currKey = getKey(currMassPoints[i_m], cateIndex);
    double currN = (*m_ws->data(Form("data_%s",currKey.Data()))).sumEntries();
    (*m_ws->data(Form("data_%s",currKey.Data()))).plotOn(frame);
    (*m_ws->pdf(Form("sigPdf_%s",currKey.Data()))).plotOn(frame,
							  RooFit::LineColor(4));
    if (i_m == 0) frame->Draw();
    else frame->Draw("SAME");
  }
    
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.75, 0.88, Form("category %d", cateIndex));
  can->Print(fileName);
}

/**
   -----------------------------------------------------------------------------
   Plot a resonance PDF for one value of the resonance mass in one category.
   @param resonanceMass - the mass value in GeV.
   @param cateIndex - the index of the category.
   @param fileName - the location of the plot file (include format, .pdf etc.)
*/
void SigParam::plotSingleResonance(double resonanceMass, int cateIndex,
				   TString fileName) {
  std::cout << "SigParam: Plotting resonance at mass " << resonanceMass 
	    << " in category " << cateIndex << std::endl;
  
  TCanvas *can = new TCanvas("can","can",800,600);
  can->cd();
  RooPlot* frame = m_yy->frame(RooFit::Bins(40),
			      RooFit::Range(resonanceMass-10,resonanceMass+10));
  frame->SetYTitle("Events/0.5 GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  TString currKey = getKey(resonanceMass,cateIndex);
  (*m_ws->data(Form("data_%s",currKey.Data()))).plotOn(frame);
  (*m_ws->pdf(Form("sigPdf_%s",currKey.Data())))
    .plotOn(frame, RooFit::LineColor(2));
  
  frame->Draw();
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.18, 0.88, Form("category %d", cateIndex));
  can->Print(fileName);
}

/**
   -----------------------------------------------------------------------------
   Create a regularized mass variable that helps minimizers converge.
   @param resonanceMass - the mass value in GeV.
   @returns - the regularized mass mR = (m-100)/100;
*/
double SigParam::regularizedMass(double resonanceMass) {
  return ((resonanceMass - 100.0) / 100.0);
}

/**
   -----------------------------------------------------------------------------
   Create the resonance shape corresponding to a single mass point in a single
   analysis category. The shape will be stored in the workspace under the name
   "sigPdf_%s", where %s is given by the getKey() method of this class.
   @param resonanceMass - the truth mass of the resonance.
   @param cateIndex - the index of the category to fit.
   @param function - the functional form of the resonance.
*/
void SigParam::resonanceCreator(double resonanceMass, int cateIndex, 
				TString function) {
  std::cout << "SigParam: Adding resonance to workspace" << "\tresonanceMass: "
	    << resonanceMass << "\tcategory: " << cateIndex << std::endl;
  
  // Check that the dataset has been defined and is not empty.
  if (!dataExists(resonanceMass, cateIndex)) {
    std::cout << "SigParam: Cannot fit resonance, no dataset." << std::endl;
    exit(0);
  }
  
  TString currKey = getKey(resonanceMass, cateIndex);
  
  // Define the Crystal Ball + Gaussian shape:
  if (function.EqualTo("CBGA")) {
    // Cystal Ball component:
    m_ws->factory(Form("RooCBShape::pdfCB_%s(m_yy, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCB_%s, nCB_c%d)", currKey.Data(), currKey.Data(), currKey.Data(), m_listMSS.Data(), currKey.Data(), currKey.Data(), m_listMRS.Data(), currKey.Data(), cateIndex));
    // Gaussian component:
    m_ws->factory(Form("RooGaussian::pdfGA_%s(m_yy, prod::muGA_%s(muCBNom_%s%s), prod::sigmaGA_%s(sigmaGANom_%s%s))", currKey.Data(), currKey.Data(), currKey.Data(), m_listMSS.Data(), currKey.Data(), currKey.Data(), m_listMRS.Data()));
    // Crystal Ball + Gaussian:
    m_ws->factory(Form("SUM::sigPdf_%s(fracCB_c%d*pdfCB_%s,pdfGA_%s)", currKey.Data(), cateIndex, currKey.Data(), currKey.Data()));
  }
  
  // Define double-sided Crystal Ball shape:
  else if (function.EqualTo("DoubleCB")) {
    m_ws->factory(Form("HggTwoSidedCBPdf::sigPdf_%s(m_yy, prod::muCB_%s(muCBNom_%s%s), prod::sigmaCB_%s(sigmaCBNom_%s%s), alphaCBLo_%s, nCBLo_c%d, alphaCBHi_%s, nCBHi_c%d)", currKey.Data(), currKey.Data(), currKey.Data(), m_listMSS.Data(), currKey.Data(), currKey.Data(), m_listMRS.Data(), currKey.Data(), cateIndex, currKey.Data(), cateIndex));
  }
  
  std::cout << "SigParam: Resonance succesfully added." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Save the workspace containing the parameterization data to file.
   @param fileName - name of the file (full path, ending with .root extension.
*/
void SigParam::saveParameterization(TString fileName) {
  m_ws->writeToFile(fileName);
}

/**
   -----------------------------------------------------------------------------
   Make the parameters of a PDF free or fixed.
   @param pdf - the PDF containing the parameters to be freed/fixed.
   @param isConstant - true iff setting the parameters constant.
*/
void SigParam::setParamsConstant(RooAbsPdf* pdf, bool isConstant) {
  
  RooArgSet *currArgs = pdf->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    currIter->setConstant(isConstant);
  }
}

/**
   -----------------------------------------------------------------------------
   Set the function to use for parameterization of a variable as a function of
   the resonance mass.
   @param varName - the variable name.
   @param function - the name of the function for parameterizing the variable.
*/
void SigParam::setVarParameterization(TString varName, TString function) {
  std::cout << "SigParam: Using " << function << " to parameterize " << varName
	    << std::endl;
  m_funcList[varName] = function;
}
