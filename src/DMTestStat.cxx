////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMTestStat.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 30/07/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0, CL, and CLs based on an input //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMTestStat.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DMTestStat class. 
   @param newConfigFile - The analysis config file.
   @param newDMSignal - The Dark Matter signal to incorporate in the model.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
DMTestStat::DMTestStat(TString newConfigFile, TString newDMSignal,
		       TString newOptions, RooWorkspace *newWorkspace) {
  std::cout << "DMTestStat: Initializing...\n\t" << newConfigFile << "\n\t"
	    << newDMSignal << "\n\t" << newOptions << "\n\t" << std::endl;
  
  // Assign input variables: 
  m_DMSignal = newDMSignal;
  m_options = newOptions;
  
  // Start with a clean class:
  clearData();
  
  // Set the analysis configuration:
  m_config = new Config(newConfigFile);
  m_jobName = m_config->getStr("jobName");
  m_cateScheme = m_config->getStr("cateScheme");
  
  // Use Asimov data if the analysis is blind.
  m_dataForObs = (config->getBool("doBlind")) ? "asimovDataMu0" : "obsData";
  m_dataForExp = "asimovDataMu0";
  
  if (newWorkspace == NULL) {
    m_inputFile
      = new TFile(Form("%s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
		       (config->getStr("masterOutput")).Data(),
		       m_jobName.Data(), m_DMSignal.Data()), "read");
    
    if (m_inputFile->IsOpen()) {
      std::cout << "DMTestStat: Loading workspace." << std::endl;
      m_workspace = (RooWorkspace*)m_inputFile->Get("combinedWS");
    }
    else {
      std::cout << "DMTestStat: Error loading workspace." << std::endl;
      
      // Load the workspace from the nominal location.
      DMWorkspace *m_dmws = new DMWorkspace(newConfigFile, newDMSignal,
					    "FromFile");
      m_workspace = m_dmws->getCombinedWorkspace();
    }
  }
  // Use the workspace passed to the class constructor:
  else m_workspace = newWorkspace;
  
  m_mc = (ModelConfig*)m_workspace->obj("modelConfig");
  
  // Map storing all calculations:
  m_calculatedValues.clear();
  
  // Create output directories:
  m_outputDir = Form("%s/%s/TestStat/",
		     (config->getStr("masterOutput")).Data(), m_jobName.Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/CL/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/p0/", m_outputDir.Data()));
  
  // Make new or load old values:
  if (m_options.Contains("FromFile")) loadStatsFromFile();
  
  // Finished instantiating the test statistic class. 
  std::cout << "DMTestStat: Initialized Successfully!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Get the value of one of the test statistics.
   @param testStat - the test stat. name (p0, CL, CLs).
   @param observed - true iff observed, false if expected. 
   @param N - the standard deviation (-2, -1, 0, +1, +2). 
   @returns - the test statistic value.
*/
double DMTestStat::accessValue(TString testStat, bool observed, int N) {
  TString currMapKey = getKey(testStat, observed, N);
  // Check that corresponding entry exists:
  if (mapValueExists(currMapKey)) return m_calculatedValues[currMapKey];
  else return 0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the CL and CLs values using model fits.
*/
void DMTestStat::calculateNewCL() {
  std::cout << "DMTestStat: Calculating CLs" << std::endl;
  
  // Calculate observed qmu: 
  double muHatObs = 0.0;
  double nllMu1Obs = getFitNLL(m_dataForObs, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObs, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(m_dataForExp, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExp, 0.0, false, muHatExp);
  double expQMu = getQMuFromNLL(nllMu1Exp, nllMuHatExp, muHatExp, 1);
  
  // Calculate CL:
  double expCLn2 = getCLFromQMu(expQMu, 0, -2);
  double expCLn1 = getCLFromQMu(expQMu, 0, -1);
  double expCLp1 = getCLFromQMu(expQMu, 0, 1);
  double expCLp2 = getCLFromQMu(expQMu, 0, 2);
  double expCL = getCLFromQMu(expQMu, 0, 0);
  double obsCL = getCLFromQMu(obsQMu, 1, 0);
  
  // Write CL values to file:
  ofstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt", m_outputDir.Data(), 
		   m_DMSignal.Data()));
  textCL << m_DMSignal << " " << obsCL << " " << expCLn2 << " " << expCLn1
	 << " " << expCL << " " << expCLp1 << " " << expCLp2 << std::endl;
  textCL.close();
  
  // Print summary:
  std::cout << "  expected CL +2s = " << expCLp2 << std::endl;
  std::cout << "  expected CL +1s = " << expCLp1 << std::endl;
  std::cout << "  expected CL nom = " << expCL    << std::endl;
  std::cout << "  expected CL -1s = " << expCLn1 << std::endl;
  std::cout << "  expected CL -2s = " << expCLn2 << std::endl;
  std::cout << "  observed CL = " << obsCL << std::endl;
  std::cout << " " << std::endl;
  if (m_allGoodFits) std::cout << "All good fits? True" << std::endl;
  else std::cout << "All good fits? False" << std::endl;
  cout << " " << endl;
  if (obsQMu < 0) std::cout << "WARNING! obsQMu < 0 : " << obsQMu << std::endl;
  if (expQMu < 0) std::cout << "WARNING! expQMu < 0 : " << expQMu << std::endl;
  
  // save CL and CLs for later access:
  m_calculatedValues[getKey("CL",0,-2)] = expCLn2;
  m_calculatedValues[getKey("CL",0,-1)] = expCLn1;
  m_calculatedValues[getKey("CL",0,0)] = expCL;
  m_calculatedValues[getKey("CL",0,1)] = expCLp1;
  m_calculatedValues[getKey("CL",0,2)] = expCLp2;
  m_calculatedValues[getKey("CL",1,0)] = obsCL;
  
  m_calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(expCLn2);
  m_calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(expCLn1);
  m_calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(expCL);
  m_calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(expCLp1);
  m_calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(expCLp2);
  m_calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(obsCL);
}

/**
   -----------------------------------------------------------------------------
   Calculate the p0 value using model fits.
*/
void DMTestStat::calculateNewP0() {
  std::cout << "DMTestStat: calculating p0." << std::endl;
  
  // Calculate observed q0: 
  double muHatObs = 0.0;
  double nllMu0Obs = getFitNLL(m_dataForObs, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObs, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(m_dataForExp, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExp, 0.0, false, muHatExp);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  
  // Calculate p0 from q0:
  double expP0 = getP0FromQ0(expQ0);
  double obsP0 = getP0FromQ0(obsQ0);
  
  // Write p0 values to file:
  ofstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", m_outputDir.Data(), 
		   m_DMSignal.Data()));
  textP0 << m_DMSignal << " " << expP0 << " " << obsP0 << std::endl;
  textP0.close();
  
  // Print summary:
  std::cout << "\n  Expected p0 = " << expP0 << std::endl;
  std::cout << "  Observed p0 = " << obsP0 << std::endl;
  if (fitsAllConverged()) {
    std::cout << "All good fits? True\n" << std::endl;
  }
  else {
    std::cout << "All good fits? False\n" << std::endl;
  }
  
  // Save p0 for later access:
  m_calculatedValues[getKey("p0", 1, 0)] = obsP0;
  m_calculatedValues[getKey("p0", 0, 0)] = expP0;
}

/**
   -----------------------------------------------------------------------------
   Clears all data stored by the class, but does not modify the workspace.
*/
void DMTestStat::clearData() {
  m_allGoodFits = true;
  m_calculatedValues.clear();
  m_namesGlobs.clear();
  m_namesNP.clear();
  m_valuesGlobs.clear();
  m_valuesNP.clear();
  
  //m_nBins = 240;
  
  m_doSaveSnapshot = false;
  m_doPlot = false;
  m_plotDir = "";

  clearFitParamSettings();
}

/**
   -----------------------------------------------------------------------------
   Clears all specifications for parameter values during fits.
*/
void DMTestStat::clearFitParamSettings() {
  m_setParamConsts.clear();
  m_setParamNames.clear();
  m_setParamVals.clear();
}

/**
   -----------------------------------------------------------------------------
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param cateWS - the current category workspace.
   @param valMuDM - the value of the dark matter signal strength to use.
   @param valMuSM - the value of the SM Higgs signal strength to use.
   @returns - A RooDataSet with Asimov data, also imports in workspace. 
*/
RooDataSet* DMTestStat::createAsimovData(int valMuDM, int valMuSM) {
  std::cout << "DHWorkspace: Creating Asimov data, mu_DM = " << valMuDM
	    << " and mu_SM = " << valMuSM << std::endl;
  
  // Set mu_DH and mu_SH to the specified values:
  RooRealVar *poi = m_workspace->var("mu_DM");
  double initialMuDM = poi->getVal();
  double initialMuSM = m_workspace->var("mu_SM")->getVal();
  poi->setVal(valMuDM);
  poi->setConstant(true);
  m_workspace->var("mu_SM")->setVal(valMuSM);
  m_workspace->var("mu_SM")->setConstant(true);
    
  //RooDataSet *asimovData = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_workspace->pdf("combinedPdf"), *m_workspace->set("observables"));
  RooDataSet *asimovData = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_mc->GetPdf(), *m_mc->GetObservables());
  
  TString asimovName = Form("asimovDataMu%d_%s", valMuDM, m_DMSignal.Data());
  asimovData->SetNameTitle(asimovName, asimovName);
  
  m_workspace->import(*asimovData);
  m_workspace->var("mu_DM")->setVal(initialMuDM);
  m_workspace->var("mu_SM")->setVal(initialMuSM);
  m_workspace->var("mu_DM")->setConstant(false);
  m_workspace->var("mu_SM")->setConstant(true);
  std::cout << "DMWorkspace: Asimov data has " << asimovData->sumEntries() 
	    << " entries" << std::endl;
  return asimovData;
}

/**
   -----------------------------------------------------------------------------
   Create an Asimov dataset using the name.
   @param datasetName - the name of the Asimov data set.
   @returns - A RooDataSet with Asimov data, also imports in workspace.
*/
RooDataSet* DMTestStat::createAsimovData(TString datasetName) {
  if (datasetName.Contains("Mu1")) return createAsimovData(1, 1);
  else if (datasetName.Contains("Mu0")) return createAsimovData(0, 1);
  else {
    std::cout << "DMTestStat: Asimov data creation error for " << datasetName
	      << std::endl;
    return NULL;
  }
}

/**
   -----------------------------------------------------------------------------
   Create a pseudo-dataset with a given value of DM and SM signal strength.
   @param seed - the random seed for dataset generation.
   @param valMuDM - the value of the DM signal strength.
   @param valMuSM - the value of the SM signal strength.
   @returns - a pseudo-dataset.
*/
RooDataSet* DMTestStat::createPseudoData(int seed, int valMuDM, int valMuSM,
					 bool fixMu) {
  std::cout << "DMTestStat: Create pseudodata with seed = " << seed 
	    << "muDM = " << valMuDM << " muSM = " << valMuSM << std::endl;
  
  // Load the original parameters from profiling:
  m_workspace->loadSnapshot("paramsOrigin");
  
  RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  RooArgSet* originValsNP
    = (RooArgSet*)m_mc->GetNuisanceParameters()->snapshot();
  RooRealVar* firstPOI = (RooRealVar*)m_mc->GetParametersOfInterest()->first();
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
  
  map<string,RooDataSet*> toyDataMap; 
  RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  RooDataSet *dataTemp[20];
  
  // Loop over all channels:
  int index = 0;
  // Previously this was commented and similar line below was uncommented
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
  
  //numEventsPerCate.clear();
    
  // Set mu_DH and mu_SH to the specified values:
  RooRealVar *poi = m_workspace->var("mu_DM");
  double initialMuDM = poi->getVal();
  double initialMuSM = m_workspace->var("mu_SM")->getVal();
  poi->setVal(valMuDM);
  poi->setConstant(fixMu);
  m_workspace->var("mu_SM")->setVal(valMuSM);
  m_workspace->var("mu_SM")->setConstant(true);
  
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (int i_p = 0; i_p < (int)m_setParamNames.size(); i_p++) {
    std::cout << "i_p=" << i_p << ", " << m_setParamNames[i_p] << std::endl;
    m_workspace->var(m_setParamNames[i_p])->setVal(m_setParamVals[i_p]);
    m_workspace->var(m_setParamNames[i_p])->setConstant(m_setParamConsts[i_p]);
  }
  
  // Iterate over the categories:
  while ((cateType=(RooCatType*) cateIter->Next())) {
    RooAbsPdf *currPDF = combPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPDF->getObservables(observables);
    RooArgSet *currGlobs = currPDF->getObservables(globalObservables);
    RooRealVar *t = (RooRealVar*)currObs->first();
    
    //statistics::randomizeSet(currPDF, currGlobs, -1);
    //statistics::constSet(currGlobs, true);
    
    // If you want to bin the pseudo-data (speeds up calculation):
    if (m_options.Contains("Binned")) {
      currPDF->setAttribute("PleaseGenerateBinned");
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) {
	currObs->setBins(120);
      }
      dataTemp[index]
	= (RooDataSet*)currPDF->generate(*currObs, AutoBinned(true),
					 Extended(currPDF->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Construct unbinned pseudo-data by default:
    else {
      dataTemp[index] = (RooDataSet*)currPDF->generate(*currObs,Extended(true));
    }
    
    toyDataMap[(std::string)cateType->GetName()] = dataTemp[index];
    //numEventsPerCate.push_back((double)dataTemp[index]->sumEntries());
    index++;
  }
  
  // Import the new data into the workspace:
  RooDataSet* pseudoData = new RooDataSet("toyData", "toyData", *observables, 
					  RooFit::Index(*categories),
					  RooFit::Import(toyDataMap));
  
  // release nuisance parameters:
  statistics::constSet(nuisanceParameters, false);
  
  // Import into the workspace:
  m_workspace->import(*pseudoData);
  
  return pseudoData;
}

/**
   -----------------------------------------------------------------------------
   Check if all of the fits done by this class have converged.
   @returns - true iff. all of the fits have been successfully convergent.
*/
bool DMTestStat::fitsAllConverged() { 
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMu.
   @param x - the value of the test statistic.
   @returns - the value of the asymptotic test statistic distribution.
*/
double DMTestStat::functionQMu(double x) {
  // This corresponds to the "special case" of mu=mu'
  double result = TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMuTilde.
   @param x - the value of the test statistic.
   @param asimovTestStat - the test stat value on Asimov data with mu=0 but
                           fitting under mu=1 hypothesis.
   @returns - the value of the asymptotic test statistic distribution.
*/
double DMTestStat::functionQMuTilde(double x, double asimovTestStat) {
  // This corresponds to the "special case" of mu=mu'
  double result = 0.0;
  double cutoff = asimovTestStat; // asimov test stat...
  if (x == 0) result = 0.0;
  else if (x <= cutoff) {
    result = (TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x)));
  }
  else {
    result = (TMath::Exp(-1*TMath::Power((x+cutoff),2) / (8*cutoff))
	      / (2*sqrt(2*TMath::Pi()*cutoff)));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the CL value from CLs.
   @param CLs - the CLs value to convert to CL.
   @returns - the corresponding CL value.
*/
double DMTestStat::getCLFromCLs(double CLs) {
  return (1.0 - CLs);
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value from CL.
   @param CL - the CL value to convert to CLs.
   @returns - the corresponding CLs value.
*/
double DMTestStat::getCLsFromCL(double CL) {
  return (1.0 - CL);
}

/**
   -----------------------------------------------------------------------------
   Get the CL value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param observed - true of observed stat., false if expected result.
   @param N - the sigma value (-2,-1,0,1,2). Use 0 for median.
   @returns - the CLs value.
*/
double DMTestStat::getCLFromQMu(double qMu, bool observed, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, observed, N));
  return CL;
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param observed - true of observed stat., false if expected result.
   @param N - the sigma value (-2,-1,0,1,2). Use 0 for median.
   @returns - the CLs value.
*/
double DMTestStat::getCLsFromQMu(double qMu, bool observed, double N) {
  // N = 0 for exp and obs
  double pMu = getPMuFromQMu(qMu);
  double pB = getPbFromN(N);
  double CLs = pMu / (1.0 - pB);
  return CLs;
}

/**
   -----------------------------------------------------------------------------
   Get the negative-log-likelihood for a fit of a specified type to a specified
   dataset.
   @param datasetName - the name of the dataset in the workspace.
   @param muVal - the mu value to fix.
   @param fixMu - true if mu should be fixed to the specified value.
   @param &profiledMu - the profiled value of mu (passed by reference)
   @returns - the nll value.
*/
double DMTestStat::getFitNLL(TString datasetName, double muVal, bool fixMu,
			     double &profiledMu) { 
  std::cout << "DMTestStat: getFitNLL(" << datasetName << ", " << muVal
	    << ", " << fixMu << ")" << std::endl;
    
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  m_workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)m_workspace->getSnapshot("paramsOrigin");
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstpoi = (RooRealVar*)poi->first();
   RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
    
  // Look for dataset, and create if non-existent and Asimov. 
  RooAbsData *data = NULL;
  if (m_workspace->data(datasetName)) data = m_workspace->data(datasetName);
  else if (datasetName.Contains("asimovData")) {
    data = createAsimovData(datasetName);
  }
  else {
    std::cout << "DMTestStat: Error! Requested data not available: " 
	      << datasetName << std::endl;
    exit(0);
  }
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  firstpoi->setVal(muVal);
  firstpoi->setConstant(fixMu);
  
  // Check if parameter settings have been specified during fit:
  for (int i_p = 0; i_p < (int)m_setParamNames.size(); i_p++) {
    std::cout << "i_p=" << i_p << ", " << m_setParamNames[i_p] << std::endl;
    m_workspace->var(m_setParamNames[i_p])->setVal(m_setParamVals[i_p]);
    m_workspace->var(m_setParamNames[i_p])->setConstant(m_setParamConsts[i_p]);
  }
  
  // Iterate over SM mu values and fix all to 1:
  RooArgSet* muConstants = (RooArgSet*)m_workspace->set("muSMConstants");
  TIterator *iterMuConst = muConstants->createIterator();
  RooRealVar *currMuConst = NULL;
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    std::cout << "DMTestStat: Setting " << currMuConst->GetName()
	      << " constant." << std::endl;
    currMuConst->setVal(1.0);
    currMuConst->setConstant(true);
  }
  
  int status = 0; 
  RooNLLVar* varNLL = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize(status, varNLL, "", NULL, false);
  if (status != 0) m_allGoodFits = false;
  
  // Save a snapshot if requested:
  if (m_doSaveSnapshot) {
    TString muDMValue = fixMu ? (Form("%d",(int)muVal)) : "Free";
    m_workspace->saveSnapshot(Form("paramsProfileMu%s", muDMValue.Data()),
			      *poiAndNuis);
  }
  
  // Plot the fit result if the user has set an output directory for plots:
  if (m_doPlot) {
    if (fixMu && ((int)muVal) == 1) plotFits("Mu1", datasetName);
    else if (fixMu && ((int)muVal) == 0) plotFits("Mu0", datasetName);
    else plotFits("MuFree", datasetName);
  }
  
  // Save the NLL and mu from profiling:
  profiledMu = firstpoi->getVal();
  double nllValue = varNLL->getVal();
  delete varNLL;
  
  // Save names and values of nuisance parameters:
  m_namesNP.clear();
  m_valuesNP.clear();
  TIterator *iterNuis = nuisanceParameters->createIterator();
  RooRealVar *currNuis = NULL;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    m_namesNP.push_back((std::string)currNuis->GetName());
    m_valuesNP.push_back(currNuis->getVal());
  }
  
  // Save names and values of global observables:
  m_namesGlobs.clear();
  m_valuesGlobs.clear();
  TIterator *iterGlobs = globalObservables->createIterator();
  RooRealVar *currGlob = NULL;
  while ((currGlob = (RooRealVar*)iterGlobs->Next())) {
    m_namesGlobs.push_back((std::string)currGlob->GetName());
    m_valuesGlobs.push_back(currGlob->getVal());
  }
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  return nllValue;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of global observable names from the most recent fit.
*/
std::vector<std::string> DMTestStat::getGlobsNames() {
  return m_namesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of global observable values from the most recent fit.
*/
std::vector<double> DMTestStat::getGlobsValues() {
  return m_valuesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get the key for the value map.
   @param testStat - the test statistic.
   @param observed - true iff. observed.
   @param N - the sigma value (-2,-1,0,1,2).
*/
TString DMTestStat::getKey(TString testStat, bool observed, int N) {
  TString currKey = testStat;
  
  if (observed) currKey += "_obs";
  else currKey += "_exp";
  
  if (N > 0) currKey += Form("_p%d",N);
  else if (N < 0) currKey += Form("_n%d",N);
  
  return currKey;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of nuisance parameter names from the most recent fit.
*/
std::vector<std::string> DMTestStat::getNPNames() {
  return m_namesNP;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of nuisance parameter values from the most recent fit.
*/
std::vector<double> DMTestStat::getNPValues() {
  return m_valuesNP;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of p0 based on the test statistic q0.
   @param q0 - the test statistic q0.
   @returns - the value of p0.
*/
double DMTestStat::getP0FromQ0(double q0) {
  double p0 = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0)));
  return p0;
}

/**
   -----------------------------------------------------------------------------
   Calculate pB based on the standard deviation.
   @param N - the standard deviation.
   @returns - the value of pB.
*/
double DMTestStat::getPbFromN(double N) {
  double pB = 1 - ROOT::Math::gaussian_cdf(N);
  return pB;
}

/**
   -----------------------------------------------------------------------------
   Calculate the pB value based on qMu.
   @param qMu - the test statistic qMu.
   @param sigma - the sigma value...
   @param mu - the mu value... 
   @returns - the value of pB.
*/
double DMTestStat::getPbFromQMu(double qMu, double sigma, double mu) {
  double pB = 1 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMu));
  return pB;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of pMu.
   @param qMu - the test statistic qMu.
   @returns - the value of pMu.
*/
double DMTestStat::getPMuFromQMu(double qMu) {
  double pMu = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu)));
  return pMu;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic q0 based on the nll.
   @param nllMu0 - nll of a fit with signal strength 0;
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @returns - the value of q0.
*/
double DMTestStat::getQ0FromNLL(double nllMu0, double nllMuHat, double muHat) {
  double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
  return q0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMu based on the nll.
   @param nllMu - nll of a fit with signal strength mu.
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @param muTest - tested value of signal strength.
   @returns - the value of qMu.
*/
double DMTestStat::getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
				 double muTest) {
  double qMu = 0;
  if (muHat < muTest) qMu = 2 * (nllMu - nllMuHat);
  else qMu = 0.0;
  return qMu;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMuTilde based on the nll.
   @param nllMu - nll of a fit with signal strength mu.
   @param nllMu0 - nll of a fit with signal strength 0.
   @param nllMuHat - nll of a fit with profiled signal strength.
   @param muHat - profiled signal strength.
   @param muTest - tested value of signal strength.
   @returns - the value of qMuTilde.
*/
double DMTestStat::getQMuTildeFromNLL(double nllMu, double nllMu0,
				      double nllMuHat, double muHat,
				      double muTest) {
  double qMuTilde = 0;
  if (muHat <= 0) qMuTilde = 2 * (nllMu - nllMu0);
  else if (muHat > 0 && muHat <= muTest) qMuTilde = 2 * (nllMu - nllMuHat);
  else if (muHat > muTest) qMuTilde = 0;
  return qMuTilde;
}

/**
   -----------------------------------------------------------------------------
   Load the statistics files (p0 and CL) that were previously generated. If none
   are found, then create from scratch automatically.
*/
void DMTestStat::loadStatsFromFile() {
  
  // Load input p0 file:
  ifstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", m_outputDir.Data(),
		   m_DMSignal.Data()));
  
  // Load input CL file:
  ifstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt", m_outputDir.Data(), 
		   m_DMSignal.Data()));
  
  // If the input files don't exist, create from scratch:
  if (!textCL || !textP0) {
    calculateNewCL();
    calculateNewP0();
    return;
  }
  
  // Read p0 values:
  TString inName; double inExpP0; double inObsP0; 
  while (!textP0.eof()) {
    textP0 >> inName >> inExpP0 >> inObsP0;
    textP0.close();
  }
  m_calculatedValues[getKey("p0",1,0)] = inObsP0;
  m_calculatedValues[getKey("p0",0,0)] = inExpP0;
  
  // Read CL values:
  double inObsCL, inExpCLn2, inExpCLn1, inExpCL, inExpCLp1, inExpCLp2;
  while (!textCL.eof()) {
    textCL >> inName >> inObsCL >> inExpCLn2 >> inExpCLn1 >> inExpCL
	   >> inExpCLp1 >> inExpCLp2;
  }
  textCL.close();
  
  // save CL and CLs for later access:
  m_calculatedValues[getKey("CL", 0, -2)] = inExpCLn2;
  m_calculatedValues[getKey("CL", 0, -1)] = inExpCLn1;
  m_calculatedValues[getKey("CL", 0, 0)] = inExpCL;
  m_calculatedValues[getKey("CL", 0, 1)] = inExpCLp1;
  m_calculatedValues[getKey("CL", 0, 2)] = inExpCLp2;
  m_calculatedValues[getKey("CL", 1, 0)] = inObsCL;
  
  m_calculatedValues[getKey("CLs", 0, -2)] = getCLsFromCL(inExpCLn2);
  m_calculatedValues[getKey("CLs", 0, -1)] = getCLsFromCL(inExpCLn1);
  m_calculatedValues[getKey("CLs", 0, 0)] = getCLsFromCL(inExpCL);
  m_calculatedValues[getKey("CLs", 0, 1)] = getCLsFromCL(inExpCLp1);
  m_calculatedValues[getKey("CLs", 0, 2)] = getCLsFromCL(inExpCLp2);
  m_calculatedValues[getKey("CLs", 1, 0)] = getCLsFromCL(inObsCL);
}

/**
   -----------------------------------------------------------------------------
   Plot the fits produced by the specified model.
   @param fitType - the type of fit.
   @param datasetName - the name of the profiled dataset.
*/
void DMTestStat::plotFits(TString fitType, TString datasetName) {
  
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  
  // loop over categories:
  for (int i_c = 0; i_c < config->getInt("nCategories"); i_c++) {
    TString currCateName = Form("%s_%d", m_cateScheme.Data(), i_c);
    RooPlot* frame =  (*m_workspace->var("m_yy_"+currCateName)).frame(50);
    m_workspace->data(Form("%s_%s", datasetName.Data(), currCateName.Data()))
      ->plotOn(frame);
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfDM_"+currCateName))), 
	      LineColor(6));
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfSM_"+currCateName))), 
	      LineColor(3));
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame, Components((*m_workspace->pdf("bkgPdf_"+currCateName))), 
	      LineColor(4));
    (*m_workspace->pdf("model_"+currCateName)).plotOn(frame, LineColor(2));
    
    //double chi2 = frame->chiSquare();
    frame->SetYTitle("Events / GeV");
    frame->SetXTitle("M_{#gamma#gamma} [GeV]");
    frame->Draw();
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.2, 0.81, Form("Category %d", i_c));
    text.DrawLatex(0.2, 0.87, Form("Signal %s", m_DMSignal.Data()));
    TH1F *histDM = new TH1F("histDM", "histDM", 1, 0, 1);
    TH1F *histSM = new TH1F("histSM", "histSM", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDM->SetLineColor(6);
    histSM->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.61, 0.63, 0.89, 0.77);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDM, "Dark matter", "l");
    leg.AddEntry(histSM, "SM Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sig. + bkg.", "l");
    leg.Draw("SAME");
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_DMSignal.Data(), fitType.Data(), currCateName.Data()));
  }
  delete can;
}

/*
   -----------------------------------------------------------------------------
   Plot the fits produced by the specified model.
   @param fitType - the type of fit.
   @param datasetName - the name of the dataset to plot.
void DMTestStat::plotFits(TString fitType, TString datasetName) {
  std::cout << "DMTestStat: Plot final fits for " << fitType << std::endl;
  std::cout << "DEBUG1" << std::endl;
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  // loop over categories:
  int nCategories = DMAnalysis::getNumCategories(m_cateScheme, m_anaType);
  for (int i_c = 0; i_c < nCategories; i_c++) {
    can->Clear();
    std::cout << "DEBUG2" << std::endl;
    TString currCateName
      = DMAnalysis::cateIndexToName(m_cateScheme, m_anaType, i_c);
    TString obsName = m_anaType.EqualTo("NonRes") ? 
      Form("m_yy_%s", currCateName.Data()) : 
      Form("m_bbyy_%s", currCateName.Data());
    RooPlot* frame =  (*m_workspace->var(obsName)).frame(50);
    TString cutName = Form("categories_%s==categories_%s::%s", m_anaType.Data(),
			   m_anaType.Data(), currCateName.Data());
    std::cout << "DEBUG3" << std::endl;
    RooDataSet *currData
      =(RooDataSet*)(m_workspace->data(Form("%s", datasetName.Data())));
    std::cout << "DEBUG3.1" << std::endl;
    //currData->plotOn(frame, RooFit::Cut(cutName));
    
    // Everything below here was commented out for debugging.
    //RooCategory *categories
    //=(RooCategory*)(m_workspace->var(Form("categories_%s",m_anaType.Data())));
    m_workspace->Print("v");
    std::cout << "DEBUG3.2" << std::endl;
    RooArgSet tempSet = m_workspace->allCats();
    tempSet.Print("v");
    std::cout << "DEBUG3.3" << std::endl;
    RooCategory *categories
      =(RooCategory*)(m_workspace->cat(Form("categories_%s",m_anaType.Data())));

    std::cout << "DEBUG4" << std::endl;
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfDH_"+currCateName))),
	      LineColor(6));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("sigPdfSH_"+currCateName))),
	      LineColor(3));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, Components((*m_workspace->pdf("bkgPdf_"+currCateName))), 
	      LineColor(4));
    (*m_workspace->pdf(Form("model_%s",m_anaType.Data())))
      .plotOn(frame, LineColor(2));
    
    std::cout << "DEBUG5" << std::endl;
    //double chi2 = frame->chiSquare();
    frame->SetYTitle("Events / GeV");
    frame->SetXTitle("M_{#gamma#gamma} [GeV]");
    frame->Draw();
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.2, 0.81, currCateName);
    TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
    TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDH->SetLineColor(6);
    histSH->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.61, 0.63, 0.89, 0.77);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDH, "Di-Higgs", "l");
    leg.AddEntry(histSH, "Single Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sum", "l");
    leg.Draw("SAME");
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_DHSignal.Data(), fitType.Data(), currCateName.Data()));
    
  }
  delete can;
}
*/

/**
   -----------------------------------------------------------------------------
   Check whether the specified value has been stored in the value map.
   @param mapKey - the key for the map of values.
   @returns - true iff the categorization has been defined. 
*/
bool DMTestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool nonExistent = (m_calculatedValues.find(mapKey) ==
		      m_calculatedValues.end());
  if (nonExistent) {
    std::cout << "DMTestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Choose whether or not to save snapshots from profiling data.
   @param doSaveSnapshot - true iff you want to save snapshots in future fits.
*/
void DMTestStat::saveSnapshots(bool doSaveSnapshot) {
  m_doSaveSnapshot = doSaveSnapshot;
}

/**
   -----------------------------------------------------------------------------
   Set an output directory and enable plotting.
   @param directory - the output directory path.
*/
void DMTestStat::setPlotDirectory(TString directory) {
  m_plotDir = directory;
  m_doPlot = true;
}

/**
   -----------------------------------------------------------------------------
   Set the named parameter to a certain value and either fix or free it.
   @param paramName - The name of the fit parameter.
   @param paramVal - The new value of the fit parameter.
   @param doSetConstant - True iff the parameter should be set constant. 
*/
void DMTestStat::setParams(TString paramName, double paramVal,
			   bool doSetConstant) {
  m_setParamConsts.push_back(doSetConstant);
  m_setParamNames.push_back(paramName);
  m_setParamVals.push_back(paramVal);
}
