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
  m_dataForObsQ0 = (m_config->getBool("doBlind")) ? "asimovDataMu1" : "obsData";
  m_dataForExpQ0 = "asimovDataMu1";
  m_dataForObsQMu =(m_config->getBool("doBlind")) ? "asimovDataMu0" : "obsData";
  m_dataForExpQMu = "asimovDataMu0";
  
  if (!newWorkspace) {
    m_inputFile
      = new TFile(Form("%s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
		       (m_config->getStr("masterOutput")).Data(),
		       m_jobName.Data(), m_DMSignal.Data()), "read");
    
    if (m_inputFile->IsOpen()) {
      std::cout << "DMTestStat: Loading workspace." << std::endl;
      m_workspace = (RooWorkspace*)m_inputFile->Get("combinedWS");
    }
    else {
      std::cout << "DMTestStat: ERROR loading workspace. "
		<< "Trying to load from DMWorkspace..." << std::endl;
      
      // Load the workspace from the nominal location.
      DMWorkspace *m_dmws = new DMWorkspace(newConfigFile, newDMSignal,
					    "FromFile");
      m_workspace = m_dmws->getCombinedWorkspace();
    }
    
    if (!m_workspace) {
      std::cout << "DMTestStat: ERROR! Workspace was not found..." << std::endl;
      exit(0);
    }
    else {
      std::cout << "DMTestStat: Printing the loaded workspace." << std::endl;
      m_workspace->Print("v");
    }
  }
  // Use the workspace passed to the class constructor:
  else m_workspace = newWorkspace;
  
  // Then load the ModelConfig if it exists, exit otherwise:
  if (m_workspace->obj("modelConfig")) {
    m_mc = (ModelConfig*)m_workspace->obj("modelConfig");
  }
  else {
    std::cout << "DMTestStat: ERROR! ModelConfig was not found." << std::endl;
    exit(0);
  }
  
  // Map storing all calculations:
  m_calculatedValues.clear();
  
  // Create output directories:
  m_outputDir = Form("%s/%s/DMTestStat/", 
		     (m_config->getStr("masterOutput")).Data(), 
		     m_jobName.Data());
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
  double nllMu1Obs = getFitNLL(m_dataForObsQMu, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQMu, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(m_dataForExpQMu, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQMu, 0.0, false, muHatExp);
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
  double nllMu0Obs = getFitNLL(m_dataForObsQ0, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQ0, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(m_dataForExpQ0, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQ0, 0.0, false, muHatExp);
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

void DMTestStat::createAsimovData(int valMuDM, int valMuSM) {
  std::cout << "DHWorkspace: Creating Asimov data, mu_DM = " << valMuDM
	    << " and mu_SM = " << valMuSM << std::endl;
  
  // Set mu_DH and mu_SH to the specified values:
  double initialMuDM = m_workspace->var("mu_DM")->getVal();
  double initialMuSM = m_workspace->var("mu_SM")->getVal();
  bool initialMuDMConst = m_workspace->var("mu_DM")->isConstant();
  bool initialMuSMConst = m_workspace->var("mu_SM")->isConstant();
  m_workspace->var("mu_DM")->setVal(valMuDM);
  m_workspace->var("mu_DM")->setConstant(true);
  m_workspace->var("mu_SM")->setVal(valMuSM);
  m_workspace->var("mu_SM")->setConstant(true);
  
  // Dataset map for combined dataset creation:
  RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  map<string,RooDataSet*> asimovDataMap; asimovDataMap.clear();
  RooArgSet *allObs = new RooArgSet();
  
  // Loop over categories, create Asimov dataset in each.
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    TString currCateName = Form("%s_%d", m_cateScheme.Data(), i_c);
    std::cout << "Check0" << std::endl;
    m_workspace->pdf(Form("model_%s",currCateName.Data()))->Print("v");
    std::cout << "Check0.1" << std::endl;
    m_workspace->set("observables")->Print("v");
    std::cout << "Check0.2" << std::endl;
    //RooDataSet *data = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_workspace->pdf(Form("model_%s",currCateName.Data())), *m_workspace->set(Form("observables_%s",currCateName.Data())));
    RooDataSet *data = (RooDataSet*)AsymptoticCalculator::GenerateAsimovData(*m_workspace->pdf(Form("model_%s",currCateName.Data())), *m_workspace->pdf(Form("model_%s",currCateName.Data()))->getObservables(*m_workspace->set("observables")));
    std::cout << "Check1" << std::endl;
    data->SetNameTitle(Form("asimovDataMu%d_%s",valMuDM,currCateName.Data()),
		       Form("asimovDataMu%d_%s",valMuDM,currCateName.Data()));
    m_workspace->import(*data);
    asimovDataMap[(std::string)currCateName] = data;
    //allObs->add(*m_workspace->set(Form("observables_%s",currCateName.Data())));
  }
  
  // Import the new data into the workspace:
  TString asimovName = Form("asimovDataMu%d", valMuDM);
  RooDataSet* asimovData
    = new RooDataSet(asimovName, asimovName, 
		     *m_workspace->set(Form("observables")),
		     RooFit::Index(*categories), RooFit::Import(asimovDataMap));
  
  m_workspace->import(*asimovData);
  m_workspace->var("mu_DM")->setVal(initialMuDM);
  m_workspace->var("mu_SM")->setVal(initialMuSM);
  m_workspace->var("mu_DM")->setConstant(initialMuDMConst);
  m_workspace->var("mu_SM")->setConstant(initialMuSMConst);
  std::cout << "DMTestStat: Asimov data has " << asimovData->sumEntries() 
	    << " entries" << std::endl;
}


   -----------------------------------------------------------------------------
   Create an Asimov dataset using the name.
   @param datasetName - the name of the Asimov data set.
   @returns - A RooDataSet with Asimov data, also imports in workspace.

void DMTestStat::createAsimovData(TString datasetName) {
  if (datasetName.Contains("Mu1")) createAsimovData(1, 1);
  else if (datasetName.Contains("Mu0")) createAsimovData(0, 1);
  else {
    std::cout << "DMTestStat: Asimov data creation error for " << datasetName
	      << std::endl;
    exit(0);
  }
}
*/

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
  while ((cateType = (RooCatType*)cateIter->Next())) {
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
    
  // Check that dataset exists:
  if (!m_workspace->data(datasetName)) {
    std::cout << "DMTestStat: Error! Requested data not available: " 
	      << datasetName << std::endl;
    exit(0);
  }
  // Check PDF also exists:
  else if (!combPdf) {
    std::cout << "DMTestStat: ERROR! Requested PDF not found..." << std::endl;
    exit(0);
  }
  
  // Free nuisance parameters before fit:
  statistics::constSet(nuisanceParameters, false, origValNP);
  // Fix the global observables before the fit:
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
   
  // The actual fit command:
  RooNLLVar* varNLL = (RooNLLVar*)combPdf
    ->createNLL(*m_workspace->data(datasetName), Constrain(*nuisanceParameters),
  		Extended(combPdf->canBeExtended()));
    
  RooFitResult *fitResult = statistics::minimize(varNLL, "", NULL, true);
  if (!fitResult || fitResult->status() != 0) m_allGoodFits = false;
  
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
   Plot the fits produced by the specified model.
   @param fitType - the type of fit.
   @param datasetName - the name of the profiled dataset.
*/
void DMTestStat::plotFits(TString fitType, TString datasetName) {
  std::cout << "DMTestStat: Plot fit " << fitType << ", " << datasetName
	    << std::endl;
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  TPad *pad1 = new TPad( "pad1", "pad1", 0.00, 0.33, 1.00, 1.00 );
  TPad *pad2 = new TPad( "pad2", "pad2", 0.00, 0.00, 1.00, 0.33 );
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->SetBorderMode(0);
  
  can->cd();
  pad1->Draw();
  pad2->Draw();
  
  int xBins = 50; 
  int xMin = m_config->getNum("DMMyyRangeLo");
  int xMax = m_config->getNum("DMMyyRangeHi");
  // Loop over categories:
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    can->cd();
    //pad1->Clear();
    pad1->cd();
    TString currCateName = Form("%s_%d", m_cateScheme.Data(), i_c);
    RooPlot* frame =  (*m_workspace->var("m_yy_"+currCateName)).frame(xBins);
    (*m_workspace->data(Form("%s_%s", datasetName.Data(), currCateName.Data())))
      .plotOn(frame);
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame,Components((*m_workspace->pdf("sigPdfDM_"+currCateName))), 
	      LineColor(kRed+1), LineStyle(2));
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame,Components((*m_workspace->pdf("sigPdfSM_"+currCateName))), 
	      LineColor(kGreen+2));
    (*m_workspace->pdf("model_"+currCateName))
      .plotOn(frame,Components((*m_workspace->pdf("bkgPdf_"+currCateName))), 
	      LineColor(kMagenta+1));
    (*m_workspace->pdf("model_"+currCateName)).plotOn(frame,LineColor(kBlue+1));
    frame->GetYaxis()->SetRangeUser(0.00001,frame->GetYaxis()->GetXmax());
    //frame->GetYaxis()->SetRangeUser(0.00001,3.0);
    frame->SetYTitle("Events / GeV");
    frame->SetXTitle("M_{#gamma#gamma} [GeV]");
    frame->Draw();
  
    // Also draw the ATLAS text:
    TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
    l.SetTextFont(72); l.SetTextSize(0.05); l.DrawLatex(0.5,0.88,"ATLAS");
    l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.62,0.88,"Internal");
    l.DrawLatex(0.5, 0.83, Form("#scale[0.8]{#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}}",(m_config->getNum("analysisLuminosity")/1000.0)));
    
    TLatex text; text.SetNDC(); text.SetTextColor(1); text.SetTextFont(42);
    text.SetTextSize(0.05);
    text.DrawLatex(0.5, 0.78, Form("Category %d", i_c));
    TH1F *histDM = new TH1F("histDM", "histDM", 1, 0, 1);
    TH1F *histSM = new TH1F("histSM", "histSM", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    TH1F *histData = new TH1F("histData", "histData", 1, 0, 1);
    histDM->SetLineColor(kRed+1);
    histDM->SetLineStyle(2);
    histSM->SetLineColor(kGreen+2);
    histBkg->SetLineColor(kMagenta+1);
    histSig->SetLineColor(kBlue+1);
    histData->SetLineColor(kBlack);
    histData->SetMarkerColor(kBlack);
    TLegend leg(0.5, 0.55, 0.89, 0.75);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDM, 
		 DMAnalysis::getPrintSampleName(m_config,m_DMSignal),"l");
    leg.AddEntry(histSM, 
		 DMAnalysis::getPrintSampleName(m_config,"SMHiggs"), "l");
    leg.AddEntry(histBkg, "Non-resonant background", "l");
    leg.AddEntry(histSig, "Signal + background", "l");
    leg.AddEntry(histData, "Data", "lep");
    leg.Draw("SAME");
    
    // Division plot:
    pad2->cd();
    TGraphErrors* subData = plotDivision((m_workspace->data(Form("%s_%s", datasetName.Data(), currCateName.Data()))), (m_workspace->pdf("model_"+currCateName)), currCateName, xMin, xMax, xBins);
    subData->GetYaxis()->SetTitle("Data / Fit");
    subData->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    subData->GetXaxis()->SetTitleOffset(0.95);
    subData->GetYaxis()->SetTitleOffset(0.7);
    subData->GetXaxis()->SetTitleSize(0.1);
    subData->GetYaxis()->SetTitleSize(0.1);
    subData->GetXaxis()->SetLabelSize(0.1);
    subData->GetYaxis()->SetLabelSize(0.1);
    subData->GetYaxis()->SetNdivisions(4);
    subData->SetMarkerColor(1);
    subData->GetXaxis()->SetRangeUser(xMin, xMax);
    //subData->GetYaxis()->SetRangeUser(0.000001, xMax);
    subData->Draw("AEP");
    TLine *line = new TLine();
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
    line->DrawLine(xMin, 1.0, xMax, 1.0);
    subData->Draw("EPSAME");
    
    // Print the canvas:
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_DMSignal.Data(), fitType.Data(), currCateName.Data()));
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot (or subtraction plot, for the moment...)
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param cateName - The name of the category.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @returns - A TGraphErrors to plot.
*/
TGraphErrors* DMTestStat::plotDivision(RooAbsData *data, RooAbsPdf *pdf, 
				       TString cateName, double xMin,
				       double xMax, double xBins) {
  RooRealVar *m_yy = (m_workspace->var("m_yy_"+cateName));
  double minOrigin = m_yy->getMin();
  double maxOrigin = m_yy->getMax();
  double nEvents = data->sumEntries();
  
  m_yy->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *m_yy,
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy),
				       RooFit::NormSet(*m_yy), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  int pointIndex = 0; int pointIndexNonZero = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    m_yy->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*m_yy), 
					 RooFit::NormSet(*m_yy), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = m_yy->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight / currPdfWeight;
    result->SetPoint(pointIndex, currMass, currWeight);
    
    double currError = originHist->GetBinError(pointIndex+1) / currPdfWeight;
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  m_yy->setMin(minOrigin);
  m_yy->setMax(maxOrigin);
  return result;
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
