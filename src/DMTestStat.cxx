////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMTestStat.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 19/04/2015                                                          //
//                                                                            //
//  This class allows the user to calculate p0, CL, and CLs based on an input //
//  workspace.                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMTestStat.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DMTestStat class. 
   @param newJobName - The name of the job
   @param newDMSignal - The Dark Matter signal to incorporate in the model.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
DMTestStat::DMTestStat(TString newJobName, TString newDMSignal,
		       TString newCateScheme, TString newOptions,
		       RooWorkspace *newWorkspace) {
  std::cout << "DMTestStat: Initializing...\n\t" << newJobName << "\n\t"
	    << newDMSignal << "\n\t" << newCateScheme << "\n\t" 
	    << newOptions << "\n\t" << std::endl;
  
  // Assign input variables: 
  jobName = newJobName;
  DMSignal = newDMSignal;
  options = newOptions;
  allGoodFits = true;
  
  // Use Asimov data if the analysis is blind.
  dataForObs = (DMAnalysis::doBlind) ? "asimovDataMu0" : "obsData";
  dataForExp = "asimovDataMu0";
  
  TFile inputFile(Form("%s/%s/DMWorkspace/rootfiles/workspaceDM_%s.root",
		       DMAnalysis::masterOutput.Data(), jobName.Data(),
		       DMSignal.Data()), "read");
  
  if (newWorkspace == NULL) {
    if (inputFile.IsOpen()) {
      std::cout << "DMTestStat: Loading workspace." << std::endl;
      workspace = (RooWorkspace*)inputFile.Get("combinedWS");
    }
    else {
      std::cout << "DMTestStat: Error loading workspace." << std::endl;
      
      // Load the workspace from the nominal location.
      dmw = new DMWorkspace(newJobName, newDMSignal, newCateScheme, "FromFile");
      workspace = dmw->getCombinedWorkspace();
    }
  }
  // Use the workspace passed to the class constructor:
  else {
    workspace = newWorkspace;
  }
  
  mc = (ModelConfig*)workspace->obj("modelConfig");
  
  // Map storing all calculations:
  calculatedValues.clear();
  
  // Create output directories:
  outputDir = Form("%s/%s/TestStat/", DMAnalysis::masterOutput.Data(), 
		   jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/CL/",outputDir.Data()));
  system(Form("mkdir -vp %s/p0/",outputDir.Data()));

  // Make new or load old values:
  if (options.Contains("FromFile")) {
    loadStatsFromFile();
  }
  //else {
  //  calculateNewCL();
  //  calculateNewP0();
  //}
  
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
  if (mapValueExists(currMapKey)) return calculatedValues[currMapKey];
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
  double nllMu1Obs = getFitNLL(dataForObs, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(dataForObs, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(dataForExp, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(dataForExp, 0.0, false, muHatExp);
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
  textCL.open(Form("%s/CL/CL_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  textCL << DMSignal << " " << obsCL << " " << expCLn2 << " " << expCLn1
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
  if (allGoodFits) std::cout << "All good fits? True" << std::endl;
  else std::cout << "All good fits? False" << std::endl;
  cout << " " << endl;
  if (obsQMu < 0) std::cout << "WARNING! obsQMu < 0 : " << obsQMu << std::endl;
  if (expQMu < 0) std::cout << "WARNING! expQMu < 0 : " << expQMu << std::endl;
  
  // save CL and CLs for later access:
  calculatedValues[getKey("CL",0,-2)] = expCLn2;
  calculatedValues[getKey("CL",0,-1)] = expCLn1;
  calculatedValues[getKey("CL",0,0)] = expCL;
  calculatedValues[getKey("CL",0,1)] = expCLp1;
  calculatedValues[getKey("CL",0,2)] = expCLp2;
  calculatedValues[getKey("CL",1,0)] = obsCL;
  
  calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(expCLn2);
  calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(expCLn1);
  calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(expCL);
  calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(expCLp1);
  calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(expCLp2);
  calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(obsCL);
}

/**
   -----------------------------------------------------------------------------
   Calculate the p0 value using model fits.
*/
void DMTestStat::calculateNewP0() {
  std::cout << "DMTestStat: calculating p0." << std::endl;
  
  // Calculate observed q0: 
  double muHatObs = 0.0;
  double nllMu0Obs = getFitNLL(dataForObs, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(dataForObs, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(dataForExp, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(dataForExp, 0.0, false, muHatExp);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  
  // Calculate p0 from q0:
  double expP0 = getP0FromQ0(expQ0);
  double obsP0 = getP0FromQ0(obsQ0);
  
  // Write p0 values to file:
  ofstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  textP0 << DMSignal << " " << expP0 << " " << obsP0 << std::endl;
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
  calculatedValues[getKey("p0",1,0)] = obsP0;
  calculatedValues[getKey("p0",0,0)] = expP0;
}

/**
   -----------------------------------------------------------------------------
   Clears all data stored by the class, but does not modify the workspace.
*/
void DMTestStat::clearData() {
  allGoodFits = true;
  calculatedValues.clear();
  namesGlobs.clear();
  namesNP.clear();
  valuesGlobs.clear();
  valuesNP.clear();
}

/**
   -----------------------------------------------------------------------------
   Check if all of the fits done by this class have converged.
   @returns - true iff. all of the fits have been successfully convergent.
*/
bool DMTestStat::fitsAllConverged() { 
  return allGoodFits;
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
    
  RooAbsPdf* combPdf = mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)workspace->getSnapshot("paramsOrigin");
  RooArgSet* poi = (RooArgSet*)mc->GetParametersOfInterest();
  RooRealVar* firstpoi = (RooRealVar*)poi->first();
  RooAbsData *data = workspace->data(datasetName);
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  firstpoi->setVal(muVal);
  firstpoi->setConstant(fixMu);
  
  // Iterate over SM mu values and fix all to 1:
  RooArgSet* muConstants = (RooArgSet*)workspace->set("muSMConstants");
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
  if (status != 0) {
    allGoodFits = false;
  }
  //profiledMu = fixMu ? muVal : firstpoi->getVal();
  profiledMu = firstpoi->getVal();
  double nllValue = varNLL->getVal();
  delete varNLL;
  
  // Save names and values of nuisance parameters:
  namesNP.clear();
  valuesNP.clear();
  TIterator *iterNuis = nuisanceParameters->createIterator();
  RooRealVar *currNuis = NULL;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    namesNP.push_back((std::string)currNuis->GetName());
    valuesNP.push_back(currNuis->getVal());
  }
  
  // Save names and values of global observables:
  namesGlobs.clear();
  valuesGlobs.clear();
  TIterator *iterGlobs = globalObservables->createIterator();
  RooRealVar *currGlob = NULL;
  while ((currGlob = (RooRealVar*)iterGlobs->Next())) {
    namesGlobs.push_back((std::string)currGlob->GetName());
    valuesGlobs.push_back(currGlob->getVal());
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
  return namesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of global observable values from the most recent fit.
*/
std::vector<double> DMTestStat::getGlobsValues() {
  return valuesGlobs;
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
  return namesNP;
}

/**
   -----------------------------------------------------------------------------
   Get a vector of nuisance parameter values from the most recent fit.
*/
std::vector<double> DMTestStat::getNPValues() {
  return valuesNP;
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
  textP0.open(Form("%s/p0/p0_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  
  // Load input CL file:
  ifstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  
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
  calculatedValues[getKey("p0",1,0)] = inObsP0;
  calculatedValues[getKey("p0",0,0)] = inExpP0;
  
  // Read CL values:
  double inObsCL, inExpCLn2, inExpCLn1, inExpCL, inExpCLp1, inExpCLp2;
  while (!textCL.eof()) {
    textCL >> inName >> inObsCL >> inExpCLn2 >> inExpCLn1 >> inExpCL
	   >> inExpCLp1 >> inExpCLp2;
  }
  textCL.close();
  
  // save CL and CLs for later access:
  calculatedValues[getKey("CL",0,-2)] = inExpCLn2;
  calculatedValues[getKey("CL",0,-1)] = inExpCLn1;
  calculatedValues[getKey("CL",0,0)] = inExpCL;
  calculatedValues[getKey("CL",0,1)] = inExpCLp1;
  calculatedValues[getKey("CL",0,2)] = inExpCLp2;
  calculatedValues[getKey("CL",1,0)] = inObsCL;
  
  calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(inExpCLn2);
  calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(inExpCLn1);
  calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(inExpCL);
  calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(inExpCLp1);
  calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(inExpCLp2);
  calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(inObsCL);
}

/**
   -----------------------------------------------------------------------------
    Check whether the specified category has been defined.
    @param cateScheme - the name of the categorization.
    @returns - true iff the categorization has been defined. 
*/
bool DMTestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool nonExistent = (calculatedValues.find(mapKey) == calculatedValues.end());
  if (nonExistent) {
    std::cout << "DMTestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return !nonExistent;
}
