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
   Calculate the test statistic, using the nominal (expected) location of the
   workspace.
   @param newJobName - The name of the job
   @param newDMSignal - The Dark Matter signal to incorporate in the model.
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile"), etc.
*/
DMTestStat::DMTestStat(TString newJobName, TString newDMSignal,
		       TString newCateScheme, TString newOptions) {
  
  // Load the workspace from the nominal location.
  DMWorkspace *dmw = new DMWorkspace(newJobName, newDMSignal, newCateScheme,
				     "FromFile");
  RooWorkspace *newWorkspace = dmw->getCombinedWorkspace();
  DMTestStat(newJobName, newDMSignal, newOptions, newWorkspace);
  return;
}

/**
   Constructor for the DMTestStat class. 
   @param newJobName - The name of the job
   @param newDMSignal - The Dark Matter signal to incorporate in the model.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
DMTestStat::DMTestStat(TString newJobName, TString newDMSignal,
		       TString newOptions, RooWorkspace *newWorkspace) {

  // Assign input variables: 
  jobName = newJobName;
  DMSignal = newDMSignal;
  options = newOptions;
  allGoodFits = true;
  workspace = newWorkspace;
  
  // Map storing all calculations:
  calculatedValues.clear();
  
  // Create output directories:
  outputDir = Form("%s/TestStatistics/",masterOutput.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/CL/",outputDir.Data()));
  system(Form("mkdir -vp %s/p0/",outputDir.Data()));

  // Make new or load old values:
  if (options.Contains("FromFile")) {
    loadStatsFromFile();
  }
  else {
    calculateNewCL();
    calculateNewP0();
  }
  return;
}

/**
   Get the value of one of the test statistics.
   @param testStat - the test stat. name (p0, CL, CLs).
   @param observed - true iff observed, false if expected. 
   @param N - the standard deviation (-2, -1, 0, +1, +2). 
   @returns - the test statistic value.
*/
double DMTestStat::accessValue(TString testStat, bool observed, int N) {
  // Construct the map key:
  TString currMapKey = testStat;
  if (observed) currMapKey += "_obs";
  else currMapKey += "_exp";
  if (N < 0) currMapKey += Form("_n%d",N);
  else if (N > 0) currMapKey += Form("_p%d",N);
  
  // Check that corresponding entry exists:
  if (mapValueExists(currMapKey)) {
    return calculatedValues[currMapKey];
  }
  else {
    return 0;
  }
}

/**
   Calculate the CL and CLs values using model fits.
*/
void DMTestStat::calculateNewCL() {
  // Calculate observed qmu: 
  double muHatObs = 0.0;
  double nllMu1Obs = getFitNLL("obsData", 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL("obsData", 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL("asimovDataMu0", 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL("asimovDataMu0", 0.0, false, muHatExp);
  double expQMu = getQMuFromNLL(nllMu1Exp, nllMuHatExp, muHatExp, 1);
  
  // Calculate CL:
  double expCLn2 = getCLFromQMu(expQMu, "exp_n2");
  double expCLn1 = getCLFromQMu(expQMu, "exp_n1");
  double expCLp1 = getCLFromQMu(expQMu, "exp_p1");
  double expCLp2 = getCLFromQMu(expQMu, "exp_p2");
  double expCL = getCLFromQMu(expQMu, "exp");
  double obsCL = getCLFromQMu(obsQMu, "obs");
  
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
  calculatedValues["CL_exp_n2"] = expCLn2;
  calculatedValues["CL_exp_n1"] = expCLn1;
  calculatedValues["CL_exp"] = expCL;
  calculatedValues["CL_exp_p1"] = expCLp1;
  calculatedValues["CL_exp_p2"] = expCLp2;
  calculatedValues["CL_obs"] = obsCL;
  
  calculatedValues["CLs_exp_n2"] = getCLsFromCL(expCLn2);
  calculatedValues["CLs_exp_n1"] = getCLsFromCL(expCLn1);
  calculatedValues["CLs_exp"] = getCLsFromCL(expCL);
  calculatedValues["CLs_exp_p1"] = getCLsFromCL(expCLp1);
  calculatedValues["CLs_exp_p2"] = getCLsFromCL(expCLp2);
  calculatedValues["CLs_obs"] = getCLsFromCL(obsCL);
}

/**
   Calculate the p0 value using model fits.
*/
void DMTestStat::calculateNewP0() {
  // Calculate observed q0: 
  double muHatObs = 0.0;
  double nllMu0Obs = getFitNLL("obsData", 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL("obsData", 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL("asimovDataMu1", 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL("asimovDataMu1", 0.0, false, muHatExp);
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
  calculatedValues["p0_obs"] = obsP0;
  calculatedValues["p0_exp"] = expP0;
}

/**
   Check if all of the fits done by this class have converged.
   @returns - true iff. all of the fits have been successfully convergent.
*/
bool DMTestStat::fitsAllConverged() { 
  return allGoodFits;
}

/**
   Get the CLs value from CL.
   @param CL - the CL value to convert to CLs.
   @returns - the corresponding CLs value.
*/
double DMTestStat::getCLsFromCL(double CL) {
  return (1.0 - CL);
}

/**
   Get the CL value from CLs.
   @param CLs - the CLs value to convert to CL.
   @returns - the corresponding CL value.
*/
double DMTestStat::getCLFromCLs(double CLs) {
  return (1.0 - CLs);
}

/**
   Get the CLs value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param type - "exp_n2,... exp_p2, "exp", "obs"
   @returns - the CLs value.
*/
double DMTestStat::getCLsFromQMu(double qMu, TString type) {
  double N = 0.0;// default for exp or obs
  if (type.Contains("exp_n2")) N = -2.0;
  else if (type.Contains("exp_n1")) N = -1.0;
  else if (type.Contains("exp_p1")) N = 1.0;
  else if (type.Contains("exp_p2")) N = 2.0;
  double pMu = getPMuFromQMu(qMu);
  double pB = getPbfromN(N);
  double CLs = pMu / (1.0 - pB);
  return CLs;
}

/**
   Get the CL value using qMu and the type.
   @param qMu - the value of the test statistic.
   @param type - "exp_n2,... exp_p2, "exp", "obs"
   @returns - the CLs value.
*/
double DMTestStat::getCLFromQMu(double qMu, TString type) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, type));
  return CL;
}

// Formerly calculate_Q0:
/**
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
  if (muHat < muTest) {
    qMu = 2 * (nllMu - nllMuHat);
  }
  else {
    qMu = 0.0;
  }
  return qMu;
}

/**
   Calculate the value of p0 based on the test statistic q0.
   @param q0 - the test statistic q0.
   @returns - the value of p0.
*/
double DMTestStat::getP0FromQ0(double q0) {
  double p0 = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0)));
  return p0;
}

/**
   Calculate the value of pMu.
   @param qMu - the test statistic qMu.
   @returns - the value of pMu.
*/
double DMTestStat::getPMuFromQMu(double qMu) {
  double pMu = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu)));
  return pMu;
}

/**
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
   Calculate pB based on the standard deviation.
   @param N - the standard deviation.
   @returns - the value of pB.
*/
double DMTestStat::getPbfromN(double N) {
  double pB = 1 - ROOT::Math::gaussian_cdf(N);
  return pB;
}

/**
   Get the negative-log-likelihood for a fit of a specified type to a specified
   dataset.
   @param datasetName - the name of the dataset in the workspace.
   @param muVal - the mu value to fix.
   @param fixMu - true if mu should be fixed to the specified value.
   @param &profiledMu - the profiled value of mu (passed by reference)
   @returns - the nll value.
*/
double DMTestStat::getFitNLL(TString datasetName, double muVal, bool fixMu,
			     double& profiledMu) { 
  std::cout << "getFitNLL( " << datasetName << ", " << muVal << ", " << fixMu
	    << " )" << endl;
  
  // Load the relevant quantities from the ModelConfig:
  ModelConfig* mc = (ModelConfig*)workspace->obj("ModelConfig");
  RooAbsPdf* combPdf = mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  RooArgSet* Observables = (RooArgSet*)mc->GetObservables();
  workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)mc->GetNuisanceParameters()->snapshot();
  RooArgSet* poi = (RooArgSet*)mc->GetParametersOfInterest();
  RooRealVar* firstpoi = (RooRealVar*)poi->first();
  RooAbsData *data = workspace->data(datasetName);
  
  // WARNING! IT IS VERY IMPORTANT THAT NUIS AND GLOBOBS ARE SET PROPERLY HERE!
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet( nuisanceParameters, false, origValNP );
  // the global observables should be fixed to the nominal values...
  statistics::constSet( globalObservables, true );
  
  firstpoi->setVal(muVal);
  firstpoi->setConstant(fixMu);
  
  /*
  // loop over nuisance parameters, set theory to -1 sigma in accordance with ATLAS SUSY WG policy for 2D contours.
  int number_params = nuisanceParameters->getSize();
  TIterator *iter_nuis = nuisanceParameters->createIterator();
  RooRealVar* parg_nuis = NULL;
  while( (parg_nuis = (RooRealVar*)iter_nuis->Next()) )
  {
    TString name = parg_nuis->GetName();
    if( name.Contains("scale_PDF") )
    {
      cout << "  Setting scale_PDF NP to 0" << endl;
      parg_nuis->setVal(0.0);
      parg_nuis->setConstant(true);
    }
  }
  */

  int status = 0; 
  RooNLLVar* varNLL = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize(status, varNLL, "", NULL, false);
  if (status != 0) {
    allGoodFits = false;
  }
  profiledMu = fixMu ? muVal : firstpoi->getVal();
  double nll_value = varNLL->getVal();
  delete varNLL;
  
  // release nuisance parameters after fit and recovery the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  return nll_value;
}

/**
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
  calculatedValues["p0_obs"] = inObsP0;
  calculatedValues["p0_exp"] = inExpP0;
  
  // Read CL values:
  double inObsCL, inExpCLn2, inExpCLn1, inExpCL, inExpCLp1, inExpCLp2;
  while (!textCL.eof()) {
    textCL >> inName >> inObsCL >> inExpCLn2 >> inExpCLn1 >> inExpCL
	   >> inExpCLp1 >> inExpCLp2;
  }
  textCL.close();
  
  // save CL and CLs for later access:
  calculatedValues["CL_exp_n2"] = inExpCLn2;
  calculatedValues["CL_exp_n1"] = inExpCLn1;
  calculatedValues["CL_exp"] = inExpCL;
  calculatedValues["CL_exp_p1"] = inExpCLp1;
  calculatedValues["CL_exp_p2"] = inExpCLp2;
  calculatedValues["CL_obs"] = inObsCL;
  calculatedValues["CLs_exp_n2"] = getCLsFromCL(inExpCLn2);
  calculatedValues["CLs_exp_n1"] = getCLsFromCL(inExpCLn1);
  calculatedValues["CLs_exp"] = getCLsFromCL(inExpCL);
  calculatedValues["CLs_exp_p1"] = getCLsFromCL(inExpCLp1);
  calculatedValues["CLs_exp_p2"] = getCLsFromCL(inExpCLp2);
  calculatedValues["CLs_obs"] = getCLsFromCL(inObsCL);
}

/** 
    Check whether the specified category has been defined.
    @param cateScheme - the name of the categorization.
    @returns - true iff the categorization has been defined. 
*/
bool DMTestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool exists = (calculatedValues.find(Form("%s_0",mapKey.Data()))
		 == calculatedValues.end());
  if (!exists) {
    std::cout << "DMTestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return exists;
}
