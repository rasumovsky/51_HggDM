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
    Should probably just hand this thing a workspace. The file location etc.
    should probably go in another wrapper file.
*/
DMTestStat::DMTestStat(TString newJobName, TString newDMSignal,
		       TString newOptions) {
  
  jobName = newJobName;
  DMSignal = newDMSignal;
  options = newOptions;
  allGoodFits = true;
  
  // Map storing all calculations:
  calculatedValues.clear();
  
  // Define the input file, then make a local copy (for remote jobs):
  TString originFile = Form("%s/%s/workspaces/rootfiles/workspaceDM_%s.root",
			    masterOutput.Data(), jobName.Data(), 
			    DMSignal.Data());
  TString copiedFile = Form("workspaceDM_%s.root",DMSignal.Data());
  system(Form("cp %s %s",originFile.Data(),copiedFile.Data()));
  
  // Define and create the output directory:
  outputDir = Form("%s/%s/p0_values/single_files", masterOutput.Data(), 
		   jobName.Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Load the RooWorkspace and ModelConfig:
  TFile input_file(copiedFile,"read");
  workspace = (RooWorkspace*)input_file.Get("combined");
  system(Form("rm %s",input_filename.Data()));
}

double DMTestStat::accessValue(TString testStat, bool observed, int N) {
  // Construct the map key:
  TString currMapKey = testStat;
  if (observed) currMapKey += "_obs";
  else currMapKey += "_exp";
  if (N < 0) currMapKey += Form("_n%d",N);
  else if (N > 0) currMapKey += Form("_p%d",N);
  
  // Check that corresponding entry exists:
  if (mapValueExists(currMapKey)) {
    return calculatedValue[currMapKey];
  }
  else {
    return 0;
  }
}

void DMTestStat::calculateNewCL() {
  // Calculate observed qmu: 
  double muHatObs = 0.0;
  double nllMu1Obs = getFitNLL("obsData", 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL("obsData", 1.0, false, muHatObs);
  double qMuObs = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL("asimovDataMu0", 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL("asimovDataMu0", 0.0, false, muHatExp);
  double qMuExp = getQMuFromNLL(nll_mu1_exp, nll_muhat_exp, muHatExp, 1);
  
  // Calculate CL:
  double expCLn2 = getCLFromQMu(qMuExp, "exp_n2");
  double expCLn1 = getCLFromQMu(qMuExp, "exp_n1");
  double expCLp1 = getCLFromQMu(qMuExp, "exp_p1");
  double expCLp2 = getCLFromQMu(qMuExp, "exp_p2");
  double expCL = getCLFromQMu(qMuExp, "exp");
  double obsCL = getCLFromQMu(qMuObs, "obs");
  
  // Write CL values to file:
  ofstream textCL;
  textCL.open(Form("%s/CL_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  textCL << lambda << " " << lifetime << " " << obsCL << " " << expCLn2
	 << " " << expCLn1 << " " << expCL << " " << expCLp1 << " "
	 << expCLp2 << std::endl;
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
  if (qmu_obs < 0) std::cout << "WARNING! qMuObs < 0 : " << qMuObs << std::endl;
  if (qmu_exp < 0) std::cout << "WARNING! qMuExp < 0 : " << qMuExp << std::endl;
  
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
  textP0.open(Form("%s/p0_values_%s.txt",outputDir.Data(),DMSignal.Data()));
  textP0 << DMSignal << " " << expP0 << " " << obsP0 << std::endl;
  textP0.close();
  
  // Print summary:
  std::cout << "\n  Expected p0 = " << p0_exp << std::endl;
  std::cout << "  Observed p0 = " << p0_obs << std::endl;
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

bool DMTestStat::fitsAllConverged() { 
  return allGoodFits;
}

double getCLsFromCL(double CL) {
  double CLs = 1 - CL;
}

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

double DMTestStat::getCLFromQMu(double qMu, TString type) {
  double CLs = getCLsFromQMu(qMu, type);
  double CL = 1 - CLs;
  return CL;
}

// Formerly calculate_Q0:
double DMTestStat::getQ0FromNLL(double nllMu0, double nllMuhat, double muhat) {
  double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
  return q0;
}

// Formerly calculate_Qmu:
double DMTestStat::getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
				 double muTest) {
  double qmu = 0;
  if (muhat < mutest) {
    qmu = 2 * (nll_mu - nll_muhat);
  }
  else {
    qmu = 0.0;
  }
  return qmu;
}

// Formerly calculate_P0:
double DMTestStat::getP0FromQ0(double q0) {
  double p0 = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0)));
  return p0;
}

// Formerly calculate_PMu:
double DMTestStat::getPMuFromQMu(double qMu) {
  double pMu = 1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu)));
  return pMu;
}

// Formerly calculate_Pb:
double DMTestStat::getPbFromQMu(double qMu, double sigma, double mu) {
  double pB = 1 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMu));
  return pB;
}

// Formerly Pb_fromN:
double DMTestStat::getPbfromN(double N) {
  double pB = 1 - ROOT::Math::gaussian_cdf(N);
  return pB;
}

// Formerly get_fit_nll:
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
    Check whether the specified category has been defined.
    @param cateScheme - the name of the categorization.
    @returns - true iff the categorization has been defined. 
*/
bool DMTestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool exists = (calculatedValue.find(Form("%s_0",mapKey.Data()))
		 == calculatedValue.end());
  if (!exists) {
    std::cout << "DMTestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return exists;
}
