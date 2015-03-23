////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMSigParam.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/03/2015                                                          //
//                                                                            //
//  This class parameterizes the resonance shape of the SM Higgs including    //
//  the SM and DM production modes. For now, the program uses a single mass   //
//  point (125 GeV), and only has the SM production modes.                    //
//                                                                            //
//  NOTE: functionality for a RooCategory has been provided, but has not yet  //
//  been utilized. Consider implementing for easy access to combined PDF.     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMSigParam.h"

/**
   Initialize the DMSigParam class with new observable RooRealVar and
   RooCategory.
   classes (instead of importing them).
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @returns void.
*/
DMSigParam::DMSigParam(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions) {
  RooRealVar *newObservable = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,
					     DMMyyRangeHi);
  DMSigParam(newJobName, newSampleName, newCateScheme, newOptions,
	     newObservable);
}

/**
   Initialize the DMSigParam class with a new RooCategory.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @returns void.
*/
DMSigParam::DMSigParam(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions,
		       RooRealVar *newObservable) {
  
  // Load the selector to get category information.
  selector = new DMEvtSelect();
  nCategories = selector->getNCategories(newCateScheme);
  
  // Define a new RooCategory for the dataset, since none was provided:
  RooCategory *newCategories = new RooCategory(Form("categories_%s",
						    newCateScheme.Data()),
					       Form("categories_%s",
						    newCateScheme.Data()));
  // Loop over categories to define categories:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    newCategories->defineType(Form("%s_%d",newCateScheme.Data(),i_c));
    //newCategories->setRange(Form("rangeName_",i_b,i_r),Form("%s_%d",cateScheme.Data(),i_c));
  }
  
  // Then call the full initializer:
  DMSigParam(newJobName, newSampleName, newCateScheme, newOptions,
	     newObservable, newCategories);
}

/**
   Initialize the DMSigParam class using previously defined observable
   RooRealVar and RooCategory classes.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @param newCategories = The RooCategory to be used in the combined PDF.
   @returns void.
*/
DMSigParam::DMSigParam(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions,
		       RooRealVar *newObservable, RooCategory *newCategories) {
  std::cout << std::endl << "DMSigParam::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign the observable and categorization based on inputs:
  setMassObservable(newObservable);
  setRooCategory(newCategories);
  
  // Get the number of analysis categories if not already done:
  if (!selector) {
    selector = new DMEvtSelect();
    nCategories = selector->getNCategories(newCateScheme);
  }
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/SigParam",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/Plots",outputDir.Data())); 
  system(Form("mkdir -vp %s/all",outputDir.Data()));
  for (int i_p = 0; i_p < nProdModes; i_p++) {
    system(Form("mkdir -vp %s/%s",outputDir.Data(),(sigProdModes[i_p]).Data()));
  }
  
  // Load the signal parameterization from file or start from scratch:
  for (int i_p = 0; i_p < nProdModes; i_p++) {
    createSigParam(sigProdModes[i_p], (!options.Contains("FromFile")));
  }
  return;
}

/**
   Get a pointer to the fitted Crystal Ball component for a particular category
   and production process.
   @param cateIndex - The index of the category for which we want the PDF.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h.
   @returns A pointer to the crystal-ball component of the fit.
*/
RooCBShape* DMSigParam::getCateCrystalBall(int cateIndex, TString process) {
  return (sigCB[process])[cateIndex];
}

/**
   Get a pointer to the fitted Gaussian component for a particular category
   and production process.
   @param cateIndex - The index of the category for which we want the PDF.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h.
   @returns A pointer to the gaussian component of the fit.
 */
RooGaussian* DMSigParam::getCateGaussian(int cateIndex, TString process) {
  return (sigGA[process])[cateIndex];
}

/**
   Get the combined resonance shape for a category and production process.
   @param cateIndex - The index of the category for which we want the PDF.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h.
   @returns A pointer to the total signal shape.
*/
RooAddPdf* DMSigParam::getCateSigPDF(int cateIndex, TString process) {
  return (sigPDF[process])[cateIndex];
}

/**
   Get the signal yield for a particular process in a particular category.
   @param cateIndex - The index of the category for which we want the PDF.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h.
   @returns The signal yield for the specified process in the given category.
*/
double DMSigParam::getCateSigYield(int cateIndex, TString process) {
  return (sigYield[process])[cateIndex];
}

/**
   Get the signal yield for a particular process in all categories.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h.
   @returns The signal yield in all of the categories for the specified process.
*/
double DMSigParam::getCombSigYield(TString process) {
  double sum = 0;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    sum += getCateSigYield(i_c, process);
  }
  return sum;
}

/**
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RooRealVar* DMSigParam::getMassObservable() {
  return m_yy;
}

/**
   Returns a pointer to the RooCategory used in the combined dataset.
   @returns pointer to the RooCategory object.
*/
RooCategory* DMSigParam::getRooCategory() {
  return categories;
}

/**
   Get the value of a particular parameter of the signal PDF. 
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h
   @param param - The fit parameter. Options are: "mu", "sigmaCB", "sigmaGA", 
   "alpha", "nCB", "frac"
   @param cateIndex - The index of the category for which we want the PDF.
   @returns The value of the specified signal parameter. 
*/
double DMSigParam::getSigParam(TString process, TString param, int cateIndex) {
  RooArgSet *currArgs = ((sigPDF[process])[cateIndex])->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    if (((TString)currIter->GetName()).Contains(param)) {
      return currIter->getVal();
      break;
    }
  }
  std::cout << "DMSigParam: requested signal parameter not found." << std::endl;
  return 0.0;
}

/**
   Get the name of the output textfile for the given category index. fileType
   can either be "fit" or "yield".
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h
   @param fileType - The file type (either "yield" or "fit").
   @returns The filename for loading/saving signal parameters (yield and fit).
*/
TString DMSigParam::getSigParamFileName(TString process, TString fileType) {
  TString name = Form("%s/%s/%s_%s.txt",outputDir.Data(),process.Data(),
		      fileType.Data(),cateScheme.Data());
  return name;
}

/**
   Set the pointer to the observable. 
   @param newObservable - The new RooRealVar observable to use for datasets. 
   @returns void.
 */
void DMSigParam::setMassObservable(RooRealVar *newObservable) {
  m_yy = newObservable;
}

/**
   Set the pointer to the RooCategory object. 
   @param newCategories - The new RooCategory to use for the combined dataset. 
   @returns void.
 */
void DMSigParam::setRooCategory(RooCategory *newCategories) {
  categories = newCategories;
}

/**
   Create new signal parameterization using a call to DMMassPoints for input MC.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h
   @param makeNew - Set true if make parameterization from scratch. Else false.
   @returns void.
*/
void DMSigParam::createSigParam(TString process, bool makeNew) {
  std::cout << "DMSigParam: creating new signal fit from tree." << std::endl;
  
  // Create output file or load input file.
  ofstream outputFitFile;
  ofstream outputYieldFile;
  ifstream inputFitFile;
  ifstream inputYieldFile;
  if (makeNew) {
    outputFitFile.open(getSigParamFileName(process,"fit"));
    outputYieldFile.open(getSigParamFileName(process,"yield"));
  }
  else {
    inputFitFile.open(getSigParamFileName(process,"fit"));
    inputYieldFile.open(getSigParamFileName(process,"yield"));
  }
  
  // Vectors to store fitted PDFs
  std::vector<RooCBShape*> vectorCB; vectorCB.clear();
  std::vector<RooGaussian*> vectorGA; vectorGA.clear();
  std::vector<RooAddPdf*> vectorSignal; vectorSignal.clear();
  std::vector<double> vectorYield; vectorYield.clear();
  
  // Load the RooDataSet corresponding to the sample
  TString sampleName = nameToSample[process];
  DMMassPoints *dmmp;
  // Important to provide pointer to m_yy and categories!
  if (makeNew) dmmp = new DMMassPoints(jobName,sampleName,cateScheme,"New",
				       m_yy,categories);
  
  // Loop over categories and process modes:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    
    // WARNING: ALL THE PARAMETER RANGES MUST BE SET:
    
    // Define the fit variables (Can't avoid using >80 char per line...):
    RooRealVar *mu = new RooRealVar(Form("mu_%s_%d",process.Data(),i_c),
				    Form("mu_%s_%d",process.Data(),i_c),
				    0,0,0);
    RooRealVar *sigmaCB = new RooRealVar(Form("sigmaCB_%s_%d",process.Data(),
					      i_c),
					 Form("sigmaCB_%s_%d",process.Data(),
					      i_c),
					 0,0,0);
    RooRealVar *sigmaGA = new RooRealVar(Form("sigmaGA_%s_%d",process.Data(),
					      i_c),
					 Form("sigmaGA_%s_%d",process.Data(),
					      i_c),
					 0,0,0);
    RooRealVar *alpha = new RooRealVar(Form("alpha_%s_%d",process.Data(),i_c),
				       Form("alpha_%s_%d",process.Data(),i_c),
				       0,0,0);
    RooRealVar *nCB = new RooRealVar(Form("nCB_%s_%d",process.Data(),i_c),
				     Form("nCB_%s_%d",process.Data(),i_c),
				     0,0,0);
    RooRealVar *frac = new RooRealVar(Form("frac_%s_%d",process.Data(),i_c),
				      Form("frac_%s_%d",process.Data(),i_c),
				      0,0,0);
    
    // Define the PDFs:
    RooCBShape *currCB = new RooCBShape(Form("CB_%s_%d",process.Data(),i_c),
					Form("CB_%s_%d",process.Data(),i_c),
					*m_yy, mu, sigmaCB, alpha,
					nCB);
    
    RooGaussian *currGA = new RooGaussian(Form("GA_%s_%d",process.Data(),i_c),
					  Form("GA_%s_%d",process.Data(),i_c),
					  *m_yy, mu, sigmaGA);
    
    RooAddPdf *currSignal = new RooAddPdf(Form("sig_%s_%d",process.Data(),i_c),
					  Form("sig_%s_%d",process.Data(),i_c),
					  currCB, currGA, frac);
    
    // If making from scratch, se DMMassPoints to construct the RooDataSet:
    if (makeNew) {
      RooDataSet *currData = dmmp->getCateDataSet(i_c);
      
      // Store the signal yields in memory:
      vectorYield.push_back(currData->sumEntries());
      
      // Perform the fit:
      statistics::setDefaultPrintLevel(0);
      RooNLLVar *nLL = (RooNLLVar*)currSignal->createNLL(currData);
      statistics::minimize(nLL);
      
      // After the fit, set parameters constant:
      mu->setConstant(true);
      sigmaCB->setConstant(true);
      sigmaGA->setConstant(true);
      alpha->setConstant(true);
      nCB->setConstant(true);
      frac->setConstant(true);
      
      // Save the fitted parameters to file:
      outputFitFile << i_c << " " << mu->getVal() << " " << sigmaCB->getVal()
		    << " " << alpha->getVal() << " " << nCB->getVal() << " "
		    << sigmaGA->getVal() << " " << frac->getVal() << std::endl;
      outputYieldFile << i_c << " " << currData->sumEntries() << " " 
		      << currData->numEntries() << std::endl;
    }
    // If using previous parameterization, just load params from .txt file.
    else {
      double rC, rMu, rSigmaCB, rAlpha, rNCB, rSigmaGA, rFrac;
      while (!inputFitFile.eof()) {
	inputFitFile >> rC >> rMu >> rSigmaCB >> rAlpha >> rNCB >> rSigmaGA
		     >> rFrac;
	mu->setVal(rMu);
	sigmaCB->setVal(rSigmaCB);
	alpha->setVal(rAlpha);
	nCB->setVal(rNCB);
	sigmaGA->setVal(rSigmaGA);
	frac->setVal(rFrac);
	break;
      }
      double rSum, rNum;
      while (!inputYieldFile.eof()) {
	inputYieldFile >> rC >> rSum >> rNum;
	vectorYield.push_back(rSum);
	break;
      }
    }
    
    // Save the fitted PDFs:
    vectorCB.push_back(currCB);
    vectorGA.push_back(currGA);
    vectorSignal.push_back(currSignal);
  }
  
  if (makeNew) {
    outputFitFile.close();
    outputYieldFile.close();
  }
  else {
    inputFitFile.close();
    inputYieldFile.close();
  }
  
  // Add signal shapes in all categories to the map.
  sigCB[process] = vectorCB;
  sigGA[process] = vectorGA;
  sigPDF[process] = vectorSignal;
  sigYield[process] = vectorYield;
}

