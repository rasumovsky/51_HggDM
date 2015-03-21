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
////////////////////////////////////////////////////////////////////////////////

#include "DMSigParam.h"

/**
   Initialize the SigParam class.
*/
DMSigParam::DMSigParam(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions) {
  std::cout << std::endl << "DMSigParam::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobname;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
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
    createSigParam(sigProdModes[i_p], (!option.Contains("FromFile")));
  }
  return;
}

/**
   Get a pointer to the fitted Crystal Ball component for a particular category
   and production process.
*/
RooCBShape* DMSigParam::getCateCrystalBall(int cateIndex, TString process) {
  return (sigCB[process])[cateIndex];
}

/** Get a pointer to the fitted Gaussian component for a particular category
    and production process.
 */
RooGaussian* DMSigParam::getCateGaussian(int cateIndex, TString process) {
  return (sigGA[process])[cateIndex];
}

/**
   Get the combined resonance shape for a category and production process.
*/
RooAddPdf* DMSigParam::getCateSigPDF(int cateIndex, TString process) {
  return (sigPDF[process])[cateIndex];
}

/**
   Get the signal yield for a particular process in a particular category.
*/
double DMSigParam::getCateSigYield(int cateIndex, TString process) {
  return (sigYield[process])[cateIndex];
}

/**
   Get the signal yield for a particular process in all categories.
*/
double DMSigParam::getCombSigYield(TString process) {
  double sum = 0;
  for (int i_c = 0; i_c < ncategories; i_c++) {
    sum += getCateSigYield(i_c, process);
  }
  return sum;
}

/**
   Get the value of a particular parameter of the signal PDF. Options for the param argument are: "mu", "sigmaCB", "sigmaGA", "alpha", "nCB", "frac"

*/
double DMSigParam::getSignalParameter(TString process, TString param,
				      int cateIndex) {
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
*/
TString DMSigParam::getSigParamFileName(TString process, TString fileType) {
  TString name = Form("%s/%s/%s_%s_%d.txt",outputDir.Data(),process.Data(),
		      fileType.Data(),cateScheme.Data());
  return name;
}

/**
   Create new masspoints by looping over the TTree.
*/
void DMSigParam::createNewSigParam(TString process, bool makeNew) {
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
  std::vector<double> vectorYield; vectorYields.clear();
  
  // Load the RooDataSet corresponding to the sample
  TString sampleName = prodToSample[process];
  DMMassPoints *dmmp;
  if (makeNew) dmmp = new DMMassPoints(jobName,sampleName,cateScheme,"New");
  
  // Loop over categories and process modes:
  for (int i_c = 0; i_c < ncategories; i_c++) {
    
    // Use DMMassPoints class to construct the RooDataSet:
    RooRealVar *m_yy;
    RooDataSet *currData;
    if (makeNew) {
      currData = dmmp->getCateDataSet(i_c);
      
      // Save the signal yields:
      vectorYield.push_back(currData->sumEntries());
    
      // Get the observable from the dataset:
      /*
	THIS IS VERY WRONG. HONGTAO SAID THERE IS A SIMPLE WAY TO GET OBS.

	TIterator *iterArgs = ((RooArgSet*)currData->get())->createIterator();
	RooRealVar* currIter = NULL;
	while ((currIter = (RooRealVar*)iterArgs->Next())) {
	if (((TString)currIter->GetName()).EqualTo("m_yy")) {
	m_yy = currIter;
	break;
	}
	}
      */
    }
    else {
      m_yy = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,DMMyyRangeHi);
    }

    // Define the fit variables (Can't avoid using >80 char per line...):
    RooRealVar *currMu = new RooRealVar(Form("mu_%s_%d",process.Data(),i_c),Form("mu_%s_%d",process.Data(),i_c));
    RooRealVar *currSigmaCB = new RooRealVar(Form("sigmaCB_%s_%d",process.Data(),i_c),Form("sigmaCB_%s_%d",process.Data(),i_c));
    RooRealVar *currSigmaGA = new RooRealVar(Form("sigmaGA_%s_%d",process.Data(),i_c),Form("sigmaGA_%s_%d",process.Data(),i_c));
    RooRealVar *currAlpha = new RooRealVar(Form("alpha_%s_%d",process.Data(),i_c),Form("alpha_%s_%d",process.Data(),i_c));
    RooRealVar *currNCB = new RooRealVar(Form("nCB_%s_%d",process.Data(),i_c),Form("nCB_%s_%d",process.Data(),i_c));
    RooRealVar *currFrac = new RooRealVar(Form("frac_%s_%d",process.Data(),i_c),Form("frac_%s_%d",process.Data(),i_c));
    
    // Define the PDFs:
    RooCBShape *currCB = new RooCBShape(Form("CB_%s_%d",process.Data(),i_c),
					Form("CB_%s_%d",process.Data(),i_c),
					m_yy, mu, sigmaCB, alpha, nCB);
    
    RooGaussian *currGA = new RooGaussian(Form("GA_%s_%d",process.Data(),i_c),
					  Form("GA_%s_%d",process.Data(),i_c),
					  m_yy, mu, sigmaGA);
    
    RooAddPdf *currSignal = new RooAddPdf(Form("Sig_%s_%d",process.Data(),i_c),
					  Form("Sig_%s_%d",process.Data(),i_c),
					  currCB, currGA, frac);
    
    if (makeNew) {
      // Perform the fits:
      statistics::setDefaultPrintLevel(0);
      RooNLLVar *nLL = (RooNLLVar*)currSignal.createNLL(currData);
      statistics::minimize(nLL);
      
      // Then save the fitted parameters to file:
      outputFitFile << i_c << " " << mu->getVal() << " " << sigmaCB->getVal()
		    << " " << alpha->getVal() << " " << nCB->getVal() << " "
		    << sigmaGA->getVal() << " " << frac->getVal() << std::endl;
      outputYieldFile << i_c << " " << currData->sumEntries() << " " 
		      << currData->numEntries() << std::endl;
    }
    else {
      // THIS MUST BE FIXED ASAP!
      while (!inputFitFile.eof()) {
	inputFitFile >> rC >> mu->getVal() >> sigmaCB->getVal() 
		     >> alpha->getVal() >> nCB->getVal() >> sigmaGA->getVal() 
		     >> frac->getVal();
      
	mu->setVal(rMu);
	sigmaCB->setVal(rSigmaCB);
	alpha->setVal(rAlpha);
	nCB->setVal(rNCB);
	sigmaGA->setVal(rSigmaGA);
	frac->setVal(rFrac);
	break;
      }

      while (!inputYieldFile.eof()) {
	inputYieldFile >> rC >> currData->sumEntries() >> currData->numEntries() >> std::endl;
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
  
  // Add signal shapes in all categories to the map.
  sigCB[process] = vectorCB;
  sigGA[process] = vectorGA;
  sigPDF[process] = vectorSignal;
  sigYield[process] = vectorYield;
}

