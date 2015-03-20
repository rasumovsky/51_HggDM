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
  
  // Either load the signal parameterization from file or create new ones:
  if (options.Contains("FromFile")) loadSigParamFromFile();
  else createNewSigParam();
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
   Get the name of the output textfile for the given category index. fileType
   can either be "fit" or "yield".
*/
TString DMSigParam::getSigParamFileName(TString production, TString fileType) {
  TString name = Form("%s/%s/%s_%s_%d.txt",outputDir.Data(),production.Data(),
		      fileType.Data(),cateScheme.Data());
  return name;
}

/**
   Create new masspoints by looping over the TTree.
*/
void DMSigParam::createNewSigParam(TString production) {
  std::cout << "DMSigParam: creating new signal fit from tree." << std::endl;
  
  // Create output file for the fit parameters of this process:
  ofstream outputFile(getSigParamFileName(production,"fit"));
  
  // Vectors to store fitted PDFs
  std::vector<RooCBShape*> vectorCB; vectorCB.clear();
  std::vector<RooGaussian*> vectorGA; vectorGA.clear();
  std::vector<RooAddPdf*> vectorSignal; vectorSignal.clear();
  
  // Load the RooDataSet corresponding to the sample
  TString sampleName = prodToSample[production];
  DMMassPoints *dmmp = new DMMassPoints(jobName,sampleName,cateScheme,"New");
  
  // Loop over categories and production modes:
  for (int i_c = 0; i_c < ncategories; i_c++) {
    
    // Use DMMassPoints class to construct the RooDataSet:
    RooDataSet *currData = dmmp->getCateDataSet(i_c);
    
    // Get the observable from the dataset:
    RooRealVar *m_yy;
    TIterator *iterArgs = ((RooArgSet*)currData->get())->createIterator();
    RooRealVar* currIter = NULL;
    while ((currIter = (RooRealVar*)iterArgs->Next())) {
      if (((TString)currIter->GetName()).EqualTo("m_yy")) {
	m_yy = currIter;
	break;
      }
    }
    
    // Define the fit variables:
    RooRealVar *mu = new RooRealVar("mu","mu",);
    RooRealVar *sigmaCB = new RooRealVar("sigmaCB","sigmaCB");
    RooRealVar *sigmaGA = new RooRealVar("sigmaGA","sigmaGA");
    RooRealVar *alpha = new RooRealVar("alpha","alpha");
    RooRealVar *nCB = new RooRealVar("nCB","nCB");
    RooRealVar *frac = new RooRealVar("frac","frac");
    
    // Define the PDFs:
    RooCBShape *currCB = new RooCBShape(Form("CB_%s_%d",cateScheme.Data(),i_c),
					Form("CB_%s_%d",cateScheme.Data(),i_c),
					m_yy, mu, sigmaCB, alpha, nCB);
    
    RooGaussian *currGA = new RooGaussian(Form("GA_%s_%d",cateScheme.Data(),
					       i_c),
					  Form("GA_%s_%d",cateScheme.Data(),
					       i_c),
					  m_yy, mu, sigmaGA);
    
    RooAddPdf *currSignal = new RooAddPdf(Form("Signal_%s_%d",cateScheme.Data(),
					       i_c),
					  Form("Signal_%s_%d",cateScheme.Data(),
					       i_c),
					  currCB, currGA, frac);
    
    // Perform the fits.
    statistics::setDefaultPrintLevel(0);
    RooNLLVar *nLL = (RooNLLVar*)currSignal.createNLL(currData);
    statistics::minimize(nLL);
    
    // Save the parameters to file:
    outputFile << i_c << " " << mu->getVal() << " " << sigmaCB->getVal() << " " 
	       << alpha->getVal() << " " << nCB->getVal() << " "
	       << sigmaGA->getVal() << " " << frac->getVal() << endl;
    
    // Add to vector of shapes:
    vectorCB.push_back(currCB);
    vectorGA.push_back(currGA);
    vectorSignal.push_back(currSignal);
  }
  
  outputFile.close();
  
  // Add signal shapes in all categories to the map.
  sigCB[production] = vectorCB;
  sigGA[production] = vectorGA;
  sigPDF[production] = vectorSignal;
}
