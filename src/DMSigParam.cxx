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
//  NOTE: Fit now uses the same mean value for the Crystal Ball and Gaussian  //
//  components. 
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMSigParam.h"

/**
   Initialize the DMSigParam class with a new RooCategory.
   @param newJobName - The name of the job 
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @returns void.
*/
DMSigParam::DMSigParam(TString newJobName, TString newCateScheme,
		       TString newOptions, RooRealVar *newObservable) {
  std::cout << "\nDMSigParam::Initializing..."
	    << "\n\tjobName = " << newJobName
	    << "\n\tcateScheme = " << newCateScheme 
	    << "\n\toptions = " << newOptions << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  CommonFunc::SetAtlasStyle();
  
  if (newObservable == NULL) {
    m_yy = new RooRealVar("m_yy", "m_yy", DMMyyRangeLo, DMMyyRangeHi);
  }
  else {
    // Assign the observable and categorization based on inputs:
    setMassObservable(newObservable);
  }
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMSigParam", masterOutput.Data(), jobName.Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  system(Form("mkdir -vp %s/Plots", outputDir.Data())); 
  
  // Load the SM signal parameterization from file or start from scratch:
  for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
    system(Form("mkdir -vp %s/%s",outputDir.Data(),(sigSMModes[i_SM]).Data()));
    createSigParam(sigSMModes[i_SM], (!options.Contains("FromFile")));
  }
  
  // Load the DM signal parameterization from file or start from scratch:
  for (int i_DM = 0; i_DM < nDMModes; i_DM++) {
    system(Form("mkdir -vp %s/%s",outputDir.Data(),(sigDMModes[i_DM]).Data()));
    createSigParam(sigDMModes[i_DM], (!options.Contains("FromFile")));
  }
  
  // Also create the total SM parameterization:
  system(Form("mkdir -vp %s/SM",outputDir.Data()));
  createSigParam("SM", (!options.Contains("FromFile")));
  
  std::cout << "DMSigParam: Successfully initialized!" << std::endl;
  return;
}

/**
   Add the signal PDF for a particular process and category to the workspace.
   Also include the important energy scale and energy resolution systematic
   uncertainties. 
   @param workspace - The workspace to which we are adding the signal PDF.
   @param namesPES - The names of the photon energy scale systematics.
   @param namesPER - The names of the photon resolution systematics.
   @param process - The signal production process.
   @param cateIndex - The index of the current analysis category in cateScheme. 
   @returns void. 
*/
void DMSigParam::addSigToCateWS(RooWorkspace *&workspace,
				std::vector<TString> namesPES,
				std::vector<TString> namesPER, TString process,
				int cateIndex) {
  std::cout << "DMSigParam: Adding " << process << " to workspace" << std::endl;
  
  // A variable that converts any DM process name into "DM" for simplicity:
  TString processType = (isDMSample(process)) ? "DM" : process;
  
  // Create list of ess to multiply:
  TString listPES = ",";//added comma
  for (int i_e = 0; i_e < (int)namesPES.size(); i_e++) {
    TString atlasExpNamePES = Form("atlas_expected_%s", namesPES[i_e].Data());
    if (!(bool)workspace->obj(atlasExpNamePES)) {
      workspace->factory(Form("%s[1]",atlasExpNamePES.Data()));
    }
    if (i_e < ((int)namesPES.size()-1)) {
      listPES.Append(Form("%s,",atlasExpNamePES.Data()));//added comma
    }
    else {
      listPES.Append(Form("%s",atlasExpNamePES.Data()));//no comma
    }
  }
  if (listPES.EqualTo(",")) listPES = "";
  
  // Create list of res to multiply (process dependent implementation):
  TString listPER = ",";//added comma
  for (int i_r = 0; i_r < (int)namesPER.size(); i_r++) {
    TString atlasExpNamePER = Form("atlas_expected_%s",namesPER[i_r].Data());
    if (!(bool)workspace->obj(atlasExpNamePER)) {
      workspace->factory(Form("%s%s[1]",atlasExpNamePER.Data(),
			      processType.Data()));
    }
    if (i_r < ((int)namesPER.size()-1)) {
      listPER.Append(Form("%s%s,",atlasExpNamePER.Data(),processType.Data()));
    }
    else {
      listPER.Append(Form("%s%s",atlasExpNamePER.Data(),processType.Data()));
    }
  }
  if (listPER.EqualTo(",")) listPER = "";
  
  TString meanCB = Form("%f", getSigParam(process, "meanCB", cateIndex));
  TString massResCB = Form("%f", getSigParam(process, "sigmaCB", cateIndex));
  TString alphaCB = Form("%f", getSigParam(process, "alphaCB", cateIndex));
  TString nCB = Form("%f", getSigParam(process, "nCB", cateIndex));
  //TString meanGA = Form("%f", getSigParam(process, "meanGA", cateIndex));
  TString meanGA = Form("%f", getSigParam(process, "meanCB", cateIndex));
  TString massResGA = Form("%f", getSigParam(process, "sigmaGA", cateIndex));
  TString frac = Form("%f", getSigParam(process, "frac", cateIndex));
  
  workspace->factory(Form("RooCBShape::pdfCB%s(m_yy, prod::meanCB%s(meanCBNom%s[%s]%s), prod::massResCB%s(massResNomCB%s[%s]%s), alphaCB%s[%s], nCB%s[%s])", processType.Data(), processType.Data(), processType.Data(), meanCB.Data(), listPES.Data(), processType.Data(), processType.Data(), massResCB.Data(), listPER.Data(), processType.Data(), alphaCB.Data(), processType.Data(), nCB.Data()));
  
  //workspace->factory(Form("RooGaussian::pdfGA%s(m_yy, prod::meanGA%s(meanGANom%s[%s]%s), prod::massResGA%s(massResNomGA%s[%s]%s))", processType.Data(), processType.Data(), processType.Data(), meanGA.Data(), listPES.Data(), processType.Data(), processType.Data(), massResGA.Data(), listPER.Data()));
  workspace->factory(Form("RooGaussian::pdfGA%s(m_yy, prod::meanGA%s(meanCBNom%s[%s]%s), prod::massResGA%s(massResNomGA%s[%s]%s))", processType.Data(), processType.Data(), processType.Data(), meanCB.Data(), listPES.Data(), processType.Data(), processType.Data(), massResGA.Data(), listPER.Data()));
  
  workspace->factory(Form("SUM::sigPdf%s(frac%s[%s]*pdfCB%s,pdfGA%s)", processType.Data(), processType.Data(), frac.Data(), processType.Data(), processType.Data())); 
  
  std::cout << "DMSigParam: Finished adding " << process << " to workspace."
	    << std::endl;
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
  for (int i_c = 0; i_c < DMAnalysis::getNumCategories(cateScheme); i_c++) {
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
   Get the value of a particular parameter of the signal PDF. 
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h
   @param param - The fit parameter. Options are: "meanCB", "sigmaCB", "meanGA",
   "sigmaGA", "alphaCB", "nCB", "frac"
   @param cateIndex - The index of the category for which we want the PDF.
   @returns The value of the specified signal parameter. 
*/
double DMSigParam::getSigParam(TString process, TString param, int cateIndex) {
  // Require the Gaussian to have the same mean as the CB.
  if (param.EqualTo("meanGA")) param = "meanCB";
  RooArgSet *currArgs = ((sigPDF[process])[cateIndex])->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    if (((TString)currIter->GetName()).Contains(param)) {
      //std::cout << "DMSigParam::getSigParam(" << process << ", " << param
      //<< ", " << cateIndex << ")=" << currIter->getVal() << std::endl;
      return currIter->getVal();
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
TString DMSigParam::getSigParamFileName(TString process, TString fileType,
					int cateIndex) {
  TString name = Form("%s/%s/%s_%s_%d.txt", outputDir.Data(), process.Data(),
		      fileType.Data(), cateScheme.Data(), cateIndex);
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
   Create new signal parameterization using a call to DMMassPoints for input MC.
   @param process - The signal production process of interest. Possibilities
   are listed in DMHeader.h
   @param makeNew - Set true if make parameterization from scratch. Else false.
   @returns - void.
*/
void DMSigParam::createSigParam(TString process, bool makeNew) {
  if (makeNew) {
    std::cout << "DMSigParam: Make new " << process << " param." << std::endl;
  }
  else {
    std::cout << "DMSigParam: Load " << process << " param." << std::endl;
  }
  
  // Vectors to store fitted PDFs
  std::vector<RooCBShape*> vectorCB; vectorCB.clear();
  std::vector<RooGaussian*> vectorGA; vectorGA.clear();
  std::vector<RooAddPdf*> vectorSignal; vectorSignal.clear();
  std::vector<double> vectorYield; vectorYield.clear();
  
  // Load the RooDataSet corresponding to the sample:
  DMMassPoints *mp;
  DMMassPoints *mps[nSMModes];
  
  if (makeNew) {
    // For total SM, load all SM mass points.
    if (process.EqualTo("SM")) {
      for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
	mps[i_SM] = new DMMassPoints(jobName, sigSMModes[i_SM], cateScheme,
				     "FromFile", m_yy);
      }
    }
    // Otherwise, just load a particular mode.
    else {
      mp = new DMMassPoints(jobName, process, cateScheme, "FromFile", m_yy);
    }
  }
  
  // Loop over categories and process modes:
  for (int i_c = 0; i_c < DMAnalysis::getNumCategories(cateScheme); i_c++) {
        
    // Create output file or load input file.
    ofstream outputFitFile;
    ofstream outputYieldFile;
    ifstream inputFitFile;
    ifstream inputYieldFile;
    if (makeNew) {
      outputFitFile.open(getSigParamFileName(process,"fit",i_c));
      outputYieldFile.open(getSigParamFileName(process,"yield",i_c));
    }
    else {
      inputFitFile.open(getSigParamFileName(process,"fit",i_c));
      inputYieldFile.open(getSigParamFileName(process,"yield",i_c));
    
      // Check that the input files exist, and if they don't, make new param. 
      if (!inputFitFile || !inputYieldFile) {
	createSigParam(process, true);
	return;
      }
    }
    
    // Define the fit variables (Can't avoid using >80 char per line...):
    RooRealVar *meanCB = new RooRealVar(Form("meanCB_%s_%d",process.Data(),i_c),
					Form("meanCB_%s_%d",process.Data(),i_c),
					125.4,123.0,128.0);
    RooRealVar *sigmaCB = new RooRealVar(Form("sigmaCB_%s_%d",process.Data(),
					      i_c),
					 Form("sigmaCB_%s_%d",process.Data(),
					      i_c),
					 1.5,0.1,10.0);
    RooRealVar *alpha = new RooRealVar(Form("alphaCB_%s_%d",process.Data(),i_c),
				       Form("alphaCB_%s_%d",process.Data(),i_c),
				       1.4,0.1,10.0);
    RooRealVar *nCB = new RooRealVar(Form("nCB_%s_%d",process.Data(),i_c),
				     Form("nCB_%s_%d",process.Data(),i_c),
				     10,0.001,20);
    //RooRealVar *meanGA = new RooRealVar(Form("meanGA_%s_%d",process.Data(),i_c),
    //					Form("meanGA_%s_%d",process.Data(),i_c),
    //					125.4,123.0,128.0);
    RooRealVar *sigmaGA = new RooRealVar(Form("sigmaGA_%s_%d",process.Data(),
					      i_c),
					 Form("sigmaGA_%s_%d",process.Data(),
					      i_c),
					 3,0.1,20.0);
    RooRealVar *frac = new RooRealVar(Form("frac_%s_%d",process.Data(),i_c),
				      Form("frac_%s_%d",process.Data(),i_c),
				      0.9,0.001,1.0);
    
    // Define the PDFs:
    RooCBShape *currCB = new RooCBShape(Form("pdfCB_%s_%d",process.Data(),i_c),
					Form("pdfCB_%s_%d",process.Data(),i_c),
					*m_yy, *meanCB, *sigmaCB, *alpha, *nCB);
    RooGaussian *currGA = new RooGaussian(Form("pdfGA_%s_%d",process.Data(),
					       i_c),
					  Form("pdfGA_%s_%d",process.Data(),
					       i_c),
					  *m_yy, *meanCB, *sigmaGA);
					  //*m_yy, *meanGA, *sigmaGA);
    RooAddPdf *currSignal = new RooAddPdf(Form("sigPdf_%s_%d",process.Data(),
					       i_c),
					  Form("sigPdf_%s_%d",process.Data(),
					       i_c),
					  *currCB, *currGA, *frac);
    
    // If making from scratch, use DMMassPoints to construct the RooDataSet:
    if (makeNew) {
      RooDataSet *currData = NULL;
      // For the SM dataset, add up datasets for all production modes:
      if (process.EqualTo("SM")) {
	for (int i_SM = 0; i_SM < nSMModes; i_SM++) {
	  if (i_SM == 0) {
	    currData = mps[i_SM]->getCateDataSet(i_c);
	  }
	  else {
	    RooDataSet *tempData = mps[i_SM]->getCateDataSet(i_c);
	    currData->append(*tempData);
	  }
	}
      }
      // Otherwise just load the process-specific dataset:
      else {
	currData = mp->getCateDataSet(i_c);
      }
      
      std::cout << "DMSigParam: Starting signal fit on dataset: " << std::endl;
      currData->Print("v");
      
      // Store the signal yields in memory:
      vectorYield.push_back(currData->sumEntries());
      
      // Perform the fit:
      statistics::setDefaultPrintLevel(0);
      RooNLLVar *nLL = (RooNLLVar*)currSignal->createNLL(*currData);
      statistics::minimize(nLL);
      
      // After the fit, set parameters constant:
      meanCB->setConstant(true);
      sigmaCB->setConstant(true);
      alpha->setConstant(true);
      nCB->setConstant(true);
      //meanGA->setConstant(true);
      sigmaGA->setConstant(true);
      frac->setConstant(true);
      
      // Save the fitted parameters to file:
      outputFitFile << i_c << " "
		    << meanCB->getVal()  << " "
		    << sigmaCB->getVal() << " "
		    << alpha->getVal()   << " "
		    << nCB->getVal()     << " "
	//<< meanGA->getVal()  << " "
		    << meanCB->getVal()  << " "
		    << sigmaGA->getVal() << " "
		    << frac->getVal() << std::endl;
      outputYieldFile << i_c << " "
		      << currData->sumEntries() << " " 
		      << currData->numEntries() << std::endl;
      
      // Then plot the parameterization and the data used for the fit:
      TCanvas *can = new TCanvas();
      RooPlot* frame = m_yy->frame(Bins(40), Range(higgsMass-10, higgsMass+10));
      currData->plotOn(frame);
      currSignal->plotOn(frame, LineColor(2));
      currSignal->plotOn(frame, Components(*currCB), LineColor(4));
      currSignal->plotOn(frame, Components(*currGA), LineColor(3));
      frame->SetYTitle("Events/0.5 GeV");
      frame->SetXTitle("M_{#gamma#gamma} [GeV]");
      frame->Draw();
      TLatex text; text.SetNDC(); text.SetTextColor(1);
      text.DrawLatex(0.2, 0.76, Form("Category %d", i_c));
      text.DrawLatex(0.2, 0.82, Form("Signal %s", process.Data()));
      TH1F *histGA = new TH1F("histGA", "histGA", 1, 0, 1);
      TH1F *histCB = new TH1F("histCB", "histCB", 1, 0, 1);
      TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
      histGA->SetLineColor(3);
      histCB->SetLineColor(4);
      histSig->SetLineColor(2);
      TLegend leg(0.2, 0.6, 0.4, 0.74);
      leg.SetFillColor(0);
      leg.SetTextSize(0.04);
      leg.SetBorderSize(0);
      leg.AddEntry(histGA, "Gaussian", "l");
      leg.AddEntry(histCB, "Crystal-Ball", "l");
      leg.AddEntry(histSig, "Total Signal", "l");
      leg.Draw("SAME");
      
      can->Print(Form("%s/Plots/%s_%s_%d.eps", outputDir.Data(), process.Data(),
		      cateScheme.Data(), i_c));
    }
    // If using previous parameterization, just load params from .txt file.
    else {
      double rC, rMeanCB, rSigmaCB, rAlpha, rNCB, rMeanGA, rSigmaGA, rFrac;
      while (!inputFitFile.eof()) {
	inputFitFile >> rC >> rMeanCB >> rSigmaCB >> rAlpha >> rNCB
		     >> rMeanGA >> rSigmaGA >> rFrac;
	
	meanCB->setVal(rMeanCB);
	sigmaCB->setVal(rSigmaCB);
	alpha->setVal(rAlpha);
	nCB->setVal(rNCB);
	//meanGA->setVal(rMeanGA);
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
    
    // Save the fitted PDFs for this category:
    vectorCB.push_back(currCB);
    vectorGA.push_back(currGA);
    vectorSignal.push_back(currSignal);

    // Close input/output files for this category.
    if (makeNew) {
      outputFitFile.close();
      outputYieldFile.close();
    }
    else {
      inputFitFile.close();
      inputYieldFile.close();
    }
  }// end of loop over categories.
  
  // Add signal shapes in all categories to the map.
  sigCB[process] = vectorCB;
  sigGA[process] = vectorGA;
  sigPDF[process] = vectorSignal;
  sigYield[process] = vectorYield;
}

