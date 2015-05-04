////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMBkgModel.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/03/2015                                                          //
//                                                                            //
//  This class implements a broad range of fit functions to use for the       //
//  background model in the H->yy analysis.                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMBkgModel.h"

/**
   Initialize the DMBkgModel class and make a new RooCategory.
   @param newJobName - The name of the job 
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @returns void.
*/
DMBkgModel::DMBkgModel(TString newJobName, TString newCateScheme,
		       TString newOptions, RooRealVar *newObservable) {
  std::cout << "\nDMBkgModel::Initializing..." 
	    << "\n\tjobName = " << newJobName
	    << "\n\tcateScheme = " << newCateScheme 
	    << "\n\toptions = " << newOptions << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign the observable based on inputs:
  if (newObservable == NULL) {
    m_yy = new RooRealVar("m_yy", "m_yy", DMMyyRangeLo, DMMyyRangeHi);
  }
  else {
    setMassObservable(newObservable);
  }
  return;
}

/**
   Add the chosen background model to the workspace provided. Also add the 
   associated nuisance parameters to the nuisParams set.
   @param workspace - The workspace to which the PDFs will be added.
   @param nuisParams - The set of nuisance parameters to which the background
                       parameters will be added.
   @param cateIndex - The index of the current analysis category. 
*/
void DMBkgModel::addBkgToCateWS(RooWorkspace *&workspace, 
				RooArgSet *&nuisParams, int cateIndex) {
  std::cout << "DMBkgModel: Adding category " << cateIndex
	    << " background model to the workspace" << std::endl;
  
  // First, call the function to get the background PDF:
  RooAbsPdf* currBkgModel = getCateBkgPDF(cateIndex);
  
  // Then add it to the workspace:
  workspace->import(*currBkgModel);
  
  // Then add the parameters to the workspace:
  RooArgSet *currArgs = currBkgModel->getVariables();
  TIterator *iterArgs = currArgs->createIterator();
  RooRealVar* currIter = NULL;
  while ((currIter = (RooRealVar*)iterArgs->Next())) {
    nuisParams->add(*currIter);
  }
  
  // Finally, include a normalization parameter for the background:
  workspace->factory("nBkg[100,0,1000000]");
  nuisParams->add(*workspace->var("nBkg"));
  
  std::cout << "DMBkgModel: Finished adding category " << cateIndex
	    << " background model to the workspace" << std::endl;
}

/**
   Get the PDF for a given class.
   @param cateIndex - The category of cateScheme used in class initialization.
   @returns The corresponding background PDF for the analysis.
*/
RooAbsPdf* DMBkgModel::getCateBkgPDF(int cateIndex) {
  TString cateName = Form("%s_%d",cateScheme.Data(),cateIndex);
  TString currFunction = DMAnalysis::cateToBkgFunc(cateName);
  TString currFuncName = Form("bkg_%d",cateIndex);
  return getBkgPDFByName(currFuncName, currFunction);
}

/**
   Get a PDF for the background of a specified type and order. 
   @param name - The name of the PDF
   @param fitFunc - The function type.
                    "Bern" = Bernstein polynomial
		    "Exppol" = Exponential polynomial
   @param order - The order of the function to use.
   @returns The background PDF as a pointer to a RooAbsPdf. 
*/
RooAbsPdf* DMBkgModel::getBkgPDFByName(TString fitName, TString fitFunc) {
  std::cout << "DMBkgModel: Constructing " << fitFunc << std::endl;
  int order = getOrderFromFunc(fitFunc);
  std::cout << "DMBkgModel: Function of order " << order << std::endl;
  
  // Set the range of the m_yy variable from DMHeader.h:
  RooConstVar min("min", "min", DMMyyRangeLo);
  RooConstVar max("max", "max", DMMyyRangeHi);
  
  // Background fit variables:
  RooArgList *bkgArgs = new RooArgList();
  RooRealVar *pVar[10];
  RooRealVar *cVar[10];
  TString expFitFormat = "exp(";
  
  // Pointers to the background PDFs:
  RooAbsPdf *background;
  RooBernsteinM *bern;
  RooGenericPdf *exppol;
  
  if (fitFunc.Contains("Exppol")) {
    bkgArgs->add(*m_yy);
  }
  
  // Loop over the order of the function to define parameters:
  for (int i_p = 0; i_p <= order; i_p++) {
    // Parameters for Bernstein polynomial:
    if (fitFunc.Contains("Bern")) {
      if (i_p == 0) {
	pVar[i_p] = new RooRealVar(Form("pVar%d",i_p), 
				   Form("pVar%d",i_p), 1);
	pVar[i_p]->setConstant(true);
      }
      else {
	pVar[i_p] = new RooRealVar(Form("pVar%d",i_p), Form("pVar%d",i_p),
				   0.1, 0.0, 10.0);
	pVar[i_p]->setConstant(false);
      }
      bkgArgs->add(*pVar[i_p]);
    }
    // Parameters for exponential polynomial:
    else if (fitFunc.Contains("Exppol") && i_p < order) {
      cVar[i_p] = new RooRealVar(Form("cVar%d",i_p), Form("cVar%d",i_p),
				 0.0, -1.0, 1.0 );
      cVar[i_p]->setConstant(false);
      bkgArgs->add(*cVar[i_p]);
      expFitFormat += Form("@%d",i_p+1);
      for (int i_t = 0; i_t <= i_p; i_t++) {
	expFitFormat += "*(@0-100)";
      }
    }
  }
  expFitFormat += ")";
  
  std::cout << "Printing bkgArgs:" << std::endl;
  bkgArgs->Print("v");
  // Construct the desired PDF:
  if (fitFunc.Contains("Bern")) {
    bern = new RooBernsteinM("bkgPdf", "bkgPdf", *m_yy, *bkgArgs, &min, &max);
    background = bern;
  }
  else if (fitFunc.Contains("Exppol")) {
    exppol = new RooGenericPdf("bkgPdf", expFitFormat, *bkgArgs);
    background = exppol;
  }
  
  // Returns the constructed PDF:
  std::cout << "DMBkgModel: Finished constructing " << fitFunc << std::endl;
  return background;
}

/**
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RooRealVar* DMBkgModel::getMassObservable() {
  return m_yy;
}

/**
   Set the pointer to the observable. 
   @param newObservable - The new RooRealVar observable to use for datasets. 
   @returns void.
 */
void DMBkgModel::setMassObservable(RooRealVar *newObservable) {
  m_yy = newObservable;
}

/**
   Get the order of the function from the name.
   @param fitFunc - the name of the function.
   @returns - the order of the function.
*/
int DMBkgModel::getOrderFromFunc(TString fitFunc) {
  for (int i_o = 0; i_o < 10; i_o++) {
    if (fitFunc.Contains(Form("O%d",i_o))) {
      return i_o;
    }
  }
  return 0;
}
