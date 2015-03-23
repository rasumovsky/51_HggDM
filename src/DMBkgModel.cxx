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
   Initialize the DMBkgModel class and make a new observable RooRealVar and
   RooCategory.
   classes (instead of importing them).
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @returns void.
*/
DMBkgModel::DMBkgModel(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions) {
  RooRealVar *newObservable = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,
					     DMMyyRangeHi);
  DMBkgModel(newJobName, newSampleName, newCateScheme, newOptions, 
	     newObservable);
}

/**
   Initialize the DMBkgModel class and make a new RooCategory.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @returns void.
*/
DMBkgModel::DMBkgModel(TString newJobName, TString newSampleName, 
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
  DMBkgModel(newJobName, newSampleName, newCateScheme, newOptions,
	     newObservable, newCategories);
}

/**
   Initialize the DMBkgModel class using previously defined observable 
   RooRealVar and RooCategory classes.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in fits (m_yy).
   @param newCategories - The RooCategory to be used in the combined PDF.
   @returns void.
*/
DMBkgModel::DMBkgModel(TString newJobName, TString newSampleName, 
		       TString newCateScheme, TString newOptions,
		       RooRealVar *newObservable, RooCategory *newCategories) {
  std::cout << std::endl << "DMBkgModel::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign the observable and categorization based on inputs:
  m_yy = newObservable;
  categories = newCategories;
  
  // Get the number of analysis categories if not already done:
  if (!selector) {
    selector = new DMEvtSelect();
    nCategories = selector->getNCategories(newCateScheme);
  }
  return;
}

/**
   Get the PDF for a given class.
   @param cateIndex - The category of cateScheme used in class initialization.
   @returns The corresponding background PDF for the analysis.
*/
RooAbsPdf* DMBkgModel::getCateBkgPDF(int cateIndex) {
  TString cateName = Form("%s_%d",cateScheme.Data(),cateIndex);
  TString currFunction = cateToBkgFunc[cateName];
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
  
  int order = getOrderFromFunc(fitFunc);
    
  // Set the range of the m_yy variable from DMHeader.h:
  RooConstVar min("min","min",DMMyyRangeLo);
  RooConstVar max("max","max",DMMyyRangeHi);
  
  // Background fit variables:
  RooArgList *bkgArgs = new RooArgList();
  RooRealVar *pVar[10];
  RooRealVar *cVar[10];
  TString expFitFormat = "TMath::Exp(";
  
  // Pointers to the background PDFs:
  RooAbsPdf *background;
  RooBernsteinM *bern;
  RooGenericPdf *exppol;
  
  // Loop over the order of the function to define parameters:
  for (int i_p = 0; i_p <= order; i_p++) {
    // Parameters for Bernstein polynomial:
    if (fitFunc.Contains("Bern")) {
      if (i_p == 0) {
	pVar[i_p] = new RooRealVar(Form("pVar%d",order), 
				   Form("pVar%d",order), 1);
      }
      else {
	pVar[i_p] = new RooRealVar(Form("pVar%d",order), Form("pVar%d",order),
				   0.1, 0.0, 10.0);
      }
      bkgArgs->add(*pVar[i_p]);
    }
    // Parameters for exponential polynomial:
    else if (fitFunc.Contains("Exppol") && i_p < order) {
      cVar[i_p] = new RooRealVar(Form("cVar%d",order), Form("cVar%d",order),
				 0.0, -1.0, 1.0 );
      bkgArgs->add(*cVar[i_p]);
      expFitFormat += Form("@%d",i_p+1);
      for (int i_t = 0; i_t <= i_p; i_t++) {
	expFitFormat += "*(@0-100)";
      }
    }
  }
  expFitFormat += ")";
  
  // Construct the desired PDF:
  if (fitFunc.Contains("Bern")) {
    bern = new RooBernsteinM(fitName, fitName, *m_yy, *bkgArgs, &min, &max);
    background = bern;
  }
  else if (fitFunc.Contains("Exppol")) {
    //bkgArgs has m_yy included as first variable in list:
    //exppol = new RooGenericPdf(fitName, expFitFormat, bkgArgs);
    exppol = new RooGenericPdf(fitName, expFitFormat, *m_yy, bkgArgs);
    background = exppol;
  }
  
  // Returns the constructed PDF:
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
   Returns a pointer to the RooCategory used in the combined dataset.
   @returns pointer to the RooCategory object.
*/
RooRealVar* DMBkgModel::getRooCategory() {
  return categories;
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
   Set the pointer to the RooCategory object. 
   @param newCategories - The new RooCategory to use for the combined dataset. 
   @returns void.
 */
void DMBkgModel::setRooCategory(RooCategory *newCategories) {
  categories = newCategories;
}

int DMBkgModel::getOrderFromFunc(TString fitFunc) {
  for (int i_o = 0; i_o < 10; i_o++) {
    if (fitFunc.Contains(Form("O%d",i_o))) {
      return i_o;
    }
  }
  return 0;
}
