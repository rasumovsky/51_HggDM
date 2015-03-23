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
  jobName = newJobname;
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
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMBkgModel",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
      
  // Either load the masspoints from file or create new ones:
  if (options.Contains("FromFile")) loadMassPointsFromFile();
  else createNewMassPoints();
  return;
}

/**
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RooRealVar* DMMassPoints::getMassObservable() {
  return m_yy;
}

/**
   Returns a pointer to the RooCategory used in the combined dataset.
   @returns pointer to the RooCategory object.
*/
RooRealVar* DMMassPoints::getRooCategory() {
  return categories;
}

/**
   Set the pointer to the observable. 
   @param newObservable - The new RooRealVar observable to use for datasets. 
   @returns void.
 */
void DMMassPoints::setMassObservable(RooRealVar *newObservable) {
  m_yy = newObservable;
}

/**
   Set the pointer to the RooCategory object. 
   @param newCategories - The new RooCategory to use for the combined dataset. 
   @returns void.
 */
void DMMassPoints::setRooCategory(RooCategory *newCategories) {
  categories = newCategories;
}
