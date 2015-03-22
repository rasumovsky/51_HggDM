////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMassPoints.cxx                                                    //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This class uses a ROOT file with a TTree to implement the analysis event  //
//  categorization and produce text files or RooDataSets.                     //
//                                                                            //
//  Can either run over the TTrees to create new mass points, or load the     //
//  mass points from a previously generated text file, using newOptions =     //
//  "FromFile" or "New".                                                      //
//                                                                            //
//  Still need to get input TTree name convention. Solve in DMHeader file?    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMassPoints.h"

/**
   Initialize the MassPoint class without an observable RooRealVar.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @returns void.
*/
DMMassPoints::DMMassPoints(TString newJobName, TString newSampleName, 
			   TString newCateScheme, TString newOptions) {
  RooRealVar *newObservable = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,
					     DMMyyRangeHi);
  DMMassPoints(newJobName, newSampleName, newCateScheme, newOptions, 
	       newObservable);
}

/**
   Initialize the MassPoint class.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in the datasets (m_yy).
   @returns void.
*/
DMMassPoints::DMMassPoints(TString newJobName, TString newSampleName, 
			   TString newCateScheme, TString newOptions,
			   RooRealVar *newObservable) {
  
  // Load the selector to get category information.
  DMEvtSelect *selector = new DMEvtSelect();
  
  // Define a new RooCategory for the dataset, since none was provided:
  RooCategory newCategories = new RooCategory(Form("categories_%s",
						   newCateScheme.Data()),
					      Form("categories_%s",
						   newCateScheme.Data()));
  // Loop over categories to define categories:
  for (int i_c = 0; i_c < selector->getNCategories(newCateScheme); i_c++) {
    cat->defineType(Form("%s_%d",newCateScheme.Data(),i_c));
    //categories->setRange(Form("rangeName_",i_b,i_r),Form("%s_%d",cateScheme.Data(),i_c));
  }
  
  // Then call the full initializer:
  DMMassPoints(newJobName, newSampleName, newCateScheme, newOptions, 
	       newObservable, newCategories);
}

/**
   Initialize the MassPoint class.
   @param newJobName - The name of the job 
   @param newSampleName - The name of the data/MC sample
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile")
   @param newObservable - The RooRealVar to be used in the datasets (m_yy).
   @param newCategories - The RooCategory to be used in the combined dataset.
   @returns void.
*/
DMMassPoints::DMMassPoints(TString newJobName, TString newSampleName, 
			   TString newCateScheme, TString newOptions,
			   RooRealVar *newObservable,
			   RooCategory *newCategories) {
  std::cout << std::endl << "DMMassPoints::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign the observable and categorization based on inputs:
  setMassObservable(newObservable);
  setRooCategory(newCategories);

  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMMassPoints",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  
  // Check if data should be weighted:
  isWeighted = sampleName.Contains("MC");
  
  // Either load the masspoints from file or create new ones:
  if (options.Contains("FromFile")) loadMassPointsFromFile();
  else createNewMassPoints();
  return;
}

/**
   Create a RooDataSet containing the mass points in a given category.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The RooDataSet of the data in the specified category.
*/
RooDataSet* DMMassPoints::getCateDataSet(int cateIndex) {
  cateData[cateIndex]->Print("v");
  return cateData[cateIndex];
}

/**
   Create a RooDataSet containing the mass points in all categories.
   @returns The combined RooDataSet object for the given categorization.
*/
RooDataSet* DMMassPoints::getCombDataSet() {
  combData->Print("v");
  return combData;
}

/**
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RoORealVar* DMMassPoints::getMassObservable() {
  return m_yy;
}

/**
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The full path of the mass points text file.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex) {
  TString name = Form("%s/%s_%d.txt",outputDir.Data(),cateScheme.Data(),
		      cateIndex);
  return name;
}

/**
   Returns a pointer to the RooCategory used in the combined dataset.
   @returns pointer to the RooCategory object.
*/
RoORealVar* DMMassPoints::getRooCategory() {
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

/**
   Create new mass points by looping over the TTree.
   @returns void.
*/
void DMMassPoints::createNewMassPoints() {
  
  std::cout << "DMMassPoints: creating new mass points from tree." << std::endl;

  // Load the input TFile and TTree here:
  TFile *myFile = new TFile(nameToSample[sampleName]);
  TTree *myTree = (TTree*)myFile->Get("treename");
  DMTree *dmt = new DMTree(myTree);

  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect(dmt);
  
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*m_yy);
  
  ofstream massFiles[20];
  
  // Loop over categories to define datasets and mass files:
  for (int i_c = 0; i_c < selector->getNCategories(cateScheme); i_c++) {
    
    massFiles[i_c].open(getMassPointsFileName(i_c));
    
    if (isWeighted) {
      args->add(wt);
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     RooArgSet(*m_yy,wt), 
				     RooFit::WeightVar(wt)); 
    }
    else {
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     *m_yy);
    }
    dataMap[Form("%s_%d",cateScheme.Data(),i_c)] = cateData[i_c];
  }
    
  // Loop over the input DMTree:
  Long64_t entries = dmt->fChain->GetEntries();
  for (Long64_t event = 0; event < entries; event++) {
    dmt->fChain->GetEntry(event);
    
    // Check the cutflow:
    if (!selector->passesCut("all")) continue;
    
    
    
    // Save the categories:
    int currCate = selector->getCategoryNumber(cateScheme);
    if (currCate > 0) {
      
      m_yy->setVal(dmt->EventInfoAuxDyn_m_yy);
      
      if (isWeighted) {
	// Get the event weight:
	double evtWeight = dmt->EventInfoAuxDyn_PileupWeight;
	
	wt.setVal(evtWeight);
	cateData[currCate]->add(RooArgSet(*m_yy,wt), evtWeight);
	massFiles[currCate] << dmt->EventInfoAuxDyn_m_yy << " " << evtWeight
			    << std::endl;
      }
      else {
	cateData[currCate]->add(*m_yy);
	massFiles[currCate] << dmt->EventInfoAuxDyn_m_yy << std::endl;
      }
    }
  }
  
  // Print the cutflow and category yields (weighted if MC):
  selector->printCutflow(isWeighted);
  selector->printCategorization(isWeighted);
  selector->saveCutflow(Form("%s/cutflow_%s.txt",outputDir.Data(),
			     sampleName.Data()), isWeighted);
  selector->saveCategorization(Form("%s/categorization_%s_%s.txt",
				    outputDir.Data(),
				    cateScheme.Data(),
				    sampleName.Data()), isWeighted);
  
  // Close output mass point files.
  for (int i_f = 0; i_f < selector->getNCategories(cateScheme); i_f++) {
    massFiles[i_f].close();
  }
  
  // Create combined data set from individual categories:
  combData = new RooDataSet(Form("combData_%s",cateScheme.Data()),
			    Form("combData_%s",cateScheme.Data()),
			    *args, RooFit::Index(*categories),
			    RooFit::Import(dataMap),
			    RooFit::WeightVar(wt));
}
  
/**
   Load the mass points from text files that have already been produced. This is
   much faster than producing mass points from scratch, and is preferred. 
   @returns void.
*/
void DMMassPoints::loadMassPointsFromFile() {
  
  std::cout << "DMMassPoints: loading mass points from .txt file." << std::endl;
  
  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect();
  
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  
  RooArgSet *args = new RooArgSet();
  RooRealVar wt("wt","wt",1);
  args->add(*m_yy);
  
  // Loop over categories to define datasets and mass files:
  for (int i_c = 0; i_c < selector->getNCategories(cateScheme); i_c++) {
    
    if (isWeighted) {
      args->add(wt);
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     RooArgSet(*m_yy,wt), 
				     RooFit::WeightVar(wt)); 
    }
    else {
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     *m_yy);
    }
        
    double readMass; double readWeight;
    ifstream massFile(getMassPointsFileName(i_c));
    if (isWeighted) {
      while (massFile >> readMass >> readWeight) {
	wt.setVal(readWeight);
	m_yy->setVal(readMass);
	cateData[i_c]->add(RooArgSet(*m_yy,wt), readWeight);
      }
    }
    else {
      while (massFile >> readMass) {
	m_yy->setVal(readMass);
	cateData[i_c]->add(*m_yy);
      }
    }
    // Add the category dataset to the dataset map:
    dataMap[Form("%s_%d",cateScheme.Data(),i_c)] = cateData[i_c];
  }
    
  // Create combined data set from individual categories:
  combData = new RooDataSet(Form("combData_%s",cateScheme.Data()),
			    Form("combData_%s",cateScheme.Data()),
			    *args, RooFit::Index(*categories),
			    RooFit::Import(dataMap), RooFit::WeightVar(wt));
}
