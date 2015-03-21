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
   Initialize the MassPoint class.
*/
DMMassPoints::DMMassPoints(TString newJobName, TString newSampleName, 
			   TString newCateScheme, TString newOptions) {
  std::cout << std::endl << "DMMassPoints::Initializing..." << std::endl;
  
  // Assign member variables:
  jobName = newJobName;
  sampleName = newSampleName;
  cateScheme = newCateScheme;
  options = newOptions;
  
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
   Get the name of the output textfile for the given category index.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex) {
  TString name = Form("%s/%s_%d.txt",outputDir.Data(),cateScheme.Data(),
		      cateIndex);
  return name;
}

/**
   Create a RooDataSet containing the mass points in a given category.
*/
RooDataSet* DMMassPoints::getCateDataSet(int cateIndex) {
  cateData[cateIndex]->Print("v");
  return cateData[cateIndex];
}

/**
   Create a RooDataSet containing the mass points in all categories.
*/
RooDataSet* DMMassPoints::getCombDataSet() {
  combData->Print("v");
  return combData;
}
   
/**
   Create new masspoints by looping over the TTree.
*/
void DMMassPoints::createNewMassPoints() {
  
  std::cout << "DMMassPoints: creating new mass points from tree." << std::endl;

  // Load the input TFile and TTree here:
  TFile *myFile = new TFile(prodToSample[sampleName]);
  TTree *myTree = (TTree*)myFile->Get("treename");
  DMTree *dmt = new DMTree(myTree);

  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect(dmt);
  
  // RooFit stuff:
  RooCategory *Categories = new RooCategory(Form("Categories_%s",cateScheme.Data()),
					  Form("Categories_%s",cateScheme.Data())
					  );
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  RooRealVar *m_yy = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,DMMyyRangeHi);
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
				     RooArgSet(*m_yy,wt), WeightVar(wt)); 
    }
    else {
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),i_c),
				     *m_yy);
    }
    
    Categories->defineType(Form("%s_%d",cateScheme.Data(),i_c));
    //data_category->setRange(Form("rangeName_",i_b,i_r),Form("%s_%d",cateScheme.Data(),i_c));
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
      
      m_yy->setVal(dmt->EventInfoAuxDyn.m_yy);
      if (isWeighted) {
	wt.setVal(dmt->EventInfoAuxDyn.weight);
	cateData[currCate]->add(RooArgSet(*m_yy,wt), readWeight);
	massFiles[currCate] << dmt->EventInfoAuxDyn.m_yy << " "
			    << dmt->EventInfoAuxDyn.weight << std::endl;
      }
      else {
	cateData[currCate]->add(*m_yy);
	massFiles[currCate] << dmt->EventInfoAuxDyn.m_yy << std::endl;
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
			    *args, Index(*Categories),Import(dataMap),
			    RooFit::WeightVar(wt));
}
  
/**
   Load the mass points from text files that have already been produced.
*/
void DMMassPoints::loadMassPointsFromFile() {
  
  std::cout << "DMMassPoints: loading mass points from .txt file." << std::endl;
  
  // RooFit stuff:
  RooCategory *Categories = new RooCategory(Form("Categories_%s",
						 cateScheme.Data()),
					    Form("Categories_%s",
						 cateScheme.Data())
					    );
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  RooRealVar *m_yy = new RooRealVar("m_yy","m_yy",DMMyyRangeLo,DMMyyRangeHi);
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*m_yy);
  
  // Loop over categories to define datasets and mass files:
  for (int i_c = 0; i_c < selector->getNCategories(cateScheme); i_c++) {
    
    if (weighted) {
      args->add(wt);
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),cateIndex),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),cateIndex),
				     RooArgSet(*m_yy,wt), WeightVar(wt)); 
    }
    else {
      cateData[i_c] = new RooDataSet(Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),cateIndex),
				     Form("%s_%s_%d",sampleName.Data(),
					  cateScheme.Data(),cateIndex),
				     *m_yy);
    }
    
    data_category->defineType(Form("%s_%d",cateScheme.Data(),i_c));
    //data_category->setRange(Form("rangeName_",i_b,i_r),Form("%s_%d",cateScheme.Data(),i_c));
    
    double readMass; double readWeight;
    ifstream massFile(getMassPointsFileName(i_c));
    if (weighted) {
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
			    *args, Index(*Categories),Import(dataMap),
			    WeightVar(wt));
}
