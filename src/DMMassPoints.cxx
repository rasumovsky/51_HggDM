////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMassPoints.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
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
//  Currently, gg_gjet sample has a 'loose selection' applied. Also, the      //
//  normalization should be hard-coded as an extrapolation of the 8TeV        //
//  analysis background.                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMassPoints.h"

using namespace DMAnalysis;

/**
   Initialize the DMMassPoint class and make a new RooCategory.
   @param newConfigFile - The name of the config file.
   @param newSampleName - The name of the data/MC sample.
   @param newOptions - The job options ("New", "FromFile").
   @param newObservable - The RooRealVar to be used in the datasets (m_yy).
   @returns void.
*/
DMMassPoints::DMMassPoints(TString newConfigFile, TString newSampleName, 
			   TString newOptions, RooRealVar *newObservable) {
  std::cout << "\nDMMassPoints::Initializing..." 
	    << "\n\tconfigFile = " << newConfigFile
	    << "\n\tsampleName = " << newSampleName 
	    << "\n\toptions = " << newOptions << std::endl;
  
  // Assign member variables:
  m_sampleName = newSampleName;
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  
  // Assign the observable based on inputs:
  if (newObservable == NULL) {
    m_yy = new RooRealVar("m_yy", "m_yy", DMMyyRangeLo, DMMyyRangeHi);
  }
  else {
    setMassObservable(newObservable);
  }
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DMMassPoints", masterOutput.Data(), 
		   (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Check if data should be weighted:
  m_isWeighted = isWeightedSample(m_sampleName);
  
  // Either load the masspoints from file or create new ones:
  if (newOptions.Contains("FromFile")) loadMassPointsFromFile();
  else createNewMassPoints();
  
  std::cout << "DMMassPoints: Successfully initialized!" << std::endl;
  return;
}

/**
   Create a RooDataSet containing the mass points in a given category.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The RooDataSet of the data in the specified category.
*/
RooDataSet* DMMassPoints::getCateDataSet(int cateIndex) {
  m_cateData[cateIndex]->Print("v");
  return m_cateData[cateIndex];
}

/**
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RooRealVar* DMMassPoints::getMassObservable() {
  return m_yy;
}

/**
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The full path of the mass points text file.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex) {
  TString name = Form("%s/%s_%d_%s.txt", m_outputDir.Data(), 
		      (m_config->getStr("cateScheme")).Data(),
		      cateIndex, m_sampleName.Data());
  return name;
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
   Create new mass points by looping over the TTree.
   @returns void.
*/
void DMMassPoints::createNewMassPoints() {
  std::cout << "DMMassPoints: creating new mass points from tree." << std::endl;
  
  // Alternative: use file list:
  TString listName = nameToFileList(m_sampleName);
  TChain *chain = CommonFunc::MakeChain("CollectionTree", listName, "badfile");
  DMTree *dmt = new DMTree(chain);
  
  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect(dmt);
  
  // Tool to get the total number of events at the generator level.
  DMxAODCutflow *dmx
    = new DMxAODCutflow(DMAnalysis::nameToxAODCutFile(m_sampleName));
  double nGeneratedEvt = dmx->getEventsPassingCut(1);
  
  // Tool to load cross sections and branching ratios:
  BRXSReader *brxs = new BRXSReader(Form("%s/XSBRInputs/",masterInput.Data()));
  
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*m_yy);
  
  ofstream massFiles[20];
  
  // Loop over categories to define datasets and mass files:
  std::cout << "DMMassPoints: Define datasets & files." << std::endl;
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    
    massFiles[i_c].open(getMassPointsFileName(i_c));
    
    if (m_isWeighted) {
      args->add(wt);
      m_cateData[i_c]
	= new RooDataSet(Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 RooArgSet(*m_yy, wt), RooFit::WeightVar(wt)); 
    }
    else {
      m_cateData[i_c]
	= new RooDataSet(Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 *m_yy);
    }
    dataMap[Form("%s_%d",(m_config->getStr("cateScheme")).Data(),i_c)]
      = m_cateData[i_c];
  }
    
  // Loop over the input DMTree:
  Long64_t entries = dmt->fChain->GetEntries();
  std::cout << "DMMassPoints: Loop over DMTree with " << entries
	    << " entries." << std::endl;
  for (Long64_t event = 0; event < entries; event++) {
    dmt->fChain->GetEntry(event);
        
    // Calculate the weights for the cutflow first!
    double evtWeight = 1.0;
    if (m_isWeighted) {
      evtWeight *= (DMAnalysis::analysisLuminosity * 
		    dmt->EventInfoAuxDyn_PileupWeight / nGeneratedEvt);
      
      // Multiply by the appropriate luminosity, xsection & branching ratio.
      if (isSMSample(m_sampleName)) {
	evtWeight *= ((brxs->getSMBR(higgsMass, "gammagamma", "BR")) *
		      (brxs->getSMXS(higgsMass, m_sampleName, "XS")));
      }
      // Dark matter XSBR includes cross-section and branching ratio.
      // WARNING!!! GETTING RID OF XSBR
      else if (isDMSample(m_sampleName)) {
	/*
	evtWeight *= 10000 * ((brxs->getDMXSBR(getMediatorMass(m_sampleName),
					       getDarkMatterMass(m_sampleName),
					       getMediatorName(m_sampleName),
					       "XS")));
	*/
	evtWeight *= 0.006566;// cross-section for mX=1GeV,mZ'=100GeV
      }
      else if (m_sampleName.EqualTo("gg_gjet")) {
	evtWeight *= 57.24;//xsection*filter-eff for Sherpa gg+gj
	evtWeight *= 0.2729;// To scale to the expected bkg per fb-1
      }
    }
    
    // Check the cutflow (loose for background sample):
    if (m_sampleName.EqualTo("gg_gjet") && 
	!selector->passesCut("looseCuts",evtWeight)) {
      continue;
    }
    else if (!m_sampleName.EqualTo("gg_gjet") && 
	     !selector->passesCut("allCuts",evtWeight)) {
      continue;
    }
    
    // Save the categories:
    int currCate = selector->getCategoryNumber(m_config->getStr("cateScheme"),
					       evtWeight);
    if (currCate >= 0) {
      
      m_yy->setVal(dmt->EventInfoAuxDyn_m_yy);
      
      if (m_isWeighted) {
	wt.setVal(evtWeight);
	m_cateData[currCate]->add(RooArgSet(*m_yy,wt), evtWeight);
	massFiles[currCate] << dmt->EventInfoAuxDyn_m_yy << " " << evtWeight
			    << std::endl;
      }
      else {
	m_cateData[currCate]->add(*m_yy);
	massFiles[currCate] << dmt->EventInfoAuxDyn_m_yy << std::endl;
      }
    }
  }
  std::cout << "DMMassPoints: End of loop over input DMTree." << std::endl;
  
  // Print the cutflow and category yields (weighted if MC):
  selector->printCutflow(m_isWeighted);
  selector->printCategorization(m_isWeighted);
  selector->saveCutflow(Form("%s/cutflow_%s.txt", m_outputDir.Data(),
			     m_sampleName.Data()), m_isWeighted);
  selector->saveCategorization(Form("%s/categorization_%s_%s.txt",
				    m_outputDir.Data(),
				    (m_config->getStr("cateScheme")).Data(),
				    m_sampleName.Data()), m_isWeighted);
  
  // Close output mass point files.
  for (int i_f = 0; i_f < m_config->getInt("nCategories"); i_f++) {
    massFiles[i_f].close();
  }
    
  std::cout << "DMMassPoints: Finished creating new mass points!" << std::endl;
}
  
/**
   Load the mass points from text files that have already been produced. This is
   much faster than producing mass points from scratch, and is preferred. 
   @returns void.
*/
void DMMassPoints::loadMassPointsFromFile() {
  std::cout << "DMMassPoints: loading mass points from .txt file." << std::endl;
  if (m_isWeighted) std::cout << "\tMass points will be weighted." << std::endl;
  else std::cout << "\tMass points will be un-weighted." << std::endl;

  std::map<string,RooDataSet*> dataMap; dataMap.clear();
  RooArgSet *args = new RooArgSet();
  RooRealVar wt("wt","wt",1);
  args->add(*m_yy);
  
  // Loop over categories to define datasets and mass files:
  for (int i_c = 0; i_c < m_config->getInt("nCategories"))); i_c++) {
    
    if (m_isWeighted) {
      args->add(wt);
      m_cateData[i_c]
	= new RooDataSet(Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 RooArgSet(*m_yy,wt), RooFit::WeightVar(wt)); 
    }
    else {
      m_cateData[i_c]
	= new RooDataSet(Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 Form("%s_%s_%d", m_sampleName.Data(),
			      (m_config->getStr("cateScheme")).Data(), i_c),
			 *m_yy);
    }
        
    double readMass; double readWeight;
    ifstream massFile(getMassPointsFileName(i_c));
    std::cout << "DMMassPoints: opening " << getMassPointsFileName(i_c)
	      << std::endl;
    
    // First check that file exists. If it does not, we need to create inputs.
    if (!massFile) {
      std::cout << "DMMassPoints: Error! Cannot load from file." << std::endl;
      createNewMassPoints();
      return;
    }
    if (m_isWeighted) {
      while (massFile >> readMass >> readWeight) {
	wt.setVal(readWeight);
	m_yy->setVal(readMass);
	m_cateData[i_c]->add(RooArgSet(*m_yy, wt), readWeight);
      }
    }
    else {
      while (massFile >> readMass) {
	m_yy->setVal(readMass);
	m_cateData[i_c]->add(*m_yy);
      }
    }
    // Add the category dataset to the dataset map:
    dataMap[Form("%s_%d", (m_config->getStr("cateScheme")).Data(),i_c)]
      = m_cateData[i_c];
    std::cout << "\tEntries[" << i_c << "] = " << m_cateData[i_c]->sumEntries()
	      << std::endl;
  }
  std::cout << "DMMassPoints: Finished loading data set. " << std::endl;
}
