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
////////////////////////////////////////////////////////////////////////////////

#include "DMMassPoints.h"

/**
   -----------------------------------------------------------------------------
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
  m_configFileName = newConfigFile;
  m_options = newOptions;
  m_hists.clear();
  
  // Load the config file:
  m_config = new Config(m_configFileName);
  
  // Assign the observable based on inputs:
  if (newObservable == NULL) {
    m_yy = new RooRealVar("m_yy", "m_yy", m_config->getNum("DMMyyRangeLo"),
			  m_config->getNum("DMMyyRangeHi"));
  }
  else {
    setMassObservable(newObservable);
  }
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DMMassPoints", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Check if data should be weighted:
  m_isWeighted = DMAnalysis::isWeightedSample(m_config, m_sampleName);
  
  // Either load the masspoints from file or create new ones:
  if (m_options.Contains("FromFile")) loadMassPointsFromFile();
  else createNewMassPoints();
  
  std::cout << "DMMassPoints: Successfully initialized!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Create local copies of files and a new file list.
   @param originListName - The name of the original file list.
   @returns - The name of the new file list.
*/
TString DMMassPoints::createLocalFilesAndList(TString originListName) {
  int fileIndex = 0;
  TString outputListName = "temporaryList.txt";
  TString currFile;
  ifstream originFile(originListName);
  ofstream outputFile(outputListName);
  while (!originFile.eof()) {
    originFile >> currFile;
    if (currFile.Contains("eos/atlas")) outputFile << currFile << std::endl;
    else {
      system(Form("cp %s file%d.root", currFile.Data(), fileIndex));
      outputFile << Form("file%d.root", fileIndex) << std::endl;
    }
    fileIndex++;
  }
  originFile.close();
  outputFile.close();
  return outputListName;
}

/**
   -----------------------------------------------------------------------------
   Fill a stored 1D histogram.
   @param varName - The name of the quantity in the plot.
   @param allEvents - True iff all events included.
   @param xVal - The x-axis value for histogram filling.
   @param xWeight - The weight value for histogram filling.
   @param cateIndex - The analysis category for the event.
*/
void DMMassPoints::fillHist1D(TString varName, bool allEvents, double xVal,
			      double xWeight, int cateIndex) {
  if (allEvents) {
    m_hists[Form("%s_ALL",varName.Data())]->Fill(xVal, xWeight);
  }
  else {
    m_hists[Form("%s_PASS",varName.Data())]->Fill(xVal, xWeight);
    if (cateIndex >= 0) {
      m_hists[Form("%s_c%d_PASS",varName.Data(),cateIndex)]->Fill(xVal,xWeight);
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Create a RooDataSet containing the mass points in a given category.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The RooDataSet of the data in the specified category.
*/
RooDataSet* DMMassPoints::getCateDataSet(int cateIndex) {
  m_cateData[cateIndex]->Print("v");
  return m_cateData[cateIndex];
}

/**
   -----------------------------------------------------------------------------
   Returns a pointer to the mass observable used in the dataset.
   @returns pointer to the observable (m_yy).
*/
RooRealVar* DMMassPoints::getMassObservable() {
  return m_yy;
}

/**
   -----------------------------------------------------------------------------
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @returns The full path of the mass points text file.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex) {
  return getMassPointsFileName(cateIndex, m_sampleName);
}

/**
   -----------------------------------------------------------------------------
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @param newSampleName - The name of the data/MC sample.
   @returns The full path of the mass points text file.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex, TString sampleName) {
  TString name = Form("%s/%s_%d_%s.txt", m_outputDir.Data(), 
		      (m_config->getStr("cateScheme")).Data(),
		      cateIndex, sampleName.Data());
  return name;
}

/**
   -----------------------------------------------------------------------------
   Merge two DMMassPoint objects, essentially combining the datasets. The 
   global settings will be preserved, except for the sample name, which will
   be reset. Note: Should only be used for datasets with the same number of 
   categories, the same observables, and both using weighted data.
   @param newSampleName - The name of the data/MC sample.
   @param inputMassPoints - The MassPoints object to merge into this one.
   @param saveMassPoints - True iff. you want to immediately save to text file.
                           Preferrable to save until last merge is complete.
*/
void DMMassPoints::mergeMassPoints(TString newSampleName, 
				   DMMassPoints *inputMassPoints,
				   bool saveMassPoints) {

  
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    // First rename the datasets in this object:
    m_cateData[i_c]
      ->SetNameTitle(Form("%s_%s_%d", newSampleName.Data(),
			  (m_config->getStr("cateScheme")).Data(), i_c),
		     Form("%s_%s_%d", newSampleName.Data(),
			  (m_config->getStr("cateScheme")).Data(), i_c));
    
    // Then add the input datasets:
    m_cateData[i_c]->append(*(inputMassPoints->getCateDataSet(i_c)));
    
    // Open input and output files:
    ifstream inputFile1;
    inputFile1.open(getMassPointsFileName(i_c));
    ifstream inputFile2;
    inputFile2.open(inputMassPoints->getMassPointsFileName(i_c));
    
    ofstream outputFile;
    outputFile.open(getMassPointsFileName(i_c, newSampleName));
    
    // Then save the output, if requested:
    if (saveMassPoints) {
      double currMass; double currWeight;
      while (!inputFile1.eof()) {
	inputFile1 >> currMass >> currWeight;
	outputFile << currMass << " " << currWeight << std::endl;
      }
      
      while (!inputFile2.eof()) {
	inputFile2 >> currMass >> currWeight;
	outputFile << currMass << " " << currWeight << std::endl;
      }
    }  
    // Close input and output files:
    inputFile1.close();
    inputFile2.close();
    outputFile.close();
  }
  m_sampleName = newSampleName;
}

/**
   -----------------------------------------------------------------------------
   Create a new 1D histogram.
   @param varName - The name of the quantity in the plot.
   @param nBins - The number of bins.
   @param xMin - The minimum value of the histogram.
   @param xMax - The maximum value of the histogram.
*/
void DMMassPoints::newHist1D(TString varName, int nBins, double xMin,
			     double xMax) {
  // Inclusive histograms:
  m_hists[Form("%s_ALL",varName.Data())]
    = new TH1F(Form("%s_ALL",varName.Data()), Form("%s_ALL",varName.Data()),
	       nBins, xMin, xMax);
  m_hists[Form("%s_PASS",varName.Data())]
    = new TH1F(Form("%s_PASS",varName.Data()), Form("%s_PASS",varName.Data()),
	       nBins, xMin, xMax);
  
  // Categorized histograms:
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    m_hists[Form("%s_c%d_PASS",varName.Data(),i_c)]
      = new TH1F(Form("%s_c%d_PASS",varName.Data(),i_c),
		 Form("%s_c%d_PASS",varName.Data(),i_c),nBins,xMin,xMax);
  }
}

/**
   -----------------------------------------------------------------------------
*/
void DMMassPoints::saveHists() {
  TFile *outputFile = new TFile(Form("%s/hists_%s.root", m_outputDir.Data(),
				     m_sampleName.Data()), "RECREATE");
  // Loop over the saved histograms:
  std::map<TString,TH1F*>::iterator histIter;
  for (histIter = m_hists.begin(); histIter != m_hists.end(); histIter++) {
    histIter->second->Write();
  }
  
  outputFile->Close();
  delete outputFile;
}

/**
   -----------------------------------------------------------------------------
   Set the pointer to the observable. 
   @param newObservable - The new RooRealVar observable to use for datasets. 
   @returns void.
 */
void DMMassPoints::setMassObservable(RooRealVar *newObservable) {
  m_yy = newObservable;
}

/**
   -----------------------------------------------------------------------------
   Create new mass points by looping over the TTree.
   @returns void.
*/
void DMMassPoints::createNewMassPoints() {
  std::cout << "DMMassPoints: creating new mass points from tree." << std::endl;
  
  // Alternative: use file list:
  TString listName = DMAnalysis::nameToFileList(m_config, m_sampleName);
  // If option says copy files, make local file copies and then run:
  if (m_options.Contains("CopyFile")) {
    listName = createLocalFilesAndList(listName);
  }
  TChain *chain = CommonFunc::MakeChain("CollectionTree", listName, "badfile");
  DMTree *dmt = new DMTree(chain);
  
  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect(dmt, m_configFileName);
  
  // Tool to get the total number of events at the generator level.
  DMxAODCutflow *dmx
    = new DMxAODCutflow(DMAnalysis::nameToxAODCutFile(m_config, m_sampleName));
  double nGeneratedEvt = dmx->getEventsPassingCut(1);
  delete dmx;
  
  // Tool to load cross sections and branching ratios:
  BRXSReader *brxs = new BRXSReader(m_configFileName);
  
  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*m_yy);
  
  // Define histograms to save:
  newHist1D("pTyy", 40, 0.0, 600.0);
  newHist1D("ETMiss", 40, 0.0, 600.0);
  newHist1D("ratioETMisspTyy", 40, 0.0, 4.0);
  newHist1D("aTanRatio", 40, 0.0, TMath::Pi()/2.0);
  newHist1D("myy", 40, 105.0, 160.0);
  newHist1D("sumSqrtETMisspTyy", 40, 0.0, 600);
  newHist1D("dPhiyyETMiss", 40, 0.0, TMath::Pi());
  newHist1D("njets", 8, 0, 8);
  newHist1D("nleptons", 5, 0, 5);
  
  // Define datasets and mass files in loop over categories:
  ofstream massFiles[20];
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
      double pileupWeight = 1.0;//dmt->EventInfoAuxDyn_PileupWeight;
      evtWeight *= (m_config->getNum("analysisLuminosity") * 
		    pileupWeight / nGeneratedEvt);
      
      // Multiply by the appropriate luminosity, xsection & branching ratio.
      if (DMAnalysis::isSMSample(m_config, m_sampleName)) {
	evtWeight *= ((brxs->getSMBR(m_config->getNum("higgsMass"),
				     "gammagamma", "BR")) *
		      (brxs->getSMXS(m_config->getNum("higgsMass"), 
				     m_sampleName, "XS")));
      }
      // Dark matter XSBR includes cross-section and branching ratio.
      else if (DMAnalysis::isDMSample(m_config, m_sampleName)) {
	// Multiply by cross-section in pb:
	evtWeight *=brxs->getDMXSBR(DMAnalysis::getMediatorMass(m_config,
								m_sampleName),
				    DMAnalysis::getDarkMatterMass(m_config,
								  m_sampleName),
				    DMAnalysis::getMediatorName(m_sampleName),
				    "XS");
	evtWeight *=brxs->getDMXSBR(DMAnalysis::getMediatorMass(m_config,
								m_sampleName),
				    DMAnalysis::getDarkMatterMass(m_config,
								  m_sampleName),
				    DMAnalysis::getMediatorName(m_sampleName),
				    "FEFF");
      }
      else if (DMAnalysis::isWeightedSample(m_config, m_sampleName)) {
	evtWeight *= brxs->getMCXS(m_sampleName, "XS");
	evtWeight *= brxs->getMCXS(m_sampleName, "FEFF");
      }
      else {
	std::cout << "DMMassPoint: Error! No weighting procedure defined!"
		  << std::endl;
	exit(0);
      }
    }
    
    // The mass parameter:
    double invariantMass = dmt->HGamEventInfoAuxDyn_HighMet_yy_m;
    if (m_config->getBool("RescaleAFII") && 
	DMAnalysis::isDMSample(m_config, m_sampleName)) {
      invariantMass += 1.0;
    }
    // Then commence plotting for ALL events:
    fillHist1D("pTyy", true, dmt->HGamEventInfoAuxDyn_HighMet_yy_pt,
	       evtWeight, -1);
    fillHist1D("ETMiss", true, dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST,
	       evtWeight, -1);
    fillHist1D("ratioETMisspTyy", true, 
	       (dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST/
		dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), 
	       evtWeight, -1);
    fillHist1D("sumSqrtETMisspTyy", true, sqrt(dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST*dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST + dmt->HGamEventInfoAuxDyn_HighMet_yy_pt*dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), 
	       evtWeight, -1);
    fillHist1D("aTanRatio", true, 
	       TMath::ATan(dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST/
			   dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), 
	       evtWeight, -1);
    fillHist1D("myy", true, invariantMass, evtWeight, -1);
    fillHist1D("dPhiyyETMiss", true, 
	       dmt->HGamEventInfoAuxDyn_HighMet_yy_met_deltaPhi, evtWeight, -1);
    fillHist1D("njets", true, dmt->HGamEventInfoAuxDyn_HighMet_jet_n,
	       evtWeight, -1);
    fillHist1D("nleptons", true, dmt->HGamEventInfoAuxDyn_HighMet_lep_n2,
	       evtWeight, -1);
    
    // Make sure events pass the selection:
    if (!selector->passesCut("allCuts", evtWeight)) continue;
    
    // Save the categories:
    int currCate = selector->getCategoryNumber(m_config->getStr("cateScheme"),
					       evtWeight);
    
    // Fill the datasets:
    m_yy->setVal(invariantMass);
    if (m_isWeighted) {
      wt.setVal(evtWeight);
      m_cateData[currCate]->add(RooArgSet(*m_yy,wt), evtWeight);
    }
    else {
      m_cateData[currCate]->add(*m_yy);
    }
    massFiles[currCate] << invariantMass << " " << evtWeight << std::endl;
    
    // Then commence plotting for PASSING events:
    fillHist1D("pTyy", false, dmt->HGamEventInfoAuxDyn_HighMet_yy_pt,
	       evtWeight, currCate);
    fillHist1D("ETMiss", false, dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST,
	       evtWeight, currCate);
    fillHist1D("ratioETMisspTyy", false, 
	       (dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST/
		dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), 
	       evtWeight, currCate);
    fillHist1D("sumSqrtETMisspTyy", false, sqrt(dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST*dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST + dmt->HGamEventInfoAuxDyn_HighMet_yy_pt*dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), evtWeight, currCate);
    fillHist1D("aTanRatio", false, 
	       TMath::ATan(dmt->HGamEventInfoAuxDyn_HighMet_MET_reb_TST/
			   dmt->HGamEventInfoAuxDyn_HighMet_yy_pt), 
	       evtWeight, currCate);
    fillHist1D("myy", false, invariantMass, evtWeight, currCate);
    fillHist1D("dPhiyyETMiss", false, 
	       dmt->HGamEventInfoAuxDyn_HighMet_yy_met_deltaPhi,
	       evtWeight, currCate);
    fillHist1D("njets", false, dmt->HGamEventInfoAuxDyn_HighMet_jet_n,
	       evtWeight, currCate);
    fillHist1D("nleptons", false, dmt->HGamEventInfoAuxDyn_HighMet_lep_n2,
	       evtWeight, -1);
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
  
  // If options said to run locally, remove the files.
  if (m_options.Contains("CopyFile")) {
    removeLocalFilesAndList(listName);
  }

  // Finally, save the histograms to file:
  saveHists();
  
  std::cout << "DMMassPoints: Finished creating new mass points!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
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
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    
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
    while (massFile >> readMass >> readWeight) {
      m_yy->setVal(readMass);
      if (m_isWeighted) {
	wt.setVal(readWeight);
	m_cateData[i_c]->add(RooArgSet(*m_yy, wt), readWeight);
      }
      else {
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

/**
   -----------------------------------------------------------------------------
   Remove the local files and file list.
   @param listName - The name of the local file list.
*/
void DMMassPoints::removeLocalFilesAndList(TString listName) {
  TString currFile;
  ifstream localFile(listName);
  while (!localFile.eof()) {
    localFile >> currFile;
    if (!currFile.Contains("eos/atlas")) system(Form("rm %s", currFile.Data()));
  }
  localFile.close();
  system(Form("rm %s", listName.Data()));
}
