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
//  New option: "Syst" to implement systematic variations.                    //
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
   @return - void.
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
  
  // For the full cut-flow:
  m_componentCutFlows.clear();
  m_componentNorms.clear();
  m_cutFlowHist = NULL;
  
  // Load the config file:
  m_config = new Config(m_configFileName);
  
  // Assign the observable based on inputs:
  if (newObservable == NULL) {
    m_yy = new RooRealVar("m_yy", "m_yy", m_config->getNum("DMMyyRangeLo"),
			  m_config->getNum("DMMyyRangeHi"));
  }
  else setMassObservable(newObservable);
  
  // Assign output directory, and make sure it exists:
  m_outputDir = Form("%s/%s/DMMassPoints", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/Systematics", m_outputDir.Data()));
  
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
   Combine the MxAOD cutflow histograms with the histograms that have analysis-
   specific cut information implemented.
*/
void DMMassPoints::combineCutFlowHists() {
  std::cout << "DMMassPoints: Combine the cutflow histograms." << std::endl;
  TH1F *fullCutFlowHist = new TH1F(Form("cutFlowFull_%s", m_sampleName.Data()),
				   Form("cutFlowFull_%s", m_sampleName.Data()),
				   m_cutFlowHist->GetNbinsX(), 0,
				   m_cutFlowHist->GetNbinsX());
  
  int nMxAODCuts = (int)((m_config->getStrV("MxAODCutList")).size());
  
  // Check the sizes of samples:
  if (m_componentNorms.size() != m_componentCutFlows.size()) {
    std::cout << "DMMassPoints: ERROR! Wrong normalization of samples." 
	      << std::endl;
    exit(0);
  }
  
  // Loop over component cutflows and add them to the new cutflow:
  for (int i_c = 0; i_c < (int)m_componentCutFlows.size(); i_c++) {
    for (int i_b = 1; i_b <= m_componentCutFlows[i_c]->GetNbinsX(); i_b++) {
      double addedValue = (m_componentCutFlows[i_c]->GetBinContent(i_b) * 
			   m_componentNorms[i_c]);
      double addedError = (m_componentCutFlows[i_c]->GetBinError(i_b) * 
			   m_componentNorms[i_c]);
      // Error = sqrt(prevErr^2 + newErr^2)
      double newError = sqrt((fullCutFlowHist->GetBinError(i_b) * 
			      fullCutFlowHist->GetBinError(i_b)) + 
			     (addedError * addedError));
      fullCutFlowHist->AddBinContent(i_b, addedValue);
      fullCutFlowHist->SetBinError(i_b, newError);
    }
  }
  
  // Update the bin labeling:
  for (int i_b = 1; i_b <= fullCutFlowHist->GetNbinsX(); i_b++) {
    fullCutFlowHist->GetXaxis()
      ->SetBinLabel(i_b,(m_componentCutFlows[0])->GetXaxis()->GetBinLabel(i_b));
  }
  
  // Check the results of the cut-flow merging:
  std::cout << "DMMassPoints: The MxAOD has " 
	    << fullCutFlowHist->GetBinContent(nMxAODCuts) 
	    << " while the current program has " 
	    << m_cutFlowHist->GetBinContent(nMxAODCuts) << " for " 
	    << fullCutFlowHist->GetXaxis()->GetBinLabel(nMxAODCuts) 
	    << std::endl;
  double discrepancy = fabs((fullCutFlowHist->GetBinContent(nMxAODCuts) - 
			     m_cutFlowHist->GetBinContent(nMxAODCuts)) / 
			    fullCutFlowHist->GetBinContent(nMxAODCuts));
  if (discrepancy > 0.01) {
    std::cout << "DMMassPoints: ERROR! That discrepancy of " << discrepancy 
	      << " is too large :(" << std::endl;
    exit(0);
  }
  
  // Then add the cuts from this program:
  for (int i_b = nMxAODCuts+1; i_b <= m_cutFlowHist->GetNbinsX(); i_b++) {
    fullCutFlowHist->SetBinContent(i_b, m_cutFlowHist->GetBinContent(i_b));
    fullCutFlowHist->SetBinError(i_b, m_cutFlowHist->GetBinError(i_b));
  }
  m_cutFlowHist = fullCutFlowHist;
  fullCutFlowHist->GetXaxis()->SetBinLabel(nMxAODCuts+1, "Lepton Veto");
  fullCutFlowHist->GetXaxis()
    ->SetBinLabel(nMxAODCuts+2,Form("p_{T}^{#gamma#gamma} > %d GeV",
				    (int)(m_config->getNum("AnaCutDiphotonPT")/
					  1000.0)));
  fullCutFlowHist->GetXaxis()
    ->SetBinLabel(nMxAODCuts+3,Form("#slash{E}_{T} > %d GeV", 
				    (int)(m_config->getNum("AnaCutETMiss")/
					  1000.0)));
  
  m_cutFlowHist->GetXaxis()->SetNdivisions(m_cutFlowHist->GetNbinsX());
  m_cutFlowHist->GetXaxis()->CenterLabels();
  // This TH1F is written to file in the saveHists() method.
}

/**
   -----------------------------------------------------------------------------
   Create local copies of files and a new file list.
   @param originListName - The name of the original file list.
   @return - The name of the new file list.
*/
TString DMMassPoints::createLocalFilesAndList(TString originListName) {
  std::cout << "DMMassPoints: createLocalFilesAndList." << std::endl;
  TString outputListName = "temporaryList.txt";
  TString currLine;
  std::ifstream originFile(originListName);
  std::ofstream outputFile(outputListName);
  while (!originFile.eof()) {
    originFile >> currLine;
    if (currLine.Contains("eos/atlas")) {
      TString newLine = currLine;
      TString directory = currLine.Contains("data") ? 
	Form("%s/", (m_config->getStr("MxAODDirectoryData")).Data()) :
	Form("%s/", (m_config->getStr("MxAODDirectoryMC")).Data());
      newLine.ReplaceAll(directory, "");
      system(Form("xrdcp %s %s", currLine.Data(), newLine.Data()));
      outputFile << newLine << std::endl;
    }
    else {
      outputFile << currLine << std::endl;
    }
  }
  originFile.close();
  outputFile.close();
  std::cout << "DMMassPoints: Finished copying files." << std::endl;
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
   @return - The RooDataSet of the data in the specified category.
*/
RooDataSet* DMMassPoints::getCateDataSet(int cateIndex) {
  m_cateData[cateIndex]->Print("v");
  return m_cateData[cateIndex];
}

/**
   -----------------------------------------------------------------------------
   Returns a pointer to the mass observable used in the dataset.
   @return - A pointer to the observable (m_yy).
*/
RooRealVar* DMMassPoints::getMassObservable() {
  return m_yy;
}

/**
   -----------------------------------------------------------------------------
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @return - The full path of the mass points text file.
*/
TString DMMassPoints::getMassPointsFileName(int cateIndex) {
  return getMassPointsFileName(cateIndex, m_sampleName);
}

/**
   -----------------------------------------------------------------------------
   Get the name of the output textfile for the given category index.
   @param cateIndex - The index of the category for which we want the .txt name.
   @param newSampleName - The name of the data/MC sample.
   @return - The full path of the mass points text file.
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
    std::ifstream inputFile1;
    inputFile1.open(getMassPointsFileName(i_c));
    std::ifstream inputFile2;
    inputFile2.open(inputMassPoints->getMassPointsFileName(i_c));
    
    std::ofstream outputFile;
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
  if (m_isWeighted) m_hists[Form("%s_ALL",varName.Data())]->Sumw2(true);
  
  m_hists[Form("%s_PASS",varName.Data())]
    = new TH1F(Form("%s_PASS",varName.Data()), Form("%s_PASS",varName.Data()),
	       nBins, xMin, xMax);
  if (m_isWeighted) m_hists[Form("%s_PASS",varName.Data())]->Sumw2(true);
  
  // Categorized histograms:
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    m_hists[Form("%s_c%d_PASS",varName.Data(),i_c)]
      = new TH1F(Form("%s_c%d_PASS",varName.Data(),i_c),
		 Form("%s_c%d_PASS",varName.Data(),i_c),nBins,xMin,xMax);
    if (m_isWeighted) {
      m_hists[Form("%s_c%d_PASS",varName.Data(),i_c)]->Sumw2(true);
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Save the histograms, including cutflow and kinematic distributions, to file.
*/
void DMMassPoints::saveHists() {
  TFile *outputFile = new TFile(Form("%s/hists_%s.root", m_outputDir.Data(),
				     m_sampleName.Data()), "RECREATE");
  // Loop over the saved histograms:
  std::map<TString,TH1F*>::iterator histIter;
  for (histIter = m_hists.begin(); histIter != m_hists.end(); histIter++) {
    histIter->second->Write();
  }
  // Also save the cutflow hist to this file:
  m_cutFlowHist->Write();

  outputFile->Close();
  delete outputFile;
}

/**
   -----------------------------------------------------------------------------
   Set the pointer to the observable. 
   @param newObservable - The new RooRealVar observable to use for datasets. 
   @return - void.
 */
void DMMassPoints::setMassObservable(RooRealVar *newObservable) {
  m_yy = newObservable;
}

/**
   -----------------------------------------------------------------------------
   Create new mass points by looping over the TTree.
   @return - void.
*/
void DMMassPoints::createNewMassPoints() {
  std::cout << "DMMassPoints: creating new mass points from tree." << std::endl;
  
  // Check whether experimental systematic uncertainties should be evaluated:
  bool getSystematics = m_options.Contains("Syst");
  
  // Construct file list for the TChain:
  TString listName
    = DMAnalysis::nameToFileList(m_config, m_sampleName, getSystematics);
  
  // If option says copy files, make local file copies and then run:
  if (m_options.Contains("CopyFile")) {
    listName = createLocalFilesAndList(listName);
  }
  TChain *chain = CommonFunc::MakeChain("CollectionTree", listName, "badfile");
  DMTree *dmt = new DMTree(chain);
  
  // Tool to implement the cutflow, categorization, and counting. 
  DMEvtSelect *selector = new DMEvtSelect(dmt, m_configFileName);
  
  // Also instantiate tools with systematics:
  std::map<TString, DMEvtSelect*> sysSelectors; sysSelectors.clear();
  std::vector<TString> systList = m_config->getStrV("SystematicsList");
  if (getSystematics) {
    for (int i_s = 0; i_s < (int)systList.size(); i_s++) {
      sysSelectors[systList[i_s]] = new DMEvtSelect(dmt, m_configFileName);
      sysSelectors[systList[i_s]]->setSysVariation(systList[i_s]);
    }
  }
  
  // For updating the nTotalEventsInFile:
  TString currFileName = "";
  double nTotalEventsInFile = 0.0;

  std::map<string,RooDataSet*> dataMap;
  dataMap.clear();
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*m_yy);
  
  // Define histograms to save:
  newHist1D("pTyy", 25, 0.0, 500.0);
  newHist1D("ETMiss", 20, 0.0, 400.0);
  newHist1D("ratioETMisspTyy", 20, 0.0, 4.0);
  newHist1D("aTanRatio", 20, 0.0, TMath::Pi()/2.0);
  newHist1D("myy", 20, 105.0, 160.0);
  newHist1D("sumSqrtETMisspTyy", 20, 0.0, 600);
  //newHist1D("dPhiyyETMiss", 20, 0.0, TMath::Pi());
  newHist1D("njets", 8, 0, 8);
  newHist1D("nleptons", 5, 0, 5);
  
  // Define datasets and mass files in loop over categories:
  std::ofstream massFiles[20];
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
    
    // Load event and print the progress bar:
    dmt->fChain->GetEntry(event);
    //printProgressBar(event, entries);
    
    // Check if this is a new file (which requires a different overall norm:
    if (!currFileName.EqualTo(chain->GetFile()->GetName())) {
      currFileName = chain->GetFile()->GetName();
      std::cout << "DMMassPoints: Switch to file : " << currFileName
		<< std::endl;
      // Tool to get the total number of events at the generator level.
      DMxAODCutflow *dmx = new DMxAODCutflow(currFileName, m_configFileName);
      nTotalEventsInFile = dmx->nTotalEventsInFile();
      
      // Then also get the cutflow histogram:
      double currNorm = 1.00000000;
      TH1F* currHist = (TH1F*)dmx->getHist()->Clone();
      if (m_isWeighted) {
	// Normalization of inputs to 1 fb-1:
	currNorm
	  = (1000.0 * 	   
	     dmt->HGamEventInfoAuxDyn_crossSectionBRfilterEff["Nominal"] /
	     nTotalEventsInFile);
	
	// Add branching ratio for DM samples:
	if (DMAnalysis::isDMSample(m_config, m_sampleName)) {
	  currNorm *= m_config->getNum("BranchingRatioHyy");
	}
      }
      m_componentCutFlows.push_back(currHist);
      m_componentNorms.push_back(currNorm);
    }
    
    // Calculate the weights for the cutflow first!
    double evtWeight = 1.00000000;
    if (m_isWeighted) {
      // Normalization of inputs to 1 fb-1:
      evtWeight = (1000.0 * 
		   dmt->HGamEventInfoAuxDyn_crossSectionBRfilterEff["Nominal"] *
		   dmt->HGamEventInfoAuxDyn_weight["Nominal"] / 
		   nTotalEventsInFile);
      
      // Also add in branching ratio for DM samples:
      if (DMAnalysis::isDMSample(m_config, m_sampleName)) {
	evtWeight *= m_config->getNum("BranchingRatioHyy");
      }
    }
    
    // The mass parameter:
    double invariantMass = dmt->HGamEventInfoAuxDyn_m_yy["Nominal"] / 1000.0;
    int nLeptons = ((int)(dmt->HGamElectronsAuxDyn_pt["Nominal"]->size()) +
		    (int)(dmt->HGamMuonsAuxDyn_pt["Nominal"]->size()));
    
    // Then commence plotting for events passing inclusive H->yy selection:
    double varEtMiss = dmt->HGamEventInfoAuxDyn_TST_met["Nominal"] / 1000.0;
    double varPtYY = dmt->HGamEventInfoAuxDyn_pT_yy["Nominal"] / 1000.0;
    if (dmt->HGamEventInfoAuxDyn_cutFlow["Nominal"] >= 
	(int)m_config->getStrV("MxAODCutList").size()) {
      fillHist1D("pTyy", true, varPtYY, evtWeight, -1);
      fillHist1D("ETMiss", true, varEtMiss, evtWeight, -1);
      fillHist1D("ratioETMisspTyy", true, (varEtMiss / varPtYY), evtWeight, -1);
      fillHist1D("sumSqrtETMisspTyy", true,
		 sqrt(varEtMiss*varEtMiss+varPtYY*varPtYY), evtWeight, -1);
      fillHist1D("aTanRatio",true,TMath::ATan(varEtMiss/varPtYY),evtWeight,-1);
      fillHist1D("myy", true, invariantMass, evtWeight, -1);
      fillHist1D("njets", true, dmt->HGamEventInfoAuxDyn_Njets["Nominal"],
		 evtWeight, -1);
      fillHist1D("nleptons", true, nLeptons, evtWeight, -1);
    }
    
    // For systematic variations of the selection:
    if (getSystematics) {
      for (int i_s = 0; i_s < (int)systList.size(); i_s++) {
	sysSelectors[systList[i_s]]->passesCut("AllCuts", evtWeight);
      }
    }
    
    // Make sure events pass the selection:
    if (!selector->passesCut("AllCuts", evtWeight)) continue;
    
    // Save the categories:
    int currCate = selector->getCategoryNumber(m_config->getStr("cateScheme"),
					       evtWeight);
    
    // Fill the datasets:
    m_yy->setVal(invariantMass);
    if (m_isWeighted) {
      wt.setVal(evtWeight);
      m_cateData[currCate]->add(RooArgSet(*m_yy,wt), evtWeight);
    }
    else m_cateData[currCate]->add(*m_yy);
    
    // Write mass point to file:
    massFiles[currCate] << invariantMass << " " << evtWeight << std::endl;
    
    // Then commence plotting for PASSING events:
    fillHist1D("pTyy", false, varPtYY, evtWeight, currCate);
    fillHist1D("ETMiss", false, varEtMiss, evtWeight, currCate);
    fillHist1D("ratioETMisspTyy", false,(varEtMiss/varPtYY),evtWeight,currCate);
    fillHist1D("sumSqrtETMisspTyy", false, 
	       sqrt(varEtMiss*varEtMiss+varPtYY*varPtYY), evtWeight, currCate);
    fillHist1D("aTanRatio", false, TMath::ATan(varEtMiss/varPtYY), 
	       evtWeight, currCate);
    fillHist1D("myy", false, invariantMass, evtWeight, currCate);
    fillHist1D("njets", false, dmt->HGamEventInfoAuxDyn_Njets["Nominal"],
	       evtWeight, currCate);
    fillHist1D("nleptons", false, nLeptons, evtWeight, -1);
  }
  std::cout << "DMMassPoints: End of loop over input DMTree." << std::endl;
  
  // For systematic variations of the selection:
  if (getSystematics) {
    for (int i_s = 0; i_s < (int)systList.size(); i_s++) {
      sysSelectors[systList[i_s]]
	->saveCutflow(Form("%s/Systematics/cutflow_%s_%s.txt",
			   m_outputDir.Data(), systList[i_s].Data(),
			   m_sampleName.Data()), m_isWeighted);
      sysSelectors[systList[i_s]]
	->saveCategorization(Form("%s/Systematics/categorization_%s_%s_%s.txt",
				  m_outputDir.Data(), systList[i_s].Data(),
				  (m_config->getStr("cateScheme")).Data(),
				  m_sampleName.Data()), m_isWeighted);
    }
  }
  
  // Print the cutflow and category yields (weighted if MC):
  selector->printCutflow(m_isWeighted);
  selector->printCategorization(m_isWeighted);
  selector->saveCutflow(Form("%s/cutflow_%s.txt", m_outputDir.Data(),
			     m_sampleName.Data()), m_isWeighted);
  selector->saveCategorization(Form("%s/categorization_%s_%s.txt",
				    m_outputDir.Data(),
				    (m_config->getStr("cateScheme")).Data(),
				    m_sampleName.Data()), m_isWeighted);
  
  // Then retrieve the cutflow hist and save it:
  m_cutFlowHist = selector->retrieveCutflowHist(m_isWeighted);
  combineCutFlowHists();
  
  // Close output mass point files.
  for (int i_f = 0; i_f < m_config->getInt("nCategories"); i_f++) {
    massFiles[i_f].close();
  }
  
  // If options said to run locally, remove the files.
  if (m_options.Contains("CopyFile")) removeLocalFilesAndList(listName);

  // Finally, save the histograms (including variables and cutflow) to file:
  saveHists();
  
  std::cout << "DMMassPoints: Finished creating new mass points!" << std::endl;
  delete selector;
}

/**
   -----------------------------------------------------------------------------
   Load the mass points from text files that have already been produced. This is
   much faster than producing mass points from scratch, and is preferred. 
   @return - void.
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
    std::ifstream massFile(getMassPointsFileName(i_c));
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
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void DMMassPoints::printProgressBar(int index, int total) {
  if (index%10000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

/**
   -----------------------------------------------------------------------------
   Remove the local files and file list.
   @param listName - The name of the local file list.
*/
void DMMassPoints::removeLocalFilesAndList(TString listName) {
  TString currLine;
  std::ifstream localFile(listName);
  while (!localFile.eof()) {
    localFile >> currLine;
    if (!currLine.Contains("eos/atlas") && !currLine.EqualTo("")) {
      system(Form("rm %s", currLine.Data()));
    }
  }
  localFile.close();
  system(Form("rm %s", listName.Data()));
}
