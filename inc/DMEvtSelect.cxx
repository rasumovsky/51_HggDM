////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMEvtSelect.cxx                                                     //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 16/03/2015                                                          //
//                                                                            //
//  This class is used to implement the cuts for the H->gg + DM analysis. In  //
//  addition to testing each cut, it has a built-in event counter.            //
//                                                                            //
//  To add a new cut, two modifications must be made at the locations labeled //
//  with the tag "ADD CUT HERE":                                              //
//    - add to the list of cutList in DMEvtSelect()                          //
//    - add to the implementation of cuts in passesCut()                      //
//                                                                            //
//  Similarly, you will need to update category definitions in the locations  //
//  identified with the tag "ADD CATE HERE":                                  //
//    - add to the list of cateSchemes in DMEvtSelect()                       //
//    - add to the implementation of categories in getCategory()              //
//                                                                            //
//  Note: the counter is a bit finnicky. Either check each step of the        //
//  individually OR use passesCut("all"), but don't use both. Otherwise, you  //
//  will be double-counting events, since the "all" cut recursively calls the //
//  passesCut() method.                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMEvtSelect.h"

/**
   -----------------------------------------------------------------------------
   Initializes the tool and loads XS, BR values from files. 
   @param newTree - The TTree which contains the sample.
*/
DMEvtSelect::DMEvtSelect(DMTree* newTree, TString newConfigFile) {
  std::cout << "DMEvtSelect: Initializing DMEvtSelect" << std::endl;
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  
  // Store the general and analysis specific cuts:
  m_cutList.clear();
  
  // Get the MxAOD cuts:
  std::vector<TString> cutFlowMxAODs = m_config->getStrV("MxAODCutList");
  for (int i_c = 0; i_c < (int)cutFlowMxAODs.size(); i_c++) {
    m_cutList.push_back(cutFlowMxAODs[i_c]);
  }
  
  // Then add analysis-specific cuts:
  if (m_config->getBool("LeptonVeto")) {
    m_cutList.push_back("LeptonVeto");
  }
  if (m_config->isDefined("AnaCutDiphotonPT")) {
    m_cutList.push_back("DiphotonPT");
  }
  if (m_config->isDefined("AnaCutETMiss")) {
    m_cutList.push_back("DiphotonETMiss");
  }
  m_cutList.push_back("AllCuts");
  
  // Load the category information:
  m_cateSchemesAndSizes.clear();
  m_cateSchemesAndSizes[m_config->getStr("cateScheme")]
    = m_config->getInt("nCategories");
  
  m_evtTree = newTree;
  
  // Set the systematic variation ("Nominal" by default):
  m_sysVariation = "Nominal";
  
  // Create a temporary cutflow histogram, to be replaced:
  m_cutFlowHist_weighted = new TH1F("cutflow_weighted", "cutflow_weighted",
				    m_cutList.size()-1, 0, m_cutList.size()-1);
  m_cutFlowHist_weighted->Sumw2(true);
  
  // Create a temporary cutflow histogram, to be replaced:
  m_cutFlowHist_unweighted = new TH1F("cutflow_unweighted","cutflow_unweighted",
				      m_cutList.size()-1,0,m_cutList.size()-1);
  
  // Reset event counters and initialize values to zero:
  clearCounters();
  std::cout << "DMEvtSelect: Successfully initialized!" << std::endl;
}

/** 
   -----------------------------------------------------------------------------
    Check whether the specified category has been defined.
    @param cateScheme - The name of the categorization.
    @return - True iff the categorization has been defined. 
*/
bool DMEvtSelect::cateExists(TString cateScheme) {
  // Checks if there is a key corresponding to cateScheme in the map: 
  bool nonExistent = (m_cateCount.find(Form("%s_0",cateScheme.Data()))
		      == m_cateCount.end());
  if (nonExistent) {
    std::cout << "DMEvtSelect: Category " << cateScheme << " not defined!"
	      << std::endl;
    std::cout << "DMEvtSelect: Printing m_cateCount and contents for reference."
	      << std::endl;
    for (std::map<TString,int>::iterator it = m_cateCount.begin(); 
	 it != m_cateCount.end(); it++) {
      std::cout << "\t" << it->first << "\t" << it->second << std::endl;
    }
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Clear the event counters.
*/
void DMEvtSelect::clearCounters() {
  // Clear the maps:
  m_evtCountPass.clear();
  m_evtCountPassWt.clear();
  m_evtCountTot.clear();
  m_evtCountTotWt.clear();
  m_cateCount.clear();
  m_cateCountWt.clear();
  // Then initialize all cut counters to zero:
  for (int i = 0; i < (int)m_cutList.size(); i++) {
    m_evtCountPass[m_cutList[i]] = 0;
    m_evtCountPassWt[m_cutList[i]] = 0.0;
    m_evtCountTot[m_cutList[i]] = 0;
    m_evtCountTotWt[m_cutList[i]] = 0.0;
  }
  // Then initialize all category counters to zero:
  for (std::map<TString,int>::iterator it = m_cateSchemesAndSizes.begin(); 
       it != m_cateSchemesAndSizes.end(); it++) {
    for (int j = 0; j < it->second; j++) {
      m_cateCount[Form("%s_%d", (it->first).Data(), j)] = 0;
      m_cateCountWt[Form("%s_%d", (it->first).Data(), j)] = 0.0;
    }
  }
  // Also empty the histograms:
  for (int i_b = 0; i_b <= m_cutFlowHist_weighted->GetNbinsX()+1; i_b++) {
    m_cutFlowHist_weighted->SetBinContent(i_b, 0.0);
    m_cutFlowHist_unweighted->SetBinContent(i_b, 0.0);
  }
  
}

/** 
   -----------------------------------------------------------------------------
   Check whether the specified cut has been defined.
   @param cutName - The name of the cut whose existence shall be questioned.
   @return - True iff the cut exists.
*/
bool DMEvtSelect::cutExists(TString cutName) {
  // Checks if there is a key corresponding to cutName in the map: 
  bool nonExistent = (m_evtCountTot.find(cutName) == m_evtCountTot.end());
  if (nonExistent) {
    std::cout << "DMEvtSelect: Cut " << cutName << " not defined!" << std::endl;
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Get the index of a particular cut.
   @param cutName - The name of the cut.
   @return - The index of the cut (the order in the cutflow).
*/
int DMEvtSelect::cutIndex(TString cutName) {
  for (int i_c = 0; i_c < (int)m_cutList.size(); i_c++) {
    if (cutName.EqualTo(m_cutList[i_c])) return i_c;
  }
  std::cout << "DMEvtSelect: Error! Improper cut: " << cutName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this event belongs. 
   @param cateScheme - The name of the categorization.
   @return - The category number for the event.
*/
int DMEvtSelect::getCategoryNumber(TString cateScheme) {
  return getCategoryNumber(cateScheme, 1.0);
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this weighted event belongs. 
   @param cateScheme - the name of the categorization.
   @param weight - The event weight.
   @return - The category number for the event.
*/
int DMEvtSelect::getCategoryNumber(TString cateScheme, double weight) {
  
  // check that the category is defined first. 
  if (!cateExists(cateScheme)) {
    std::cout << "DMEvtSelect: Categorization not defined: " << cateScheme 
	      << " " << weight << std::endl;
    exit(0);
  }
  
  // ADD CATE HERE:
  int currCate = -1;
  // Inclusive categorization - only 1 category.
  if (cateScheme.EqualTo("inclusive")) {
    currCate = 0;
  }
  // Split MET - low and high MET categories.
  else if (cateScheme.EqualTo("splitETMiss")) {
    if (m_evtTree->HGamEventInfoAuxDyn_TST_met[m_sysVariation] > 140000.0) {
      currCate = 1;
    }
    else currCate = 0;
  }
  // Ratio ETMiss/pT categorization:
  else if (cateScheme.EqualTo("RatioEtmPt")) {
    double ratioCut1 = m_config->getNum("RatioCut1");
    double ratioCut2 = m_config->getNum("RatioCut2");
    double currRatio = (m_evtTree->HGamEventInfoAuxDyn_TST_met[m_sysVariation] /
			m_evtTree->HGamEventInfoAuxDyn_pT_yy[m_sysVariation]);
    if (currRatio < ratioCut1) currCate = 0;
    else if (currRatio >= ratioCut1 && currRatio < ratioCut2) currCate = 1; 
    else if (currRatio >= ratioCut2) currCate = 2;
  }
  else if (cateScheme.EqualTo("combined")) {
    
    // High-MET region (ETMiss > 100 GeV):
    if (m_evtTree->HGamEventInfoAuxDyn_TST_met[m_sysVariation] >
	m_config->getNum("ETMissCut2")) {
      // Combine with high-PT to get mono-H category:
      if (m_evtTree->HGamEventInfoAuxDyn_pT_yy[m_sysVariation] > 
	  m_config->getNum("DiphotonPTCut2")) {
	currCate = 4;
      }
      else {
	currCate = 3;
      }
    }
    
    // Intermediate ETMiss and pTHard region (DEFINE pTHard!!):
    else if (m_evtTree->HGamEventInfoAuxDyn_TST_met[m_sysVariation] >
	     m_config->getNum("ETMissCut1") &&
	     m_evtTree->HGamEventInfoAuxDyn_pTHard[m_sysVariation] >
	     m_config->getNum("PTHardCut")) {
      currCate = 2;
    }
    
    // Rest category (everything else with pTyy > 15 GeV).
    else if (m_evtTree->HGamEventInfoAuxDyn_pT_yy[m_sysVariation] > 
	     m_config->getNum("DiphotonPTCut1")) {
      currCate = 1;
    }
    
    // Need to have a category for all other events.
    else {
      currCate = 0;
    }
  }
  if (currCate == -1) {
    std::cout << "DMEvtSelect: Categorization error for " << cateScheme
	      << std::endl;
    exit(0);
  }
  
  // Add to category counters:
  m_cateCount[Form("%s_%d",cateScheme.Data(),currCate)]++;
  m_cateCountWt[Form("%s_%d",cateScheme.Data(),currCate)] += weight;
  return currCate;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events in the specified category.
   @param cateScheme - The name of the categorization.
   @param cate - The category number. 
   @return - The weighted or unweighted number of events in specified category.
*/
int DMEvtSelect::getEventsPerCate(TString cateScheme, int cate) {
  if (cateExists(cateScheme)) {
    return m_cateCount[Form("%s_%d",cateScheme.Data(),cate)];
  }
  else {
    std::cout << "DMEvtSelect: ERROR! cannot retrieve events per category."
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events in the specified category.
   @param cateScheme - The name of the categorization.
   @param cate - The category number.
   @return - The number of weighted events in the category.
*/
double DMEvtSelect::getEventsPerCateWt(TString cateScheme, int cate) {
  if (cateExists(cateScheme)) {
    return m_cateCountWt[Form("%s_%d",cateScheme.Data(),cate)];
  }
  else {
    std::cout << "DMEvtSelect: ERROR! cannot retrieve events per category."
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events passing the specified cut.
   @param cutName - The name of the cut.
   @return - The integer number of events passing the cut.
*/
int DMEvtSelect::getPassingEvents(TString cutName) {
  if (cutExists(cutName)) return m_evtCountPass[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events passing the specified cut.
   @param cutName - The name of the cut.
   @return - The weighted number of events passing the cut.
*/
double DMEvtSelect::getPassingEventsWt(TString cutName) {
  if (cutExists(cutName)) return m_evtCountPassWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events tested at the specified cut.
   @param cutName - The name of the cut.
   @return - The integer number of events tested at the cut.
*/
int DMEvtSelect::getTotalEvents(TString cutName) {
  if (cutExists(cutName)) return m_evtCountTot[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events tested at the specified cut.
   @param cutName - the name of the cut.
   @return - The weighted number of events tested at the cut.
*/
double DMEvtSelect::getTotalEventsWt(TString cutName) {
  if (cutExists(cutName)) return m_evtCountTotWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   @return - The number of cuts that are defined for the analysis.
*/
int DMEvtSelect::nCuts() {
  return m_cutList.size();
}

/**
   -----------------------------------------------------------------------------
   Check whether an event passes the specified cut.
   @param cutName - The name of the cut.
   @return - True iff the event passes the cut.
*/
bool DMEvtSelect::passesCut(TString cutName) {
  return passesCut(cutName, 1.0);
}

/**
   -----------------------------------------------------------------------------
   Check whether a weighted event passes the specified cut. Adds to the 
   selection counters automatically. WARNING! Calling "AllCuts" in conjunction
   with the other counters will lead to duplication.
   @param cutName - The name of the cut.
   @param weight - The event weight.
   @return - True iff the event passes the cut.
*/
bool DMEvtSelect::passesCut(TString cutName, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutName)) return false;
  
  bool passes = true;
  
  // The mass parameter:
  double invariantMass = m_evtTree->HGamEventInfoAuxDyn_m_yy[m_sysVariation];
  
  bool isMxAODCut = false;
  int currCutIndex = 0;
  std::vector<TString> cutFlowMxAODs = m_config->getStrV("MxAODCutList");
  for (currCutIndex = 0; currCutIndex < (int)cutFlowMxAODs.size(); 
       currCutIndex++) {
    if (cutName.EqualTo(cutFlowMxAODs[currCutIndex])) {
      isMxAODCut = true;
      break;
    }
  }
  
  // MxAOD cuts:
  if (isMxAODCut) {
    passes = (m_evtTree->HGamEventInfoAuxDyn_cutFlow[m_sysVariation] >
	      currCutIndex);
  }
  // Lepton Veto Cut:
  else if (cutName.EqualTo("LeptonVeto") && m_config->getBool("LeptonVeto")) {
    passes = ((m_evtTree->HGamElectronsAuxDyn_pt)->size() == 0 &&
	      (m_evtTree->HGamMuonsAuxDyn_pt)->size() == 0);
  }
  // Cut on the diphoton transverse momentum:
  else if (cutName.EqualTo("DiphotonPT")) {
    passes = (m_evtTree->HGamEventInfoAuxDyn_pT_yy[m_sysVariation] >
	      m_config->getNum("AnaCutDiphotonPT"));
  }
  // Cut on the event missing transverse energy:
  else if (cutName.EqualTo("DiphotonETMiss")) {
    passes = (m_evtTree->HGamEventInfoAuxDyn_TST_met[m_sysVariation] > 
	      m_config->getNum("AnaCutETMiss"));
  }
  // Check whether event passes all of the cuts above:
  else if (cutName.EqualTo("AllCuts")) {
    for (int i = 0; i < (int)m_cutList.size(); i++) {
      if (m_cutList[i].EqualTo("AllCuts")) {
	continue;
      }
      else if (!passesCut(m_cutList[i], weight)) {
	passes = false;
	break;
      }
    }
  }
  
  // Add to total counters:
  m_evtCountTot[cutName]++;
  m_evtCountTotWt[cutName] += weight;
  
  // Add to passing counters:
  if (passes) {
    m_evtCountPass[cutName]++;
    m_evtCountPassWt[cutName] += weight;
    m_cutFlowHist_weighted->Fill(cutIndex(cutName), weight);
    m_cutFlowHist_unweighted->Fill(cutIndex(cutName), 1.0);
  }
  return passes;
}

/**
   -----------------------------------------------------------------------------
   Print the categories.
   @param weighted - True iff the event counts should be weighted.
*/
void DMEvtSelect::printCategorization(bool weighted) {
  std::cout << "Printing Categories: " << std::endl;
  // iterate over category names:
  for (std::map<TString,int>::iterator it = m_cateSchemesAndSizes.begin(); 
       it != m_cateSchemesAndSizes.end(); it++){
    std::cout << "\t" << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	std::cout << m_cateCountWt[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
      else {
	std::cout << m_cateCount[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
    }
    std::cout << std::endl;
  }
}

/**
   -----------------------------------------------------------------------------
   Print the cutflow.
   @param weighted - True iff the event counts should be weighted.
*/
void DMEvtSelect::printCutflow(bool weighted) {
  std::cout << "Printing Cutflow: " << std::endl;
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)m_cutList.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      std::cout << "\t" << m_cutList[i] << "\t" 
		<< m_evtCountPassWt[m_cutList[i]] << " / " 
		<< m_evtCountTotWt[m_cutList[i]] << std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      std::cout << "\t" << m_cutList[i] << "\t" << m_evtCountPass[m_cutList[i]]
		<< " / " << m_evtCountTot[m_cutList[i]] << std::endl;
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Export the histogram showing the cutflow.
   @param weighted - True iff you want the weighted cutflow histogram.
   @return - The cutflow histogram.
*/
TH1F* DMEvtSelect::retrieveCutflowHist(bool weighted) {
  if (weighted) return m_cutFlowHist_weighted;
  else return m_cutFlowHist_unweighted;
}

/**
   -----------------------------------------------------------------------------
   Save the categories.
   @param fileName - The output filename for the category yields.
   @param weighted - True iff the event counts should be weighted.
*/
void DMEvtSelect::saveCategorization(TString fileName, bool weighted) {
  ofstream outFile(fileName);
  // iterate over category names:
  for (std::map<TString,int>::iterator it = m_cateSchemesAndSizes.begin();
       it != m_cateSchemesAndSizes.end(); it++){
    outFile << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	outFile << m_cateCountWt[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
      else {
	outFile << m_cateCount[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
    }
    outFile << std::endl;
  }
  outFile.close();
}

/**
   -----------------------------------------------------------------------------
   Save the cutflow.
   @param fileName - The output filename for the cutflow.
   @param weighted - True iff the event counts should be weighted.
*/
void DMEvtSelect::saveCutflow(TString fileName, bool weighted) {
  ofstream outFile(fileName);
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)m_cutList.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      outFile << "\t" << m_cutList[i] << "\t" << m_evtCountPassWt[m_cutList[i]]
	      << " / " << m_evtCountTotWt[m_cutList[i]] << std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      outFile << "\t" << m_cutList[i] << "\t" << m_evtCountPass[m_cutList[i]]
	      << " / " << m_evtCountTot[m_cutList[i]] << std::endl;
    }
  }
  outFile.close();
}

/**
   -----------------------------------------------------------------------------
   Set which systematic uncertainty variation to use with the TTree. 
   @param sysVariation - The name of the systematic variation to use for the 
   selection ("Nominal" by default).
*/
void DMEvtSelect::setSysVariation(TString sysVariation) {
  m_sysVariation = sysVariation;
}

/**
   -----------------------------------------------------------------------------
   Give the selector class a new TTree to handle.
   @param newTree - The TTree on which the cuts shall be tested.
*/
void DMEvtSelect::setTree(DMTree *newTree) {
  m_evtTree = newTree;
  std::cout << "DMEvtSelect: A new TTree has been linked." << std::endl;
}
