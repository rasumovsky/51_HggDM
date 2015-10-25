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
  cutList.clear();
  
  // Get the MxAOD cuts:
  std::vector<TString> cutFlowMxAODs = m_config->getStrV("MxAODCutList");
  for (int i_c = 0; i_c < (int)cutFlowMxAODs.size(); i_c++) {
    cutList.push_back(cutFlowMxAODs[i_c]);
  }
  
  // Then add analysis-specific cuts:
  if (m_config->getBool("LeptonVeto")) {
    cutList.push_back("LeptonVeto");
  }
  cutList.push_back("DiphotonPT");
  cutList.push_back("DiphotonETMiss");
  cutList.push_back("AllCuts");
  
  // Load the category information:
  cateSchemesAndSizes.clear();
  cateSchemesAndSizes[m_config->getStr("cateScheme")]
    = m_config->getInt("nCategories");
  
  evtTree = newTree;
  
  // Create a temporary cutflow histogram, to be replaced:
  m_cutFlowHist_weighted = new TH1F("cutflow_weighted", "cutflow_weighted",
				    cutList.size()-1, 0, cutList.size()-1);
  
  // Create a temporary cutflow histogram, to be replaced:
  m_cutFlowHist_unweighted = new TH1F("cutflow_unweighted","cutflow_unweighted",
				      cutList.size()-1, 0, cutList.size()-1);
  
  // Reset event counters and initialize values to zero:
  clearCounters();
  std::cout << "DMEvtSelect: Successfully initialized!" << std::endl;
}

/** 
   -----------------------------------------------------------------------------
    Check whether the specified category has been defined.
    @param cateScheme - The name of the categorization.
    @returns - True iff the categorization has been defined. 
*/
bool DMEvtSelect::cateExists(TString cateScheme) {
  // Checks if there is a key corresponding to cateScheme in the map: 
  bool nonExistent = (cateCount.find(Form("%s_0",cateScheme.Data()))
		      == cateCount.end());
  if (nonExistent) {
    std::cout << "DMEvtSelect: Category " << cateScheme << " not defined!"
	      << std::endl;
    
    std::cout << "Printing cateCount and contents for reference." << std::endl;
    std::map<TString,int>::iterator it;
    for (it = cateCount.begin(); it != cateCount.end(); it++) {
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
  evtCountPass.clear();
  evtCountPassWt.clear();
  evtCountTot.clear();
  evtCountTotWt.clear();
  cateCount.clear();
  cateCountWt.clear();
  // Then initialize all cut counters to zero:
  for (int i = 0; i < (int)cutList.size(); i++) {
    evtCountPass[cutList[i]] = 0;
    evtCountPassWt[cutList[i]] = 0.0;
    evtCountTot[cutList[i]] = 0;
    evtCountTotWt[cutList[i]] = 0.0;
  }
  // Then initialize all category counters to zero:
  std::map<TString,int>::iterator it;
  for (it = cateSchemesAndSizes.begin(); 
       it != cateSchemesAndSizes.end(); it++) {
    for (int j = 0; j < it->second; j++) {
      cateCount[Form("%s_%d", (it->first).Data(), j)] = 0;
      cateCountWt[Form("%s_%d", (it->first).Data(), j)] = 0.0;
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
   @returns - True iff the cut exists.
*/
bool DMEvtSelect::cutExists(TString cutName) {
  // Checks if there is a key corresponding to cutName in the map: 
  bool nonExistent = (evtCountTot.find(cutName) == evtCountTot.end());
  if (nonExistent) {
    std::cout << "DMEvtSelect: Cut " << cutName << " not defined!" << std::endl;
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Get the index of a particular cut.
   @param cutName - The name of the cut.
   @returns - The index of the cut (the order in the cutflow).
*/
int DMEvtSelect::cutIndex(TString cutName) {
  for (int i_c = 0; i_c < (int)cutList.size(); i_c++) {
    if (cutName.EqualTo(cutList[i_c])) return i_c;
  }
  std::cout << "DMEvtSelect: Error! Improper cut: " << cutName << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this event belongs. 
   @param cateScheme - The name of the categorization.
   @returns - The category number for the event.
*/
int DMEvtSelect::getCategoryNumber(TString cateScheme) {
  return getCategoryNumber(cateScheme, 1.0);
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this weighted event belongs. 
   @param cateScheme - the name of the categorization.
   @param weight - The event weight.
   @returns - The category number for the event.
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
    if (evtTree->HGamEventInfoAuxDyn_TST_met > 140000.0) currCate = 1;
    else currCate = 0;
  }
  // Ratio ETMiss/pT categorization:
  else if (cateScheme.EqualTo("RatioEtmPt")) {
    double ratioCut1 = m_config->getNum("RatioCut1");
    double ratioCut2 = m_config->getNum("RatioCut2");
    double currRatio = (evtTree->HGamEventInfoAuxDyn_TST_met / 
			evtTree->HGamEventInfoAuxDyn_pT_yy);
    if (currRatio < ratioCut1) currCate = 0;
    else if (currRatio >= ratioCut1 && currRatio < ratioCut2) currCate = 1; 
    else if (currRatio >= ratioCut2) currCate = 2;
  }
  
  if (currCate == -1) {
    std::cout << "DMEvtSelect: Categorization error for " << cateScheme
	      << std::endl;
    exit(0);
  }
  
    // Add to category counters:
  cateCount[Form("%s_%d",cateScheme.Data(),currCate)]++;
  cateCountWt[Form("%s_%d",cateScheme.Data(),currCate)] += weight;
  return currCate;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events in the specified category.
   @param cateScheme - The name of the categorization.
   @param cate - The category number. 
*/
int DMEvtSelect::getEventsPerCate(TString cateScheme, int cate) {
  if (cateExists(cateScheme)) {
    return cateCount[Form("%s_%d",cateScheme.Data(),cate)];
  }
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events in the specified category.
   @param cateScheme - The name of the categorization.
   @param cate - The category number.
   @returns - The number of weighted events in the category.
*/
double DMEvtSelect::getEventsPerCateWt(TString cateScheme, int cate) {
  if (cateExists(cateScheme)) {
    return cateCountWt[Form("%s_%d",cateScheme.Data(),cate)];
  }
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events passing the specified cut.
   @param cutName - The name of the cut.
   @returns - The integer number of events passing the cut.
*/
int DMEvtSelect::getPassingEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountPass[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events passing the specified cut.
   @param cutName - The name of the cut.
   @returns - The weighted number of events passing the cut.
*/
double DMEvtSelect::getPassingEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountPassWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events tested at the specified cut.
   @param cutName - The name of the cut.
   @returns - The integer number of events tested at the cut.
*/
int DMEvtSelect::getTotalEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountTot[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events tested at the specified cut.
   @param cutName - the name of the cut.
   @returns - The weighted number of events tested at the cut.
*/
double DMEvtSelect::getTotalEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountTotWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   @returns - The number of cuts that are defined for the analysis.
*/
int DMEvtSelect::nCuts() {
  return cutList.size();
}

/**
   -----------------------------------------------------------------------------
   Check whether an event passes the specified cut.
   @param cutName - The name of the cut.
   @returns - True iff the event passes the cut.
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
   @returns - True iff the event passes the cut.
*/
bool DMEvtSelect::passesCut(TString cutName, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutName)) return false;
  
  bool passes = true;
  
  // The mass parameter:
  double invariantMass = evtTree->HGamEventInfoAuxDyn_m_yy;
  
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
    passes = (evtTree->HGamEventInfoAuxDyn_cutFlow > currCutIndex);
  }
  // Lepton Veto Cut:
  else if (cutName.EqualTo("LeptonVeto") && m_config->getBool("LeptonVeto")) {
    passes = ((evtTree->HGamElectronsAuxDyn_pt)->size() == 0 && (evtTree->HGamMuonsAuxDyn_pt)->size() == 0);
  }
  // Cut on the diphoton transverse momentum:
  else if (cutName.EqualTo("DiphotonPT")) {
    passes = (evtTree->HGamEventInfoAuxDyn_pT_yy >
	      m_config->getNum("AnaCutDiphotonPT"));
  }
  // Cut on the event missing transverse energy:
  else if (cutName.EqualTo("DiphotonETMiss")) {
    passes = (evtTree->HGamEventInfoAuxDyn_TST_met > 
	      m_config->getNum("AnaCutETMiss"));
  }
  // Check whether event passes all of the cuts above:
  else if (cutName.EqualTo("AllCuts")) {
    for (int i = 0; i < (int)cutList.size(); i++) {
      if (cutList[i].EqualTo("AllCuts")) {
	continue;
      }
      else if (!passesCut(cutList[i], weight)) {
	passes = false;
	break;
      }
    }
  }
  
  // Add to total counters:
  evtCountTot[cutName]++;
  evtCountTotWt[cutName] += weight;
  
  // Add to passing counters:
  if (passes) {
    evtCountPass[cutName]++;
    evtCountPassWt[cutName] += weight;
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
  std::map<TString,int>::iterator it;
  for (it = cateSchemesAndSizes.begin(); it != cateSchemesAndSizes.end(); it++){
    std::cout << "\t" << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	std::cout << cateCountWt[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
      else {
	std::cout << cateCount[Form("%s_%d",(it->first).Data(), j)] << " ";
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
  for (int i = 0; i < (int)cutList.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      std::cout << "\t" << cutList[i] << "\t" << evtCountPassWt[cutList[i]]
		<< " / " << evtCountTotWt[cutList[i]] << std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      std::cout << "\t" << cutList[i] << "\t" << evtCountPass[cutList[i]]
		<< " / " << evtCountTot[cutList[i]] << std::endl;
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Export the histogram showing the cutflow.
   @param weighted - True iff you want the weighted cutflow histogram.
   @returns - The cutflow histogram.
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
  std::map<TString,int>::iterator it;
  for (it = cateSchemesAndSizes.begin(); it != cateSchemesAndSizes.end(); it++){
    outFile << "\t" << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	outFile << cateCount[Form("%s_%d",(it->first).Data(), j)] << " ";
      }
      else {
	outFile << cateCountWt[Form("%s_%d",(it->first).Data(), j)] << " ";
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
  for (int i = 0; i < (int)cutList.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      outFile << "\t" << cutList[i] << "\t" << evtCountPassWt[cutList[i]]
	      << " / " << evtCountTotWt[cutList[i]] << std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      outFile << "\t" << cutList[i] << "\t" << evtCountPass[cutList[i]]
	      << " / " << evtCountTot[cutList[i]] << std::endl;
    }
  }
  outFile.close();
}

/**
   -----------------------------------------------------------------------------
   Give the selector class a new TTree to handle.
   @param newTree - The TTree on which the cuts shall be tested.
*/
void DMEvtSelect::setTree(DMTree *newTree) {
  evtTree = newTree;
  std::cout << "DMEvtSelect: A new TTree has been linked." << std::endl;
}
