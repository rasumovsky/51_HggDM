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
   @param newTree - the TTree which contains the sample.
*/
DMEvtSelect::DMEvtSelect(DMTree* newTree, TString newConfigFile) {
  std::cout << "DMEvtSelect: Initializing DMEvtSelect" << std::endl;
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  
  // ADD CUT HERE
  cutList.clear();
  cutList.push_back("photonPt");
  cutList.push_back("photonEta");
  cutList.push_back("photonIso");
  cutList.push_back("photonID");
  cutList.push_back("diphotonMass");
  if (m_config->getBool("LeptonVeto")) {
    cutList.push_back("leptonVeto");
  }
  cutList.push_back("diphotonPt");
  cutList.push_back("diphotonETMiss");
  cutList.push_back("looseCuts");//same as allCuts but no photonID or photonIso
  cutList.push_back("allCuts");
  
  // ADD CATE HERE ([name] = # categories):
  cateSchemesAndSizes.clear();
  cateSchemesAndSizes[m_config->getStr("cateScheme")]
    = m_config->getInt("nCategories");
  
  evtTree = newTree;
  selectionUsed = "";
  categoryUsed = "";
  
  // Reset event counters and initialize values to zero:
  clearCounters();
  std::cout << "DMEvtSelect: Successfully initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events in the specified category.
   @param cateScheme - the name of the categorization.
   @param cate - the category number. 
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
   @param cateScheme - the name of the categorization.
   @param cate - the category number.
   @returns - the number of weighted events in the category.
*/
double DMEvtSelect::getEventsPerCateWt(TString cateScheme, int cate) {
  if (cateExists(cateScheme)) {
    return cateCountWt[Form("%s_%d",cateScheme.Data(),cate)];
  }
  else return -1;
}

/** DEPRECATED
   Get the number of categories in the named categorization scheme.
   @param cateScheme - the name of the categorization.
   @returns - the number of categories in the categorization scheme.
int DMEvtSelect::getNCategories(TString cateScheme) {
  if (cateExists(cateScheme)) return cateSchemesAndSizes[cateScheme];
  else return -1;
}
*/

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events passing the specified cut.
   @param cutName - the name of the cut.
   @returns - the integer number of events passing the cut.
*/
int DMEvtSelect::getPassingEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountPass[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events passing the specified cut.
   @param cutName - the name of the cut.
   @returns - the weighted number of events passing the cut.
*/
double DMEvtSelect::getPassingEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountPassWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the (integer) number of events tested at the specified cut.
   @param cutName - the name of the cut.
   @returns - the integer number of events tested at the cut.
*/
int DMEvtSelect::getTotalEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountTot[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Get the weighted number of events tested at the specified cut.
   @param cutName - the name of the cut.
   @returns - the weighted number of events tested at the cut.
*/
double DMEvtSelect::getTotalEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountTotWt[cutName];
  else return -1;
}

/**
   -----------------------------------------------------------------------------
   Print the cutflow.
   @param weighted - true iff the event counts should be weighted.
*/
void DMEvtSelect::printCutflow(bool weighted) {
  std::cout << "Printing Cutflow: " << std::endl;
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)cutList.size(); i++) {
    if (selectionUsed.EqualTo("looseCuts") && cutList[i].EqualTo("allCuts")) {
      continue;
    }
    if (selectionUsed.EqualTo("allCuts") && cutList[i].EqualTo("looseCuts")) {
      continue;
    }
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
   Print the categories.
   @param weighted - true iff the event counts should be weighted.
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
   Save the cutflow.
   @param fileName - the output filename for the cutflow.
   @param weighted - true iff the event counts should be weighted.
*/
void DMEvtSelect::saveCutflow(TString fileName, bool weighted) {
  ofstream outFile(fileName);
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)cutList.size(); i++) {
    if (selectionUsed.EqualTo("looseCuts") && cutList[i].EqualTo("allCuts")) {
      continue;
    }
    if (selectionUsed.EqualTo("allCuts") && cutList[i].EqualTo("looseCuts")) {
      continue;
    }
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
   Save the categories.
   @param fileName - the output filename for the category yields.
   @param weighted - true iff the event counts should be weighted.
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
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this event belongs. 
   @param cateScheme - the name of the categorization.
   @returns - the category number for the event.
*/
int DMEvtSelect::getCategoryNumber(TString cateScheme) {
  return getCategoryNumber(cateScheme, 1.0);
}

/**
   -----------------------------------------------------------------------------
   Find the category in which this weighted event belongs. 
   @param cateScheme - the name of the categorization.
   @param weight - the event weight.
   @returns - the category number for the event.
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
    categoryUsed = "inclusive";
    currCate = 0;
  }
  // Split MET - low and high MET categories.
  else if (cateScheme.EqualTo("splitETMiss")) {
    categoryUsed = "splitETMiss";
    if (evtTree->HGamEventInfoAuxDyn_HighMet_MET_reb_TST > 140.0) currCate = 1;
    else currCate = 0;
  }
  // Ratio ETMiss/pT categorization:
  else if (cateScheme.EqualTo("RatioEtmPt")) {
    double ratioCut1 = m_config->getNum("RatioCut1");
    double ratioCut2 = m_config->getNum("RatioCut2");
    double currRatio = (evtTree->HGamEventInfoAuxDyn_HighMet_MET_reb_TST / 
			evtTree->HGamEventInfoAuxDyn_HighMet_yy_pt);
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
   Check whether an event passes the specified cut.
   @param cutName - the name of the cut.
   @returns - true iff the event passes the cut.
*/
bool DMEvtSelect::passesCut(TString cutName) {
  return passesCut(cutName, 1.0);
}

/**
   -----------------------------------------------------------------------------
   Check whether a weighted event passes the specified cut. Adds to the 
   selection counters automatically. WARNING! Calling "allCuts" or "looseCuts" 
   in conjunction with the other counters will lead to duplication.
   @param cutName - the name of the cut.
   @param weight - the event weight.
   @returns - true iff the event passes the cut.
*/
bool DMEvtSelect::passesCut(TString cutName, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutName)) return false;
  
  // The mass parameter:
  double invariantMass = evtTree->HGamEventInfoAuxDyn_HighMet_yy_m;
    
  // ADD CUT HERE:
  bool passes = true;
  // Cut on photon transverse momenta / diphoton mass:
  if (cutName.EqualTo("photonPt")) {
    passes = ((evtTree->HGamEventInfoAuxDyn_HighMet_y1_pt/invariantMass>0.35) &&
	      (evtTree->HGamEventInfoAuxDyn_HighMet_y2_pt/invariantMass>0.25));
  }
  // Cut on the photon pseudorapidities:
  else if (cutName.EqualTo("photonEta")) {
    passes = (evtTree->HGamEventInfoAuxDyn_HighMet_y1_eta < 2.37 && 
	      !(evtTree->HGamEventInfoAuxDyn_HighMet_y1_eta > 1.37 &&
		evtTree->HGamEventInfoAuxDyn_HighMet_y1_eta < 1.52) &&
	      evtTree->HGamEventInfoAuxDyn_HighMet_y2_eta < 2.37 && 
	      !(evtTree->HGamEventInfoAuxDyn_HighMet_y2_eta > 1.37 &&
		evtTree->HGamEventInfoAuxDyn_HighMet_y2_eta < 1.52));
  }
  // Cut on the calo/track isolation of the photons.
  else if (cutName.EqualTo("photonIso")) {
    //passes = (evtTree->HGamEventInfoAuxDyn_HighMet_y1_track_iso < 2.6 && 
    //	      evtTree->HGamEventInfoAuxDyn_HighMet_y2_track_iso < 2.6);
    passes = true;
  }
  // Cut on the ID variable of the photons.
  else if (cutName.EqualTo("photonID")) {
    //passes = (evtTree->HGamEventInfoAuxDyn_HighMet_y1_ID == 2 &&
    //	      evtTree->HGamEventInfoAuxDyn_HighMet_y2_ID == 2);
    passes = true;
  }
  // Cut on the diphoton invariant mass:
  else if (cutName.EqualTo("diphotonMass")) {
    passes = (invariantMass > m_config->getNum("DMMyyRangeLo") &&
	      invariantMass < m_config->getNum("DMMyyRangeHi"));
  }
  else if (cutName.EqualTo("leptonVeto") && m_config->getBool("LeptonVeto")) {
    passes = (evtTree->HGamEventInfoAuxDyn_HighMet_lep_n2 == 0);
  }
  // Cut on the diphoton transverse momentum:
  else if (cutName.EqualTo("diphotonPt")) {
    passes = (evtTree->HGamEventInfoAuxDyn_HighMet_yy_pt > 
	      m_config->getNum("AnaCutDiphotonPt"));
  }
  // Cut on the event missing transverse energy:
  else if (cutName.EqualTo("diphotonETMiss")) {
    passes = (evtTree->HGamEventInfoAuxDyn_HighMet_MET_reb_TST > 
	      m_config->getNum("AnaCutETMiss"));
  }
  // Check whether event passes all of the cuts above:
  else if (cutName.EqualTo("allCuts")) {
    selectionUsed = "allCuts";
    for (int i = 0; i < (int)cutList.size(); i++) {
      if (cutList[i].EqualTo("allCuts") || cutList[i].EqualTo("looseCuts")) {
	continue;
      }
      else if (!passesCut(cutList[i], weight)) {
	passes = false;
	break;
      }
    }
  }
  // Check whether event passes all of the cuts except ID and isolation.
  else if (cutName.EqualTo("looseCuts")) {
    selectionUsed = "looseCuts";
    for (int i = 0; i < (int)cutList.size(); i++) {
      if (cutList[i].EqualTo("allCuts") || cutList[i].EqualTo("looseCuts") ||
	  cutList[i].EqualTo("photonIso") || cutList[i].EqualTo("photonID")) {
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
  }
  return passes;
}

/**
   -----------------------------------------------------------------------------
   Give the selector class a new TTree to handle.
   @param newTree - the TTree on which the cuts shall be tested.
*/
void DMEvtSelect::setTree(DMTree *newTree) {
  evtTree = newTree;
  std::cout << "DMEvtSelect: A new TTree has been linked." << std::endl;
}

/** 
   -----------------------------------------------------------------------------
   Check whether the specified cut has been defined.
   @param cutName - the name of the cut whose existence shall be questioned.
   @returns - true iff the cut exists.
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
    Check whether the specified category has been defined.
    @param cateScheme - the name of the categorization.
    @returns - true iff the categorization has been defined. 
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
