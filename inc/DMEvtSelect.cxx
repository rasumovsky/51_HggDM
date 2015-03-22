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
//    - add to the list of cateNames in DMEvtSelect()                         //
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
   Initializes the tool and loads XS, BR values from files. The input Tree is 
   set to NULL, and would have to be set later if used. 
*/
DMEvtSelect::DMEvtSelect() {
  DMEvtSelect(NULL);
}

/**
   Initializes the tool and loads XS, BR values from files. 
*/
DMEvtSelect::DMEvtSelect(DMTree* newTree) {
    
  // ADD CUT HERE
  cutList.clear();
  cutList.push_back("photonPt");
  cutList.push_back("photonEta");
  cutList.push_back("diphotonMass");
  cutList.push_back("diphotonPt");
  cutList.push_back("diphotonETMiss");
  cutList.push_back("allCuts");
    
  // ADD CATE HERE ([name] = # categories):
  cateNamesAndSizes.clear();
  cateNamesAndSizes["inclusive"] = 1;
  cateNamesAndSizes["splitETMiss"] = 2;
  
  evtTree = newTree;

  // Reset event counters and initialize values to zero:
  clearCounters();
  std::cout << "DMEvtSelect: Successfully initialized!" << std::endl;
}

/**
   Get the (integer) number of events in the specified category.
*/
int DMEvtSelect::getEventsPerCate(TString cateName, int cate) {
  if (cateExists(cateName)) {
    return cateCount[Form("%s_%d",cateName.Data(),cate)];
  }
  else return -1;
}

/**
   Get the weighted number of events in the specified category.
*/
double DMEvtSelect::getEventsPerCateWt(TString cateName, int cate) {
  if (cateExists(cateName)) {
    return cateCountWt[Form("%s_%d",cateName.Data(),cate)];
  }
  else return -1;
}

/**
   Get the number of categories in the named categorization scheme.
*/
int DMEvtSelect::getNCategories(TString cateName) {
  if (cateExists(cateName)) return cateNamesAndSizes[cateName];
  else return -1;
}

/**
   Get the (integer) number of events passing the specified cut.
*/
int DMEvtSelect::getPassingEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountPass[cutName];
  else return -1;
}

/**
   Get the weighted number of events passing the specified cut.
*/
double DMEvtSelect::getPassingEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountPassWt[cutName];
  else return -1;
}

/**
   Get the (integer) number of events tested at the specified cut.
*/
int DMEvtSelect::getTotalEvents(TString cutName) {
  if (cutExists(cutName)) return evtCountTot[cutName];
  else return -1;
}

/**
   Get the weighted number of events tested at the specified cut.
*/
double DMEvtSelect::getTotalEventsWt(TString cutName) {
  if (cutExists(cutName)) return evtCountTotWt[cutName];
  else return -1;
}

/**
   Print the cutflow.
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
   Print the categories.
*/
void DMEvtSelect::printCategorization(bool weighted) {
  std::cout << "Printing Categories: " << std::endl;
  // iterate over category names:
  std::map<TString,int>::iterator it;
  for (it = cateNamesAndSizes.begin(); it != cateNamesAndSizes.end(); it++) {
    std::cout << "\t" << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	std::cout << cateCount[Form("%s_%d",(it->first).Data(),it->second)]
		  << " ";
      }
      else {
	std::cout << cateCountWt[Form("%s_%d",(it->first).Data(),it->second)]
		  << " ";
      }
    }
    std::cout << std::endl;
  }
}

/**
   Save the cutflow.
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
   Save the categories.
*/
void DMEvtSelect::saveCategorization(TString fileName, bool weighted) {
  ofstream outFile(fileName);
  // iterate over category names:
  std::map<TString,int>::iterator it;
  for (it = cateNamesAndSizes.begin(); it != cateNamesAndSizes.end(); it++) {
    outFile << "\t" << it->first << " ";
    for (int j = 0; j < it->second; j++) {
      if (weighted) {
	outFile << cateCount[Form("%s_%d",(it->first).Data(),it->second)]
		<< " ";
      }
      else {
	outFile << cateCountWt[Form("%s_%d",(it->first).Data(),it->second)]
		<< " ";
      }
    }
    outFile << std::endl;
  }
  outFile.close();
}

/**
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
  for (it = cateNamesAndSizes.begin(); it != cateNamesAndSizes.end(); it++) {
    for (int j = 0; j < it->second; j++) {
      cateCount[Form("%s_%d",(it->first).Data(),it->second)] = 0;
      cateCountWt[Form("%s_%d",(it->first).Data(),it->second)] = 0.0;
    }
  }
}

/**
   Find the category in which this event belongs. 
*/
int DMEvtSelect::getCategoryNumber(TString cateName) {
  return getCategoryNumber(cateName, 1.0);
}

/**
   Find the category in which this weighted event belongs. 
*/
int DMEvtSelect::getCategoryNumber(TString cateName, double weight) {
  
  // check that the category is defined first. 
  if (!cateExists(cateName)) return -1;
  
  // ADD CATE HERE:
  int currCate = -1;
  // Inclusive categorization - only 1 category.
  if (cateName.Contains("inclusive")) {
    return 0;
  }
  // Split MET - low and high MET categories.
  else if (cateName.Contains("splitETMiss")) {
    if (evtTree->EventInfoAuxDyn_metref_final > 180.0) currCate = 0;
    else currCate = 1;
  }
  
  // Add to category counters:
  cateCount[Form("%s_%d",cateName.Data(),currCate)]++;
  cateCountWt[Form("%s_%d",cateName.Data(),currCate)] += weight;
  
  // Print error message before returning bad category value.
  if (currCate == -1) { 
    std::cout << "DMEvtSelect: Error! category not defined!" << std::endl;
  }
  return currCate;
}

/**
   Check whether an event passes the specified cut.
*/
bool DMEvtSelect::passesCut(TString cutName) {
  return passesCut(cutName, 1.0);
}

/**
   Check whether a weighted event passes the specified cut.
*/
bool DMEvtSelect::passesCut(TString cutName, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutName)) return false;
  
  // ADD CUT HERE:
  bool passes = true;
  // Cut on photon transverse momenta / diphoton mass:
  if (cutName.Contains("photonPt")) {
    passes = ((evtTree->EventInfoAuxDyn_y1_pt /
	       evtTree->EventInfoAuxDyn_m_yy > 0.35) &&
	      (evtTree->EventInfoAuxDyn_y2_pt /
	       evtTree->EventInfoAuxDyn_m_yy > 0.25));
  }
  // Cut on the photon pseudorapidities:
  else if (cutName.Contains("photonEta")) {
    passes = (evtTree->EventInfoAuxDyn_y1_eta < 2.5 && 
	      evtTree->EventInfoAuxDyn_y2_eta < 2.5);
  }
  // Cut on the diphoton invariant mass:
  else if (cutName.Contains("diphotonMass")) {
    passes = (evtTree->EventInfoAuxDyn_m_yy > 105.0 && 
	      evtTree->EventInfoAuxDyn_m_yy < 160.0);
  }
  // Cut on the diphoton transverse momentum:
  else if (cutName.Contains("diphotonPt")) {
    passes = (evtTree->EventInfoAuxDyn_pt_yy > 120.0);
  }
  // Cut on the event missing transverse energy:
  else if (cutName.Contains("diphotonETMiss")) {
    passes = (evtTree->EventInfoAuxDyn_metref_final > 120.0);
  }
  // Check whether event passes all of the cuts above:
  else if (cutName.Contains("all")) {
   
    for (int i = 0; i < (int)cutList.size(); i++) {
      if (cutList[i].Contains("all")) continue;
      if (!passesCut(cutList[i])) {
	passes = false;
	break;
      }
    }
  }
  
  // Add to total counters:
  evtCountTot[cutName]++;
  evtCountTotWt[cutName]+=weight;
  
  // Add to passing counters:
  if (passes) {
    evtCountPass[cutName]++;
    evtCountPassWt[cutName]+=weight;
  }
  return passes;
}

/**
   Give the selector class a new TTree to deal with.
*/
void DMEvtSelect::setTree(DMTree *newTree) {
  evtTree = newTree;
  std::cout << "DMEvtSelect: A new TTree has been linked." << std::endl;
}

/** 
   Check whether the specified cut has been defined.
*/
bool DMEvtSelect::cutExists(TString cutName) {
  // Checks if there is a key corresponding to cutName in the map: 
  bool exists = (evtCountTot.find(cutName) == evtCountTot.end());
  if (!exists) {
    std::cout << "DMEvtSelect: Cut not defined!" << std::endl;
  }
  return exists;
}

/** 
    Check whether the specified category has been defined.
*/
bool DMEvtSelect::cateExists(TString cateName) {
  // Checks if there is a key corresponding to cateName in the map: 
  bool exists = (cateCount.find(Form("%s_0",cateName.Data()))
		 == cateCount.end());
  if (!exists) {
    std::cout << "DMEvtSelect: Category not defined!" << std::endl;
  }
  return exists;
}
