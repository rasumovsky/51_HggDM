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
//    - add to the list of cutnames in DMEvtSelect()                         //
//    - add to the implementation of cuts in passesCut()                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMEvtSelect.h"

/**
   Initializes the tool and loads XS, BR values from files. 
*/
DMEvtSelect::DMEvtSelect(DMTree* newTree) {
    
  // ADD CUT HERE
  cutnames = {"photonPt",
	      "photonEta",
	      "diphotonMass",
	      "diphotonPt",
	      "diphotonETMiss",
	      "allCuts"};
  
  evtTree = newTree;
  recursiveCall = false;
  // Reset event counters and initialize values to zero:
  clearCounters();
  std::cout << "DMEvtSelect: Successfully initialized!" << std::endl;
}

/**
   Get the (integer) number of events passing the specified cut.
*/
int DMEvtSelect::getPassingEvents(TString cutname) {
  if (cutExists(cutname)) return evtCountPass[cutname];
  else return 0;
}

/**
   Get the weighted number of events passing the specified cut.
*/
double DMEvtSelect::getPassingEventsWt(TString cutname) {
  if (cutExists(cutname)) return evtCountPassWt[cutname];
  else return 0;
}

/**
   Get the (integer) number of events tested at the specified cut.
*/
int DMEvtSelect::getTotalEvents(TString cutname) {
  if (cutExists(cutname)) return evtCountTot[cutname];
  else return 0;
}

/**
   Get the weighted number of events tested at the specified cut.
*/
double DMEvtSelect::getTotalEventsWt(TString cutname) {
  if (cutExists(cutname)) return evtCountTotWt[cutname];
  else return 0;
}

/**
   Print the cutflow.
*/
void DMEvtSelect::printCutflow(bool weighted) {
  std::cout << "Printing Cutflow: " << std::endl;
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)cutNames.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      std::cout << "\t" << cutNames[i] << "\t" << evtCountPassWt[cutNames[i]]
		<< " / " << evtCountTotalWt[cutNames[i]] std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      std::cout << "\t" << cutNames[i] << "\t" << evtCountPass[cutNames[i]]
		<< " / " << evtCountTotal[cutNames[i]] std::endl;
    }
  }
}

/**
   Save the cutflow.
*/
void DMEvtSelect::saveCutflow(TString filename, bool weighted) {
  ofstream outFile(filename);
  // Loop over the cuts and print the name as well as the pass ratio:
  for (int i = 0; i < (int)cutNames.size(); i++) {
    // Print the weighted cutflow (for MC):
    if (weighted) {
      outFile << "\t" << cutNames[i] << "\t" << evtCountPassWt[cutNames[i]]
	      << " / " << evtCountTotalWt[cutNames[i]] std::endl;
    }
    // Print the unweighted cutflow (for data):
    else {
      outFile << "\t" << cutNames[i] << "\t" << evtCountPass[cutNames[i]]
	      << " / " << evtCountTotal[cutNames[i]] std::endl;
    }
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
  // Then initialize all counters to zero:
  for (int i = 0; i < (int)cutNames.size(); i++) {
    evtCountPass[cutNames[i]] = 0;
    evtCountPassWt[cutNames[i]] = 0.0;
    evtCountTot[cutNames[i]] = 0;
    evtCountTotalWt[cutNames[i]] = 0.0;
  }
}

/**
   check whether an event passes the specified cut.
*/
bool DMEvtSelect::passesCut(TString cutname) {
  return passesCut(cutname, 1.0);
}

/**
   Check whether a weighted event passes the specified cut.
*/
bool DMEvtSelect::passesCut(TString cutname, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutname)) return false;
  
  // ADD CUT HERE:
  bool passes = true;
  // Cut on photon transverse momenta / diphoton mass:
  if (cutname.Contains("photonPt")) {
    passes = (EventInfoAuxDyn.y1_pt/EventInfoAuxDyn.m_yy > 0.35 &&
	      EventInfoAuxDyn.y2_pt/EventInfoAuxDyn.m_yy > 0.25);
  }
  // Cut on the photon pseudorapidities:
  else if (cutname.Contains("photonEta")) {
    passes = (EventInfoAuxDyn.y1_eta < 2.5 && 
	      EventInfoAuxDyn.y2_eta < 2.5);
  }
  // Cut on the diphoton invariant mass:
  else if (cutname.Contains("diphotonMass")) {
    passes = (EventInfoAuxDyn.m_yy > 105.0 && 
	      EventInfoAuxDyn.m_yy < 160.0);
  }
  // Cut on the diphoton transverse momentum:
  else if (cutname.Contains("diphotonPt")) {
    passes = (EventInfoAuxDyn.pt_yy > 120.0);
  }
  // Cut on the event missing transverse energy:
  else if (cutname.Contains("diphotonETMiss")) {
    passes = (EventInfoAuxDyn.metref_final > 120.0);
  }
  // Check whether event passes all of the cuts above:
  else if (cutname.Contains("all")) {
    recursiveCall = true;
    for (int i = 0; i < cutNames.size(); i++) {
      if (cutNames[i].Contains("all")) continue;
      if (!passesCut(cutNames[i])) {
	passes = false;
      }
    }
    recursiveCall = false;
  }
  
  // The recursiveCall flag makes sure the event counters are not fucked up by
  // recursive calls to this method that checks the cuts. 
  if (!recursiveCall) {
    // Add to total counters:
    evtCountTot[cutname]++;
    evtCountTotWt[cutname]+=weight;
    
    // Add to passing counters:
    if (passes) {
      evtCountPass[cutname]++;
      evtCountPassWt[cutname]+=weight;
    }
    return passes;
  }
}

/** 
    Check whether the specified cut has been defined.
*/
bool DMEvtSelect::cutExists(TString cutname) {
  // Checks if there is a key corresponding to cutname in the map: 
  bool exists = (evtCountTot.find(cutname) == evtCountTot.end());
  if (!exists) {
    std::cout << "DMEvtSelect: Cut not defined!" << std::endl;
  }
  return exists;
}
