////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: HDMEvtSelect.cxx                                                    //
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
//    - add to the list of cutnames in HDMEvtSelect()                         //
//    - add to the implementation of cuts in passesCut()                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HDMEvtSelect.h"

/**
   Initializes the tool and loads XS, BR values from files. 
*/
HDMEvtSelect::HDMEvtSelect(HDMTree* newTree) {
  evtTree = newTree;
  
  ///////////// ADD CUT HERE /////////////
  cutnames = {"photonPt",
	      "photonEta",
	      "diphotonMass",
	      "diphotonPt",
	      "diphotonETMiss"};
  ////////////////////////////////////////
  
  clearCounters();
  std::cout << "HDMEvtSelect: Successfully initialized!" << std::endl;
}

/**
   Get the (integer) number of events passing the specified cut.
*/
int HDMEvtSelect::getPassingEvents(TString cutname) {
  if (cutExists(cutname)) return evtCountPass[cutname];
  else return 0;
}

/**
   Get the weighted number of events passing the specified cut.
*/
double HDMEvtSelect::getPassingEventsWt(TString cutname) {
  if (cutExists(cutname)) return evtCountPassWt[cutname];
  else return 0;
}

/**
   Get the (integer) number of events tested at the specified cut.
*/
int HDMEvtSelect::getTotalEvents(TString cutname) {
  if (cutExists(cutname)) return evtCountTot[cutname];
  else return 0;
}

/**
   Get the weighted number of events tested at the specified cut.
*/
double HDMEvtSelect::getTotalEventsWt(TString cutname) {
  if (cutExists(cutname)) return evtCountTotWt[cutname];
  else return 0;
}

/**
   Clear the event counters.
*/
void HDMEvtSelect::clearCounters() {
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
bool HDMEvtSelect::passesCut(TString cutname) {
  return passesCut(cutname, 1.0);
}

/**
   Check whether a weighted event passes the specified cut.
*/
bool HDMEvtSelect::passesCut(TString cutname, double weight) {
  
  // check that map exists first.
  if (!cutExists(cutname)) return false;
  
  ///////////// ADD CUT HERE /////////////
  bool passes;
  if (cutname.Contains("photonPt")) {
    passes = (pt1/mgg > 0.35 && pt2/mgg > 0.25);
  }
  else if (cutname.Contains("photonEta")) {
    passes = (eta1 < 2.5 && eta2 < 2.5);
  }
  else if (cutname.Contains("diphotonMass")) {
    passes = (mgg > 105000 && mgg < 160000);
  }
  else if (cutname.Contains("diphotonPt")) {
    passes = (ptgg > 120000);
  }
  else if (cutname.Contains("diphotonETMiss")) {
    passes = (met > 120000);
  }
  ////////////////////////////////////////
  
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

/** 
    Check whether the specified cut has been defined.
*/
bool HDMEvtSelect::cutExists(TString cutname) {
  // Checks if there is a key corresponding to cutname in the map: 
  bool exists = (evtCountTot.find(cutname) == evtCountTot.end());
  if (!exists) {
    std::cout << "HDMEvtSelect: Cut not defined!" << std::endl;
  }
  return exists;
}
