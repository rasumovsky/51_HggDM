////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  AnaCollection.cxx                                                         //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
//  This program stores a collection of analyses for easy manipulation.       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AnaCollection.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the AnaCollection class.
   @param newIndex - The index of the analysis.
*/
AnaCollection::AnaCollection() {
  m_analyses.clear();
}

/**
   -----------------------------------------------------------------------------
   Get the optimal analysis according to some test statistic for a given signal.
   @param signal - The name of the signal.
   @param statistic - The name of the statistic.
   @param minimize - True iff the statistic should be minimized (e.g. p0).
   @returns - A pointer to the most optimal AnaInfo.
*/
AnaInfo* AnaCollection::getOptimalAnalysis(TString signal, TString statistic) {
  bool minimize = statistic.Contains("P0") ? true : false;
  AnaInfo *optimalAna = m_analyses[0];
  // Loop over analyses to see which has the best value of the test statistic.
  for (int i_a = 1; i_a < (int)m_analyses.size(); i_a++) {
    AnaInfo *currAna = m_analyses[i_a];
    double currStatVal = currAna->getStatVal(signal, statistic);
    double optStatVal = optimalAna->getStatVal(signal, statistic);
    // Only consider "good" analyses":
    if (!currAna->isGood()) continue;
    // Check whether the current analysis improves the sensitivity:
    if ((currStatVal <= optStatVal && minimize) || 
	(currStatVal >= optStatVal && !minimize)) {
      optimalAna = currAna;
    }
  }
  return optimalAna;
}

/**
   -----------------------------------------------------------------------------
   Get the number of analyses in this collection.
   @returns - The number of analyses.
*/
int AnaCollection::nAnalyses() {
  return (int)m_analyses.size();
}

/**
   -----------------------------------------------------------------------------
   Get the number of bad analyses in this collection (analyses which failed to
   load data properly, or didn't converge in the fitting.
   @returns - The number of bad analyses in the collection.
*/
int AnaCollection::nBadAnalyses() {
  int badCount = 0;
  for (int i_a = 0; i_a < (int)m_analyses.size(); i_a++) {
    AnaInfo *currAna = m_analyses[i_a];
    if (!currAna->isGood()) badCount++;
  }
  return badCount;
}

/**
   -----------------------------------------------------------------------------
   For iterative purposes.
*/
std::vector<AnaInfo*>::iterator AnaCollection::begin() {
  return m_analyses.begin();
}

/**
   -----------------------------------------------------------------------------
   For iterative purposes.
*/
std::vector<AnaInfo*>::iterator AnaCollection::end() {
  return m_analyses.end();
}

/**
   -----------------------------------------------------------------------------
   Add an analysis to the collection.
   @param newAnalysis - The analysis to add to this collection.
*/
void AnaCollection::addAnalysis(AnaInfo *newAnalysis) {
  m_analyses.push_back(newAnalysis);
}
