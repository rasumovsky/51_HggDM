////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  AnaInfo.cxx                                                               //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
//  This class stores basic analysis data for the optimization.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AnaInfo.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the AnaInfo class.
   @param newIndex - The index of the analysis.
*/
AnaInfo::AnaInfo(int newIndex) {
  setIndex(newIndex);
  m_cutValues.clear();
  m_signals.clear();
  m_cutNames.clear();
  m_statValues.clear();
  setGood(true);
}

/**
   -----------------------------------------------------------------------------
   Add the cut to the list of cuts if it is not already contained.
   @param cutName - The name of the cut.
*/
void AnaInfo::addCut(TString cutName) {
  bool matched = false;
  std::vector<TString>::iterator iter;
  for (iter = m_cutNames.begin(); iter != m_cutNames.end(); iter++) {
    if (cutName.EqualTo(*iter)) {
      matched = true;
      break;
    }
  }
  if (!matched) m_cutNames.push_back(cutName);
}

/**
   -----------------------------------------------------------------------------
   Add the signal to the list of signals if it is not already contained.
   @param signal - The name of the signal.
*/
void AnaInfo::addSignal(TString signal) {
  bool matched = false;
  std::vector<TString>::iterator iter;
  for (iter = m_signals.begin(); iter != m_signals.end(); iter++) {
    if (signal.EqualTo(*iter)) {
      matched = true;
      break;
    }
  }
  if (!matched) m_cutNames.push_back(signal);
}

/**
   -----------------------------------------------------------------------------
   Returns the expected or observed p0 for the given signal.
   @param signal - The name of the signal.
   @param observed - True for observed statistic, false for expected.
*/
double AnaInfo::getAnaP0(TString signal, bool observed) {
  if (observed) return m_statValues[Form("%s_ObsP0", signal.Data())];
  else return m_statValues[Form("%s_ExpP0", signal.Data())];
}

/**
   -----------------------------------------------------------------------------
   Returns the expected or observed CL exclusion for the given signal.
   @param signal - The name of the signal.
   @param observed - True for observed statistic, false for expected.
*/
double AnaInfo::getAnaCL(TString signal, bool observed) {
  if (observed) return m_statValues[Form("%s_ObsCL", signal.Data())];
  else return m_statValues[Form("%s_ExpCL", signal.Data())];
}

/**
   -----------------------------------------------------------------------------
   Get a list of cuts that were defined for this analysis.
   @returns - A vector of cut names.
*/
std::vector<TString> AnaInfo::getCutList() {
  return m_cutNames;
}

/**
   -----------------------------------------------------------------------------
   @param cutName - The name of the cut to retrieve the value of. 
*/
double AnaInfo::getCutVal(TString cutName) {
  return m_cutValues[cutName];
}

/**
   -----------------------------------------------------------------------------
   Get the index of this analysis.
   @returns - The analysis index.
*/
int AnaInfo::getIndex() { 
  return m_index;
}

/**
   -----------------------------------------------------------------------------
   Get a list of signals defined for this analysis.
   @returns - A vector of signal names.
*/
std::vector<TString> AnaInfo::getSignalList() {
  return m_signals;
}

/**
   -----------------------------------------------------------------------------
   Get the statistic value for the given signal.
   @param signal - The name of the signal.
   @param statistic - The name of the statistic.
   @returns - the value of the statistic for the signal.
*/
double AnaInfo::getStatVal(TString signal, TString statistic) {
  return m_statValues[Form("%s_%s", signal.Data(), statistic.Data())];
}

/**
   -----------------------------------------------------------------------------
   Check whether analysis attributes are good.
   @returns - True iff the analysis completed successfully.
*/
bool AnaInfo::isGood() {
  return m_isGood;
}

/**
   -----------------------------------------------------------------------------
   Print the contents of this analysis.
   @param signal - The name of the signal to use for reference test statistics.
*/
void AnaInfo::printAna(TString signal) {
  
  std::cout << "\nAnaInfo: Printing analysis with values for signal " 
	    << signal << std::endl;
  std::cout << "\tIsGood? = " << m_isGood << std::endl;
  std::cout << "Analysis index = " << m_index << std::endl;/// ADD INTEX!
  std::map<TString,double>::iterator iter;
  for (iter = m_cutValues.begin(); iter != m_cutValues.end(); iter++) {
    TString currName = iter->first;
    double currVal = iter->second;
    std::cout << "\t" << currName << " = " << currVal << std::endl;
  }
  std::cout << "Exp CL = " << m_statValues[Form("%s_ExpCL",signal.Data())]
	    << std::endl;
  std::cout << "Exp p0 = " << m_statValues[Form("%s_ExpP0",signal.Data())]
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Set the value of the named analysis cut.
   @param cutName - The name of the cut.
   @param cutVal - The value of the cut.
*/
void AnaInfo::setCutVal(TString cutName, double cutVal) {
  m_cutValues[cutName] = cutVal;
  addCut(cutName);
}

/**
   -----------------------------------------------------------------------------
   Set whether the analysis is good or bad.
   @param isGood - True iff the analysis is good.
*/
void AnaInfo::setGood(bool isGood) {
  m_isGood = isGood;
}

/**
   -----------------------------------------------------------------------------
   Set the analysis index.
   @param index - The new analysis index.
*/
void AnaInfo::setIndex(int index) {
  m_index = index;
}

/**
   -----------------------------------------------------------------------------
   Set the value of the statistic for the given signal.
   @param signal - The name of the signal.
   @param statistic - The name of the statistic.
   @param value - The value for the specified statistic for the signal.
*/
void AnaInfo::setStatVal(TString signal, TString statistic, double value) {
  m_statValues[Form("%s_%s", signal.Data(), statistic.Data())] = value;
  addSignal(signal);
}
