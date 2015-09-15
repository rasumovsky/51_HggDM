////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMxAODCutflow.cxx                                                   //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 23/04/2015                                                          //
//                                                                            //
//  This class reads the cutflow histograms created by the xAOD tool and can  //
//  return the number of events at each step.                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMxAODCutflow.h"

/**
   -----------------------------------------------------------------------------
   Initializes the tool from file.
   Initializes the tool and loads XS, BR values from files. 
   @param newTree - the TTree which contains the sample.
*/
DMxAODCutflow::DMxAODCutflow(TString fileName) {
  std::cout << "DMxAODCutflow: Initializing DMxAODCutflow" << std::endl;
  std::cout << "\tfrom file: " << fileName << std::endl;
  
  //TFile *m_inputFile = new TFile(fileName);
  //TH1F *m_histCuts = (TH1F*)m_inputFile->Get("HighMet_EventCutFlow");
  m_inputFile = new TFile(fileName);
  m_histCuts = (TH1F*)m_inputFile->Get("HighMet_EventCutFlow");
  
  if (!m_histCuts) {
    std::cout << "DMxAODCutflow: Error loading file: " << fileName << std::endl;
    exit(0);
  }
  
  // Reset event counters and initialize values to zero:
  cutList.clear();
  passCounter.clear();
  
  // Then fill with histogram contents:
  nCuts = 0;
  for (int i_b = 1; i_b <= (int)m_histCuts->GetNbinsX(); i_b++) {
    TString currName = (TString)m_histCuts->GetXaxis()->GetBinLabel(i_b);
    double currValue = m_histCuts->GetBinContent(i_b);
    if (currName.Length() > 1) {
      cutList.push_back(currName);
      passCounter[currName] = currValue;
      nCuts++;
    }
  }
  //delete m_histCuts;
  //delete m_inputFile;
  
  // Print the MxAOD cutflow.
  printxAODCutflow();
  
  std::cout << "DMxAODCutflow: Successfully initialized!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Check whether the specified cut has been defined.
   @param cutName - the name of the cut whose existence shall be questioned.
   @returns - true iff the cut exists.
*/
bool DMxAODCutflow::cutExists(TString cutName) {
  // Checks if there is a key corresponding to cutName in the map: 
  bool nonExistent = (passCounter.find(cutName) == passCounter.end());
  if (nonExistent) {
    std::cout << "DMxAODCutflow: Cut " << cutName << " not defined!"
	      << std::endl;
  }
  return !nonExistent;
}

/**
   -----------------------------------------------------------------------------
   Get the name of a cut based on its order in the xAOD cutflow, starting at 1.
   @param order - the order of the cut in the xAOD cutflow.
   @returns - the name of the cut.
*/
TString DMxAODCutflow::getCutNameByOrder(int order) {
  if (order > 0 && order <= nCuts) {
    return cutList[order-1];
  }
  else {
    std::cout << "DMxAODCutflow: Error! No cut of order " << order << std::endl;
    return "";
  }
}

/**
   -----------------------------------------------------------------------------
   Get the order of a cut based on its name.
   @param cutName - the name of the cut.
   @returns - the order of the cut in the xAOD cutflow.
*/
int DMxAODCutflow::getCutOrderByName(TString cutName) {
  if (cutExists(cutName)) {
    for (int i_c = 0; i_c < nCuts; i_c++) {
      if (cutName.EqualTo(cutList[i_c])) {
	return (i_c + 1);
      } 
    }
  }
  return 0;
}

/**
   -----------------------------------------------------------------------------
   Get the number of events passing a particular stage of the xAOD cutflow.
   @param cutName - the name of the cut.
   @returns the number of passing events.
*/
double DMxAODCutflow::getEventsPassingCut(int order) {
  return getEventsPassingCut(getCutNameByOrder(order));
}

/**
   -----------------------------------------------------------------------------
   Get the number of events passing a particular stage of the xAOD cutflow.
   @param cutName - the name of the cut.
   @returns the number of passing events.
*/
double DMxAODCutflow::getEventsPassingCut(TString cutName) {
  if (cutExists(cutName)) {
    return passCounter[cutName];
  }
  else {
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get a pointer to the cutflow histogram.
*/
TH1F* DMxAODCutflow::getHist() {
  return m_histCuts;
}

/**
   -----------------------------------------------------------------------------
   Get the percentage of events passing the named cut.
   @param cutName - the name of the cut.
   @returns - the percentage of events passing just this cut.
*/
double DMxAODCutflow::getPercentPassingCut(TString cutName) {
  // Divide current passing events by the initial events:
  int priorCutOrder = getCutOrderByName(cutName) - 1;
  if (priorCutOrder < 1) priorCutOrder = 1;
  
  TString priorCutName = getCutNameByOrder(priorCutOrder);
  double result = 100.0 * (getEventsPassingCut(cutName) / 
			   getEventsPassingCut(priorCutName));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the acceptance X efficiency of the analysis up to this point.
   @param cutName - the name of the cut.
   @returns - the percentage acceptance times efficiency total.
*/
double DMxAODCutflow::getAccXEffAtCut(TString cutName) {
  // Divide current passing events by the initial events:
  double result = 0.0;
  if (getEventsPassingCut(getCutNameByOrder(1)) > 0.0) {
    result = 100.0 * (getEventsPassingCut(cutName) /
		      getEventsPassingCut(getCutNameByOrder(1)));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Print the stored cuts and the event values.
*/
void DMxAODCutflow::printxAODCutflow() {
  std::cout << "DMxAODCutflow: Printing the xAOD cutflow." << std::endl;
  for (int i_c = 0; i_c < nCuts; i_c++) {
    std::cout << "\t" << cutList[i_c] << "\t" << passCounter[cutList[i_c]] 
	      << "\t" << getAccXEffAtCut(cutList[i_c]) << "%." << std::endl;
  }
}
