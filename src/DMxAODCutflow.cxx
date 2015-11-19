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
DMxAODCutflow::DMxAODCutflow(TString fileName, TString configFileName) {
  std::cout << "DMxAODCutflow: Initializing DMxAODCutflow" << std::endl;
  std::cout << "\tfrom file: " << fileName << std::endl;
  
  m_inputFile = TFile::Open(fileName);
  
  // Decide whether the sample is skimmed or unskimmed:
  Config *config = new Config(configFileName);
  m_unskimmed = !DMAnalysis::isSkimmed(config, fileName);
  
  // NOTE: at the moment, the weighted cutflow histogram is the 
  // weighted_nodalitz histogram. Not sure whether this is the correct choice.
  // need to check.
  
  // Find the cutflow histograms from the file based on limited name info:
  m_histCuts_weighted = NULL;
  m_histCuts_unweighted = NULL;
  //m_histCuts_weighted_NoDalitz = NULL;
  TIter next(m_inputFile->GetListOfKeys());
  TObject *currObj;
  while ((currObj = (TObject*)next())) {
    TString currName = currObj->GetName();
    if (currName.Contains("CutFlow") && currName.Contains("weighted")
    	&& currName.Contains("noDalitz")) {
      m_histCuts_weighted = (TH1F*)m_inputFile->Get(currName);
    }
    //else if (currName.Contains("CutFlow") && currName.Contains("weighted")
    //	     && currName.Contains("noDalitz")) {
    //m_histCuts_weighted_NoDalitz = (TH1F*)m_inputFile->Get(currName);
    //}
    else if (currName.Contains("CutFlow") && !currName.Contains("weighted")
	     && !currName.Contains("noDalitz")) {
      m_histCuts_unweighted = (TH1F*)m_inputFile->Get(currName);
    }
  }
  if (m_histCuts_weighted && m_histCuts_unweighted) {
    std::cout << "DMxAODCutflow: Sample is weighted." << std::endl;
    m_isWeighted = true;
  }
  else if (!m_histCuts_weighted && m_histCuts_unweighted) {
    std::cout << "DMxAODCutflow: Sample is unweighted." << std::endl;
    m_isWeighted = false;
    m_histCuts_weighted = NULL;
  }
  else if (!m_histCuts_unweighted) {
    std::cout << "DMxAODCutflow: Error loading file: " << fileName << std::endl;
    exit(0);
  }
  
  // If "NoDalitz" cutflow exists, use it instead:
  //if (m_histCuts_weighted_NoDalitz) {
  //  m_histCuts_weighted = m_histCuts_weighted_NoDalitz;
  //}
  
  // Reset event counters and initialize values to zero:
  cutList.clear();
  passCounter_weighted.clear();
  passCounter_unweighted.clear();
  
  // Then fill with histogram contents:
  nCuts = 0;
  for (int i_b = 1; i_b <= (int)m_histCuts_unweighted->GetNbinsX(); i_b++) {
    TString currName
      = (TString)m_histCuts_unweighted->GetXaxis()->GetBinLabel(i_b);
    if (currName.Length() > 1) {
      cutList.push_back(currName);
      if (m_isWeighted) {
	passCounter_weighted[currName]
	  = m_histCuts_weighted->GetBinContent(i_b);
      }
      passCounter_unweighted[currName]
	= m_histCuts_unweighted->GetBinContent(i_b);
      nCuts++;
    }
  }
  
  // Print the MxAOD cutflow.
  printxAODCutflow();
  
  std::cout << "DMxAODCutflow: Successfully initialized!" << std::endl;
  delete currObj;
  delete config;
}

/**
   -----------------------------------------------------------------------------
   Check whether the specified cut has been defined.
   @param cutName - the name of the cut whose existence shall be questioned.
   @return - true iff the cut exists.
*/
bool DMxAODCutflow::cutExists(TString cutName) {
  // Checks if there is a key corresponding to cutName in the map: 
  bool nonExistent = (passCounter_unweighted.find(cutName) == 
		      passCounter_unweighted.end());
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
   @return - the name of the cut.
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
   @return - the order of the cut in the xAOD cutflow.
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
   @param cutName - The name of the cut.
   @return - The number of passing events.
*/
double DMxAODCutflow::getEventsPassingCut(int order) {
  return getEventsPassingCut(getCutNameByOrder(order));
}

/**
   -----------------------------------------------------------------------------
   Get the number of events passing a particular stage of the xAOD cutflow.
   @param cutName - The name of the cut.
   @return - The number of passing events.
*/
double DMxAODCutflow::getEventsPassingCut(TString cutName) {
  if (cutExists(cutName)) {
    if (m_isWeighted) return passCounter_weighted[cutName];
    else return passCounter_unweighted[cutName];
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
  if (m_isWeighted) return m_histCuts_weighted;
  else return m_histCuts_unweighted;
}

/**
   -----------------------------------------------------------------------------
   Get the percentage of events passing the named cut.
   @param cutName - The name of the cut.
   @return - The percentage of events passing just this cut.
*/
double DMxAODCutflow::getPercentPassingCut(TString cutName) {
  // Divide current passing events by the initial events:
  int priorCutOrder = getCutOrderByName(cutName) - 1;
  if (priorCutOrder < 3) priorCutOrder = getCutOrderByName(cutName);
  
  TString priorCutName = getCutNameByOrder(priorCutOrder);
  double result = 100.0 * (getEventsPassingCut(cutName) / 
			   getEventsPassingCut(priorCutName));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the acceptance X efficiency of the analysis up to this point.
   @param cutName - the name of the cut.
   @return - the percentage acceptance times efficiency total.
*/
double DMxAODCutflow::getAccXEffAtCut(TString cutName) {
  // Divide current passing events by the initial events:
  double result = 0.0;
  if (getEventsPassingCut(getCutNameByOrder(3)) > 0.0) {
    result = 100.0 * (getEventsPassingCut(cutName) /
		      getEventsPassingCut(getCutNameByOrder(3)));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the total number of events in the file, for normalization of MC.
*/
double DMxAODCutflow::nTotalEventsInFile() {
  double nTotalEvents = 0.0;
  if (!m_isWeighted) nTotalEvents = m_histCuts_unweighted->GetBinContent(1);//?
  else {
    if (m_unskimmed) nTotalEvents = m_histCuts_weighted->GetBinContent(3);
    else {
      nTotalEvents = (m_histCuts_weighted->GetBinContent(3) * 
		      m_histCuts_unweighted->GetBinContent(2) /
		      m_histCuts_unweighted->GetBinContent(1));
    }
  }
  return nTotalEvents;
}

/**
   -----------------------------------------------------------------------------
   Print the stored cuts and the event values.
*/
void DMxAODCutflow::printxAODCutflow() {
  std::cout << "DMxAODCutflow: Printing the xAOD cutflow." << std::endl;
  for (int i_c = 0; i_c < nCuts; i_c++) {
    double currPassVal = m_isWeighted ? 
      passCounter_weighted[cutList[i_c]] : passCounter_unweighted[cutList[i_c]];
    std::cout << "\t" << cutList[i_c] << "\t" << currPassVal 
	      << "\t" << getAccXEffAtCut(cutList[i_c]) << "%." << std::endl;
  }
}
