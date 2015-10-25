////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DMOptAnalysis.cxx                                                         //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 21/08/2015                                                          //
//                                                                            //
//  This program loads many analyses (all with variations) and compares the   //
//  results for optimization purposes.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMOptAnalysis.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the DMOptAnalysis class.
   @param newConfigFile - The analysis configuration file.
*/
DMOptAnalysis::DMOptAnalysis(TString newConfigFile) {
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  TString jobName = m_config->getStr("jobName");
  
  // set input and output directories:
  m_outputDir = Form("%s/%s/DMOptAnalysis", 
		     (m_config->getStr("masterOutput")).Data(), jobName.Data());
  // Create output directory:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
    
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
  
  // Load the analysis optimization overview:
  m_anaCollection = new AnaCollection();
  loadOptimizationData(Form("%s/%s/DMMaster/",
			    (m_config->getStr("masterOutput")).Data(), 
			    (m_config->getStr("jobName")).Data()));
  
  // Plot all variable combinations:
  for (int i_c1 = 0; i_c1 < (int)m_optimizedCutList.size(); i_c1++) {
    for (int i_c2 = 0; i_c2 < (int)m_optimizedCutList.size(); i_c2++) {
      if (i_c1 == i_c2) continue;
      plotOptimizationPoints(m_optimizedCutList[i_c1],m_optimizedCutList[i_c2]);
    }
  }
  
  // Create optimization plots for every signal type:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_s = 0; i_s < (int)sigDMModes.size(); i_s++) {
    system(Form("mkdir -vp %s/%s",m_outputDir.Data(),(sigDMModes[i_s]).Data()));
    plotCutsAndStat(sigDMModes[i_s], m_optimizedCutList[0], 
      		    m_optimizedCutList[1], "ExpCL");
    plotCutsAndStat(sigDMModes[i_s], m_optimizedCutList[0],
		    m_optimizedCutList[1], "ExpP0");
  }
  
  std::cout << "DMOptAnalysis: Finished loading!" << std::endl;
  
  // Find the optimal analysis for each signal:
  for (int i_s = 0; i_s < (int)sigDMModes.size(); i_s++) {
    AnaInfo *optimalAna
      = m_anaCollection->getOptimalAnalysis(sigDMModes[i_s], "ExpP0");
    optimalAna->printAna(sigDMModes[i_s]);
  }
}

/**
   -----------------------------------------------------------------------------
   Checks whether the given double is contained in the vector, and adds it
   if it is not found.
   @param currVector - The current vector of doubles.
   @param newDouble - The double to check for membership in the list.
   @returns - A list of unique doubles.
*/
std::vector<double> DMOptAnalysis::checkDoubleList(std::vector<double> currList,
						   double newDouble) {
  for (int i_d = 0; i_d < (int)currList.size(); i_d++) {
    if (fabs(currList[i_d] - newDouble) <= 0.001) return currList;
  }
  currList.push_back(newDouble);
  return currList;
}

/**
   -----------------------------------------------------------------------------
   Get the axis nBins and range information for a variable.
   @param cutName - The name of the variable.
   @param bins - The number of bins.
   @param min - The minimum of the variable.
   @param max - The maximum of the variable.
   @returns - Updated values for bins, min, max passed by reference.
*/
void DMOptAnalysis::getHistBinsAndRange(TString cutName, int &bins, double &min,
					double &max) {
  std::vector<double> uniqueCutValues; uniqueCutValues.clear();
  // Loop over all analysis points:
  std::vector<AnaInfo*>::iterator iter;
  for (iter = m_anaCollection->begin(); iter != m_anaCollection->end(); iter++){
    AnaInfo *currAna = *iter;
    double currCutValue = currAna->getCutVal(cutName);
    uniqueCutValues = checkDoubleList(uniqueCutValues, currCutValue);
  }
  bins = uniqueCutValues.size();
  double tempMin = minEntry(uniqueCutValues);
  double tempMax = maxEntry(uniqueCutValues);
  double increment = (tempMax - tempMin) / ((double)bins-1.0);
  // So that the histogram bins are centered. 
  min = tempMin - (0.5 * increment);
  max = tempMax + (0.5 * increment);
}

/**
   -----------------------------------------------------------------------------
   List the contents of a directory.
   @param directory - The input directory.
   @returns - A vector containing the directory contents.
*/
std::vector<TString> DMOptAnalysis::listDirectoryContents(TString directory) {
  std::vector<TString> result; result.clear();
  system(Form("ls %s | tee templist.txt", directory.Data()));
  TString currFile;
  ifstream dirListFile("templist.txt");
  while (!dirListFile.eof()) {
    dirListFile >> currFile;
    result.push_back(Form("%s/%s", directory.Data(), currFile.Data()));
  }
  system("rm templist.txt");
  return result;
}

/**
   -----------------------------------------------------------------------------
   Load the optimization data from a given directory.
   @param directory- The input directory.
*/
void DMOptAnalysis::loadOptimizationData(TString directory) {
  std::cout << "DMOptAnalysis: Loading optimization data files." << std::endl;
    
  bool isFirstLine = true;
  //int jobIndex; double cut1Val; double cut2Val;
  std::string currStrLine;
  TString currLine;
  ifstream inputSummaryFile(Form("%s/jobSummary.txt", directory.Data()));
  if (!inputSummaryFile.is_open()) {
    std::cout << "DMOptAnalysis: Failed to open jobSummary.txt" << std::endl;
    exit(0);
  }
  while (getline(inputSummaryFile, currStrLine)) {
    currLine = TString(currStrLine);
    std::vector<TString> vectorizedLine = vectorizeTString(currLine, " ");
    
    // The first line of the file stores cut names:
    if (isFirstLine) {
      m_optimizedCutList.clear();// A list to store the cut names
      for (int i_c = 1; i_c < (int)vectorizedLine.size(); i_c++) {
	m_optimizedCutList.push_back(vectorizedLine[i_c]);
      }
      isFirstLine = false;
    }
    
    // The remaining lines store job indices and cut values:
    else {
      AnaInfo *currAna = new AnaInfo((vectorizedLine[0]).Atoi());
      for (int i_c = 1; i_c < (int)vectorizedLine.size(); i_c++) {
	double currentCutVal = (vectorizedLine[i_c]).Atof();
	if (currentCutVal > 1000.) currentCutVal = currentCutVal / 1000.0;
	currAna->setCutVal(m_optimizedCutList[i_c-1], currentCutVal);
      }
      
      // Then loop over the signals:
      std::vector<TString> signalList = m_config->getStrV("sigDMModes");
      for (int i_DM = 0; i_DM < (int)signalList.size(); i_DM++) {
	
	// Load p0:
	TString currP0FileName = Form("%s/single_files/p0_%d/p0_values_%s.txt",
				      directory.Data(), currAna->getIndex(),
				      signalList[i_DM].Data());
	ifstream currP0File(currP0FileName);
	if (!currP0File.is_open()) currAna->setGood(false);
	else {
	  TString currP0Name; double currExpP0; double currObsP0;
	  while (!currP0File.eof()) {
	    currP0File >> currP0Name >> currExpP0 >> currObsP0;
	    currAna->setStatVal(signalList[i_DM], "ExpP0", currExpP0);
	    currAna->setStatVal(signalList[i_DM], "ObsP0", currObsP0);
	  }
	}
	currP0File.close();
	
	// Load CL:
	TString currCLFileName = Form("%s/single_files/CL_%d/CL_values_%s.txt",
				      directory.Data(), currAna->getIndex(),
				      signalList[i_DM].Data());
	ifstream currCLFile(currCLFileName);
	if (!currCLFile.is_open()) currAna->setGood(false);
	else {
	  TString currCLName; double currObsCL; double currExpCLN2;
	  double currExpCLN1; double currExpCL; double currExpCLP1;
	  double currExpCLP2;
	  while (!currCLFile.eof()) {
	    currCLFile >> currCLName >> currObsCL >> currExpCLN2 >> currExpCLN1
		       >> currExpCL >> currExpCLP1 >> currExpCLP2;
	    currAna->setStatVal(signalList[i_DM], "ExpCLN2", currExpCLN2);
	    currAna->setStatVal(signalList[i_DM], "ExpCLN1", currExpCLN1);
	    currAna->setStatVal(signalList[i_DM], "ExpCL", currExpCL);
	    currAna->setStatVal(signalList[i_DM], "ExpCLP1", currExpCLP1);
	    currAna->setStatVal(signalList[i_DM], "ExpCLP2", currExpCLP2);
	    currAna->setStatVal(signalList[i_DM], "ObsCL", currObsCL);
	  }
	}
	currCLFile.close();
      }// End loop over DM signals
      // Add the current analysis to the collection:
      m_anaCollection->addAnalysis(currAna);
    }
  }// End loop over summary file listing analyses.
  inputSummaryFile.close();
  
  std::cout << "DMOptAnalysis: Finished loading analysis data." << std::endl;
  std::cout << "\tFailed to load " << m_anaCollection->nBadAnalyses() << " / " 
	    << m_anaCollection->nAnalyses() << " analyses." << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Returns the maximum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @returns - The maximum entry in the vector.
*/
double DMOptAnalysis::maxEntry(std::vector<double> currList) {
  double maximum = -999999.9;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter > maximum) maximum = *mIter;
  }
  return maximum;
}

/**
   -----------------------------------------------------------------------------
   Returns the minimum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @returns - The minimum entry in the vector.
*/
double DMOptAnalysis::minEntry(std::vector<double> currList) {
  double minimum = 999999.9;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter < minimum) minimum = *mIter;
  }
  return minimum;
}

/**
   -----------------------------------------------------------------------------

void DMOptAnalysis::plot2DScatter(TString quantity1, TString quantity2) {
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
 
  int xBins; int yBins; double xMin; double xMax; double yMin; double yMax;
  // Loop over the cuts to get the range:
  getHistBinsAndRange(cutNameX, xBins, xMin, xMax);
  getHistBinsAndRange(cutNameY, yBins, yMin, yMax);
  TH2F *hScatter = new TH2F("hScatter","hScatter",100,xMin,xMax,100,yMin,yMax);
  
  h_scat->SetMarkerSize(0.5);
  hScatter->Draw("scat=1.0");
  can->Print(Form("%s/plot_%s_vs_%s.eps", m_outputDir.Data(), quantity1.Data(),
		  quantity2.Data()));
  can->Clear();
  delete can;
  delete hScatter;
}
*/

/**
   -----------------------------------------------------------------------------
   Create a 2D surface plot with cuts on the x and y axes and a test statistic
   on the vertical axis.
   @param signal - The signal sample for plotting.
   @param cutNameX - The name of the x-axis cut.
   @param cutNameY - The name of the y-axis cut.
   @param statistic - The name of the test statistic for the z-axis.
   @param minimize - True iff the statistic should be minimized (e.g. p0).
*/
void DMOptAnalysis::plotCutsAndStat(TString signal, TString cutNameX,
				    TString cutNameY, TString statistic) {
  bool minimize = statistic.Contains("P0") ? true : false;
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  can->SetRightMargin(0.2);
  int xBins; int yBins; double xMin; double xMax; double yMin; double yMax;
  // Loop over the cuts to get the range:
  getHistBinsAndRange(cutNameX, xBins, xMin, xMax);
  getHistBinsAndRange(cutNameY, yBins, yMin, yMax);
  TH2F *hCont = new TH2F("hCont","hCont",xBins,xMin,xMax,yBins,yMin,yMax);
  hCont->GetXaxis()->SetTitle(DMAnalysis::getPrintVarName(cutNameX));
  hCont->GetYaxis()->SetTitle(DMAnalysis::getPrintVarName(cutNameY));
  hCont->GetZaxis()->SetTitle(DMAnalysis::getPrintVarName(statistic));
  hCont->GetZaxis()->SetTitleOffset(1.2);
  double xOpt=-999; double yOpt=-999; double zOpt=-999; double zNonOpt=-999;
  std::vector<AnaInfo*>::iterator iter;
  for (iter = m_anaCollection->begin(); iter != m_anaCollection->end(); iter++){
    AnaInfo *currAna = *iter;
    if (currAna->isGood()) {
      hCont->Fill(currAna->getCutVal(cutNameX), currAna->getCutVal(cutNameY),
		  currAna->getStatVal(signal, statistic));
   
    }
  }
  //hCont->Draw("surf1");
  hCont->Draw("COLZ");
  can->Print(Form("%s/%s/%s_%s_vs_%s.eps", m_outputDir.Data(), signal.Data(), 
		  statistic.Data(), cutNameX.Data(), cutNameY.Data()));
  can->Clear();
  delete can;
  delete hCont;
}

/**
   -----------------------------------------------------------------------------
   Show the points that are used for the optimization. If code ran successfully,
   use -1, otherwise failure is +1 for histogram filling.
   @param cutName1 - The name of the first cut.
   @param cutName2 - The name of the second cut.
*/
void DMOptAnalysis::plotOptimizationPoints(TString cutName1, TString cutName2) {
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  can->SetRightMargin(0.2);
  int bins1; int bins2; double min1; double min2; double max1; double max2;
  getHistBinsAndRange(cutName1, bins1, min1, max1);
  getHistBinsAndRange(cutName2, bins2, min2, max2);
  
  TH2F *hPoints = new TH2F("hPoints", "hPoints",
			   bins1, min1, max1, bins2, min2, max2);
  std::vector<AnaInfo*>::iterator iter;
  for (iter = m_anaCollection->begin(); iter != m_anaCollection->end(); iter++){
    AnaInfo *currAna = *iter;
    if (currAna->isGood()) {
      hPoints->Fill(currAna->getCutVal(cutName1), currAna->getCutVal(cutName2),
		    1.0);
    }
    else {
      //hPoints->Fill(currAna->getCutVal(cutName1),currAna->getCutVal(cutName1),
      //		    2.0);
    }
  }
  hPoints->GetXaxis()->SetTitle(DMAnalysis::getPrintVarName(cutName1));
  hPoints->GetYaxis()->SetTitle(DMAnalysis::getPrintVarName(cutName2));
  hPoints->GetZaxis()->SetTitle("fit status (0=bad)");
  hPoints->Draw("COLZ");
  can->Print(Form("%s/plot_points_%s_%s.eps", m_outputDir.Data(), 
		  cutName1.Data(), cutName2.Data()));
  can->Clear();
  delete can;
  delete hPoints;
}

/**
   -----------------------------------------------------------------------------
   Turn a list of strings into a vector.
*/
std::vector<TString> DMOptAnalysis::vectorizeTString(TString originString,
						     TString delim) {
  std::vector<TString> result; result.clear();
  TObjArray *tokenizedLine = originString.Tokenize(" ");
  for (int i_e = 0; i_e < tokenizedLine->GetEntries(); i_e++) {
    TObjString *obj = (TObjString*)tokenizedLine->At(i_e);
    TString subString = obj->GetString();
    subString.ReplaceAll("\t", "");
    subString.ReplaceAll(" ", "");
    result.push_back(subString);
  }
  return result;
}
