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
  loadOptimizationData(Form("%s/%s/DMMaster/",
			    (m_config->getStr("masterOutput")).Data(), 
			    (m_config->getStr("jobName")).Data()));
  
  plotOptimizationPoints("ATanRatio1", "ATanRatio2");

  // Create optimization plots for every signal type:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_s = 0; i_s < (int)sigDMModes.size(); i_s++) {
    system(Form("mkdir -vp %s/%s",m_outputDir.Data(),(sigDMModes[i_s]).Data()));
    plotCutsAndStat(sigDMModes[i_s],"ATanRatio1","ATanRatio2","ExpCL",false);
    plotCutsAndStat(sigDMModes[i_s],"ATanRatio1","ATanRatio2","ExpP0",true);
  }
  
  std::cout << "DMOptAnalysis: Finished!" << std::endl;
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
*/
void DMOptAnalysis::getHistBinsAndRange(TString cutName, int &bins, double &min,
					double &max) {
  std::vector<double> uniqueCutValues; uniqueCutValues.clear();
  // Loop over all analysis points:
  for (int i_a = 1; i_a < (int)m_analysisList.size(); i_a++) {
    AnalysisAttributes *currAna = m_analysisList[i_a];
    double currCutValue = currAna->cutValues[cutName];
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
   For a given signal, based upon a given test statistic, identify the most
   sensitive analysis.
   @param signal - The signal sample for plotting.
   @param statistic - The name of the test statistic for the z-axis.
   @param minimize - True iff the statistic should be minimized (e.g. p0).
   @returns - The analysis index.
*/
int DMOptAnalysis::getOptAnaIndex(TString signal, TString statistic, 
				  bool minimize) {
  int optimalIndex = 0;
  AnalysisAttributes *optimalAnalysis = m_analysisList[optimalIndex];
  for (int i_a = 1; i_a < (int)m_analysisList.size(); i_a++) {
    AnalysisAttributes *currAna = m_analysisList[i_a];
    if ((currAna->statValues[mapKey(signal, statistic)] <
	 optimalAnalysis->statValues[mapKey(signal, statistic)] && minimize) || 
	(currAna->statValues[mapKey(signal, statistic)] >
	 optimalAnalysis->statValues[mapKey(signal, statistic)] && !minimize)) {
      optimalIndex = i_a;
      optimalAnalysis = m_analysisList[optimalIndex];
    }
  }
  return optimalIndex;
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
  
  int jobIndex; double cut1Val; double cut2Val;
  ifstream inputSummaryFile(Form("%s/jobSummary.txt", directory.Data()));
  while (!inputSummaryFile.eof()) {
    inputSummaryFile >> jobIndex >> cut1Val >> cut2Val;
    AnalysisAttributes *currAna;
    currAna->isGood = true;
    currAna->index = jobIndex;
    currAna->cutValues.clear();
    currAna->cutValues["ATanRatio1"] = TMath::ATan(cut1Val);
    currAna->cutValues["ATanRatio2"] = TMath::ATan(cut2Val);
    
    // Then loop over the signals:
    std::vector<TString> signalList = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)signalList.size(); i_DM++) {
      currAna->signals.clear();
      currAna->signals.push_back(signalList[i_DM]);
      currAna->statValues.clear();
      
      // Load p0:
      TString currP0FileName = Form("%s/single_files/p0_%d/p0_values_%s.txt",
				    directory.Data(), jobIndex,
				    signalList[i_DM].Data());
      ifstream currP0File(currP0FileName);
      if (!currP0File) {
	std::cout << "Problem loading file " << currP0FileName << std::endl;
	currAna->isGood = false;
      }
      else {
	TString currP0Name; double currExpP0; double currObsP0;
	while (!currP0File.eof()) {
	  currP0File >> currP0Name >> currExpP0 >> currObsP0;
	  currAna->statValues[mapKey(signalList[i_DM], "ExpP0")] = currExpP0;
	  currAna->statValues[mapKey(signalList[i_DM], "ObsP0")] = currObsP0;
	}
      }
      currP0File.close();
      
      // Load CL:
      TString currCLFileName = Form("%s/single_files/CL_%d/CL_values_%s.txt",
				    directory.Data(), jobIndex,
				    signalList[i_DM].Data());
      ifstream currCLFile(currCLFileName);
      if (!currCLFile) {
	std::cout << "Problem loading file " << currCLFileName << std::endl;
	currAna->isGood = false;
      }
      else {
	TString currCLName; double currObsCL; double currExpCLN2;
	double currExpCLN1; double currExpCL; double currExpCLP1;
	double currExpCLP2;
	while (!currP0File.eof()) {
	  currCLFile >> currCLName >> currObsCL >> currExpCLN2 >> currExpCLN1
		     >> currExpCL >> currExpCLP1 >> currExpCLP2;
	  currAna->statValues[mapKey(signalList[i_DM],"ExpCLN2")] = currExpCLN2;
	  currAna->statValues[mapKey(signalList[i_DM],"ExpCLN1")] = currExpCLN1;
	  currAna->statValues[mapKey(signalList[i_DM],"ExpCL")]   = currExpCL;
	  currAna->statValues[mapKey(signalList[i_DM],"ExpCLP1")] = currExpCLP1;
	  currAna->statValues[mapKey(signalList[i_DM],"ExpCLP2")] = currExpCLP2;
	  currAna->statValues[mapKey(signalList[i_DM],"ObsCL")]   = currObsCL;
	}
      }
      currCLFile.close();
      
      m_analysisList.push_back(currAna);
    }// End loop over DM signals
  }// End loop over summary file listing analyses.
  inputSummaryFile.close();
}

/**
   -----------------------------------------------------------------------------
   Get the map key string for a given signal and test statistic.
   @param signal - The signal sample for plotting.
   @param statistic - The name of the test statistic for the z-axis.
   @returns - The map key.
*/
TString DMOptAnalysis::mapKey(TString signal, TString statistic) {
  return Form("%s_%s", signal.Data(), statistic.Data());
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
				    TString cutNameY, TString statistic,
				    bool minimize) {
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  int xBins; int yBins; double xMin; double xMax; double yMin; double yMax;
  // Loop over the cuts to get the range:
  getHistBinsAndRange(cutNameX, xBins, xMin, xMax);
  getHistBinsAndRange(cutNameY, yBins, yMin, yMax);
  TH2F *hCont = new TH2F("hCont","hCont",xBins,xMin,xMax,yBins,yMin,yMax);
  
  double xOpt=-999; double yOpt=-999; double zOpt=-999; double zNonOpt=-999;
  for (int i_a = 0; i_a < (int)m_analysisList.size(); i_a++) {
    AnalysisAttributes *currAna = m_analysisList[i_a];
    if (currAna->isGood) {
      hCont->Fill(currAna->cutValues[cutNameX], 
		  currAna->cutValues[cutNameY],
		  currAna->statValues[statistic]);
      if ((minimize && currAna->statValues[mapKey(signal,statistic)] < zOpt) || 
	  (!minimize && currAna->statValues[mapKey(signal,statistic)] > zOpt)) {
	xOpt = currAna->cutValues[cutNameX];
	yOpt = currAna->cutValues[cutNameY];
	zOpt = currAna->statValues[statistic];
      }
      else if ((minimize && 
		currAna->statValues[mapKey(signal,statistic)] > zNonOpt) || 
	       (!minimize && 
		currAna->statValues[mapKey(signal,statistic)] < zNonOpt)) {
	zNonOpt = currAna->statValues[mapKey(signal,statistic)];
      }
    }
  }
  hCont->Draw("surf1");
  
  TPolyLine3D *pl3d1 = new TPolyLine3D(2);
  pl3d1->SetLineColor(kRed);
  pl3d1->SetLineWidth(3);
  pl3d1->SetPoint(0, xMin, yOpt, hCont->GetZaxis()->GetXmin());
  pl3d1->SetPoint(1, xMin, yOpt, zOpt);
  pl3d1->SetPoint(2, xOpt, yOpt, zOpt);
  pl3d1->SetPoint(3, xOpt, yMin, zOpt);
  pl3d1->SetPoint(4, xOpt, yMin, hCont->GetZaxis()->GetXmin());
  pl3d1->Draw("SAME");
  
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
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  TH2F *hPoints = new TH2F("hPoints", "hPoints",
			   11, 0, (TMath::Pi()/2.0), 11, 0, (TMath::Pi()/2));
  for (int i_a = 0; i_a < (int)m_analysisList.size(); i_a++) {
    if (m_analysisList[i_a]->isGood) {
      hPoints->Fill(m_analysisList[i_a]->cutValues[cutName1], 
		    m_analysisList[i_a]->cutValues[cutName2], -1.0);
    }
    else {
      hPoints->Fill(m_analysisList[i_a]->cutValues[cutName1],
		    m_analysisList[i_a]->cutValues[cutName1], +1.0);
    }
  }
  hPoints->GetXaxis()->SetTitle(DMAnalysis::getPrintVarName(cutName1));
  hPoints->GetYaxis()->SetTitle(DMAnalysis::getPrintVarName(cutName2));
  hPoints->Draw("COLZ");
  can->Print(Form("%s/plot_optimization_points.eps", m_outputDir.Data()));
  can->Clear();
  delete can;
  delete hPoints;
}
