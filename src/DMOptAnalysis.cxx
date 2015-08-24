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
  
  std::cout << "DMOptAnalysis: Finished!" << std::endl;
}

/**
   -----------------------------------------------------------------------------
*/
TString DMOptAnalysis::getPrintName(TString originName) {
  if (originName.EqualTo("ATanRatio1")) {
    return "Lower tan^{-1}(#slash{E}_{T} / p_{T}^{#gamma#gamma}) cut";
  }
  else if (originName.EqualTo("ATanRatio1")) {
    return "Upper tan^{-1}(#slash{E}_{T} / p_{T}^{#gamma#gamma}) cut";
  }
  else {
    return originName;
  }
}

/**
   -----------------------------------------------------------------------------
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
    currAna->cutNameAndVal.clear();
    currAna->cutNameAndVal["ATanRatio1"] = TMath::ATan(cut1Val);
    currAna->cutNameAndVal["ATanRatio2"] = TMath::ATan(cut2Val);
    
    //storeCutData("ATanRatio1");
    //storeCutData("ATanRatio2");
    
    // Then loop over the signals:
    std::vector<TString> signalList = m_config->getStrV("sigDMModes");
    for (int i_DM = 0; i_DM < (int)signalList.size(); i_DM++) {
      currAna->signals.clear();
      currAna->signals.push_back(signalList[i_DM]);
      currAna->valuesExpCLN2.clear();
      currAna->valuesExpCLN1.clear();
      currAna->valuesExpCL.clear();
      currAna->valuesExpCLP1.clear();
      currAna->valuesExpCLP2.clear();
      currAna->valuesObsCL.clear();
      currAna->valuesExpP0.clear();
      currAna->valuesObsP0.clear();
      
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
	  currAna->valuesExpP0[signalList[i_DM]] = currExpP0;
	  currAna->valuesObsP0[signalList[i_DM]] = currObsP0;
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
	  currAna->valuesExpCLN2[signalList[i_DM]] = currExpCLN2;
	  currAna->valuesExpCLN1[signalList[i_DM]] = currExpCLN1;
	  currAna->valuesExpCL[signalList[i_DM]]   = currExpCL;
	  currAna->valuesExpCLP1[signalList[i_DM]] = currExpCLP1;
	  currAna->valuesExpCLP2[signalList[i_DM]] = currExpCLP2;
	  currAna->valuesObsCL[signalList[i_DM]]   = currObsCL;
	}
      }
      currCLFile.close();
      
      analysisList.push_back(currAna);
    }// End loop over DM signals
  }// End loop over summary file listing analyses.
  inputSummaryFile.close();
}

/**
   -----------------------------------------------------------------------------
*/
void DMOptAnalysis::plotOptimizationPoints(TString cutName1, TString cutName2) {
  TCanvas *can = new TCanvas("can", "can", 800, 800, 
			     11, 0, TMath::Pi()/2, 11, 0, TMath::Pi()/2);
  can->cd();
  TH2F *hPoints = new TH2F("hPoints", "hPoints");
  for (int i_a = 0; i_a < (int)analysisList.size(); i_a++) {
    if (analysisList[i_a]->isGood) {
      hPoints->Fill(analysisList[i_a]->cutNameAndVal[cutName1], 
		    analysisList[i_a]->cutNameAndVal[cutName2], -1.0);
    }
    else {
      hPoints->Fill(analysisList[i_a]->cutNameAndVal[cutName1],
		    analysisList[i_a]->cutNameAndVal[cutName1], +1.0);
    }
  }
  hPoints->GetXaxis()->SetTitle(getPrintName(cutName1));
  hPoints->GetYaxis()->SetTitle(getPrintName(cutName2));
  hPoints->Draw("COLZ");
  can->Print("%s/plot_optimization_points.eps", Form(m_outputDir.Data()));
  can->Clear();
  delete can;
  delete hPoints;
}

/**
   -----------------------------------------------------------------------------

void DMOptAnalysis::storeCutData(TString cutName) {
  // First check to see if vector for cut already exists in map
  // if so, load that. 
  // if not, create new and clear.
  // Check to see if vector contains cutname
  // if not, add then sort
}
*/
