////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotVariables.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard,                                                     //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/04/2015                                                          //
//                                                                            //
//  This class builds the workspace for the dark matter analysis fits.        //
//                                                                            //
//  Options: "StackPlot" to stack all backgrounds together.                   //
//           "LogScale" to set logarithmic y-axis.                            //
//           "CombineSM" to combine H->yy SM production modes.                //
//           "Normalize" to scale to unity                                    //
//           "Scale2Data" to scale the plot to data normalization             //
//           "Smooth" to smooth the background histograms.                    //
//                                                                            //
//  Variables: "pTyy", "ETMiss", "ratioETMisspTyy", "aTanRatio", "myy",       //
//             "sumSqrtETMisspTyy", "dPhiyyETMiss", "njets", "cutFlowFull"    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// ROOT libraries:
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TString.h"
#include "Rtypes.h"

// Package libraries:
#include "Config.h"
#include "CommonFunc.h"
#include "CommonHead.h"
#include "DMAnalysis.h"

// A map to store all of the histograms:
std::vector<TH1F*> m_hists;
std::vector<TString> m_names;
std::vector<double> m_integrals;
std::vector<TString> m_matched;

Config *m_config;
TString m_outputDir;
TString m_options;

double nonDMSum_ALL;
double nonDMSum_PASS;
double nData_ALL;
double nData_PASS;

/**
   -----------------------------------------------------------------------------
   Safely load a histogram from a file.
   @param fName - The name of the file.
   @param hName - The name of the histogram.
*/
TH1* getHistFromFile(TString fName, TString hName) {
  TFile *file = TFile::Open(fName.Data(), "READ");
  if (!file) {
    std::cout << "PlotVariables:getHistogramFromFile() : Couldn't open file "
	      << fName.Data() << std::endl;
    return NULL;
  }
  TH1 *temp = (TH1*)file->Get(hName.Data());
  
  if (!temp) {
    std::cout << "DMxAODCutflow: Couldn't find histogram "
	      << hName.Data() << " in file " << fName.Data() << std::endl;
    return NULL;
  }
  
  bool status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  hName = "cloned_" + hName;
  TH1 *hist = (TH1*)temp->Clone(hName.Data());
  SafeDelete(file);
  TH1::AddDirectory(status);
  
  return hist;
}

/**
   -----------------------------------------------------------------------------
   Check if the current histogram name has already been matched and sorted.
   @param currName - The histogram name to check for matching.
   @returns - True iff currName is included in the list of matched histograms.
*/
bool isMatched(TString currName) {
  for (std::vector<TString>::iterator iter = m_matched.begin(); iter != m_matched.end(); iter++) {
    if (iter->EqualTo(currName)) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Load the histograms for a single sample and variable from file.
   @param sampleName - The name of the sample.
   @param varName - The name of the variable in the file.
*/
void loadSampleHistograms(TString sampleName, TString varName) {
  //std::cout << "PlotVariables: loadSampleHistograms(" << sampleName << ", " 
  //	    << varName << ")..." << std::endl;
  // Get the input directory and file names:
  TString inputDir = Form("%s/%s/DMMassPoints", 
			  (m_config->getStr("masterOutput")).Data(), 
			  (m_config->getStr("jobName")).Data());
  TString fName = Form("%s/hists_%s.root", inputDir.Data(), sampleName.Data());
  
  TH1F *hAll = NULL; 
  TH1F *hPass = NULL;
  if (varName.Contains("cutFlowFull")) {
    hAll = (TH1F*)getHistFromFile(fName,Form("%s_%s", varName.Data(),
					     sampleName.Data()));
    if (m_options.Contains("CombineSM") && 
	DMAnalysis::isSMSample(m_config, sampleName)) {
      m_names.push_back(Form("SMHiggs_%s_ALL", varName.Data()));
    }
    else {
      m_names.push_back(Form("%s_%s_ALL", sampleName.Data(), varName.Data()));
    }
    m_hists.push_back(hAll);
    m_integrals.push_back(hAll->Integral());
  }
  else if (m_options.Contains("CombineSM") && 
      DMAnalysis::isSMSample(m_config, sampleName)) {
    // Load histograms before and after cuts and in each category:
    hAll = (TH1F*)getHistFromFile(fName, Form("%s_ALL",varName.Data()));
    m_names.push_back(Form("SMHiggs_%s_ALL", varName.Data()));
    m_hists.push_back(hAll);
    m_integrals.push_back(hAll->Integral());
    hPass = (TH1F*)getHistFromFile(fName, Form("%s_PASS",varName.Data()));
    m_names.push_back(Form("SMHiggs_%s_PASS", varName.Data()));
    m_hists.push_back(hPass);
    m_integrals.push_back(hPass->Integral());
    for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
      TH1F *hCate
	= (TH1F*)getHistFromFile(fName,Form("%s_c%d_PASS",varName.Data(),i_c));
      m_names.push_back(Form("SMHiggs_%s_c%d_PASS", varName.Data(), i_c));
      m_hists.push_back(hCate);
      m_integrals.push_back(hCate->Integral());
    }
  }
  else {
    // Load histograms before and after cuts and in each category:
    hAll = (TH1F*)getHistFromFile(fName, Form("%s_ALL",varName.Data()));
    m_names.push_back(Form("%s_%s_ALL", sampleName.Data(), varName.Data()));
    m_hists.push_back(hAll);
    m_integrals.push_back(hAll->Integral());
    hPass = (TH1F*)getHistFromFile(fName, Form("%s_PASS",varName.Data()));
    m_names.push_back(Form("%s_%s_PASS", sampleName.Data(), varName.Data()));
    m_hists.push_back(hPass);
    m_integrals.push_back(hPass->Integral());
    for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
      TH1F *hCate
	= (TH1F*)getHistFromFile(fName,Form("%s_c%d_PASS",varName.Data(),i_c));
      m_names.push_back(Form("%s_%s_c%d_PASS", 
			     sampleName.Data(), varName.Data(), i_c));
      m_hists.push_back(hCate);
      m_integrals.push_back(hCate->Integral());
    }
  }
  
  // Strategic exit if histograms were NULL and are about to cause segfault:
  if ((!hAll || !hPass) && !(varName.Contains("cutFlowFull") && hAll)) {
    std::cout << "ERROR! Missing histogram." << std::endl;
    exit(0);
  }
  
  // Sum up everything but the DM signal for axis scaling:
  if (varName.Contains("cutFlowFull")) {
    nonDMSum_ALL += hAll->GetBinContent(3);
  }
  else {
    if (!DMAnalysis::isDMSample(m_config, sampleName) && 
	!sampleName.Contains("Data")) {
      nonDMSum_ALL += hAll->Integral();
      nonDMSum_PASS += hPass->Integral();
    }
  }
  
  // Also count data events:
  if (sampleName.Contains("Data")) {
    nData_ALL += hAll->Integral();
    if (!varName.Contains("cutFlowFull")) nData_PASS += hPass->Integral();
  }
}

/**
   -----------------------------------------------------------------------------
   Makes a combined plot with all samples drawn for a given variable.
   @param allEvents - True iff all events are in the hist (not just selected).
   @param varName - The name of the variable in the file.
*/
void makeCombinedPlot(bool allEvents, TString varName, int cateIndex) {
  std::cout << "PlotVariables: makeCombinedPlot(" << allEvents << ", "
  	    << varName << ", " << cateIndex << ")" << std::endl;
  
  TString endTag = allEvents ? "_ALL" : "_PASS";
  
  // Begin plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  // Create a legend to which we will add items:
  TLegend leg(0.5,0.73,0.92,0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.03);
  leg.SetNColumns(2);
  // Stacked histogram for all backgrounds (including SM Higgs):
  THStack hs("hs", "stacked histogram");
  
  // Color palettes for SM Higgs, DM Higgs, and Bkg:
  Color_t colorListDM[6] = {kRed, kRed+3, kRed-2, kRed-7, kRed+4, kRed+1};
  Color_t colorListSM[6] = {kBlue,kBlue-1,kBlue+3,kBlue-3,kBlue+5,kBlue-5};
  Color_t colorListBkg[2] = {kGreen-2, kGreen+2};
  
  Color_t colorStack[14] = {kViolet+6, kTeal-9, kOrange+1, kAzure-2, kYellow-9, 
			    kRed-7, kTeal-3, kAzure-8, kMagenta-10, kOrange+7,  
			    kBlue-10, kGreen-6, kOrange-2,kRed+2};
  
  // Loop over the saved histograms:
  int index = 0; int indexSM = 0; int indexDM = 0; int indexBkg = 0;
  for (int i_s = 0; i_s < (int)m_hists.size(); i_s++) {
    TString currName = m_names[i_s];
    
    // Histogram selection:
    if (!currName.Contains(endTag)) {
      continue;
    }
    else if (cateIndex >= 0 && !currName.Contains(Form("c%d",cateIndex))) {
      continue;
    }
    else if (cateIndex < 0) {
      bool badSample = false;
      for (int i_c = 0; i_c < m_config->getStr("nCategories"); i_c++) {
	if (currName.Contains(Form("c%d",i_c))) {
	  badSample = true;
	  break;
	}
      }
      if (badSample) continue;
    }
    
    //std::cout << "\tSample index " << index << ", " << currName << std::endl;
    
    // Get pointer to the histogram and modify the name:
    TH1F *currHist = m_hists[i_s];
    currName = currName.ReplaceAll("_"+varName, "");
    currName = currName.ReplaceAll(endTag, "");
    currName = currName.ReplaceAll(Form("_c%d",cateIndex), "");
    
    // SMOOTHING?
    if (m_options.Contains("Smooth") && !DMAnalysis::isDMSample(m_config, currName) && !currName.Contains("Data")) {
      currHist->Smooth(1);
    }
    
    // Stack plot filling and adding:
    if (!DMAnalysis::isDMSample(m_config, currName) && 
	!currName.Contains("Data")) {
      double scaleFactor = 1.0;
      if (m_options.Contains("Normalize")) {
	scaleFactor = 1.0 / nonDMSum_ALL;
      }
      else if (m_options.Contains("Scale2Data")) {
	if (allEvents) scaleFactor = nData_ALL / nonDMSum_ALL;
	else scaleFactor = nData_PASS / nonDMSum_PASS;
      }
      else {
	// Scale to appropriate luminosity
	scaleFactor = 0.001 * m_config->getNum("analysisLuminosity");
      }
      currHist->Scale(scaleFactor);
      currHist->SetFillColor(colorStack[index]);
      currHist->SetLineWidth(2);
      hs.Add(currHist);
      leg.AddEntry(currHist, 
		   DMAnalysis::getPrintSampleName(m_config, currName), "F");
    }
    // DM sample:
    else if (DMAnalysis::isDMSample(m_config, currName)) {
      currHist->SetLineColor(colorListDM[indexDM]);
      currHist->SetLineStyle(2+indexDM);
      indexDM++;
      leg.AddEntry(currHist, DMAnalysis::getPrintSampleName(m_config, currName),
		   "L");
      double scaleFactor = 1.0;
      if (m_options.Contains("Normalize")) {
	if (varName.Contains("cutFlowFull")) {
	  scaleFactor = 1.0 / currHist->GetBinContent(3);
	}
	else scaleFactor = 1.0 / currHist->Integral();
      }
      else {
	// Scale to appropriate luminosity
	scaleFactor = 0.001 * m_config->getNum("analysisLuminosity");
      }
      currHist->Scale(scaleFactor);
    }
    // Data:
    else {
      leg.AddEntry(currHist, "Data", "lep");
    }
    
    // Format plot axes:
    if (!varName.Contains("cutFlowFull")) {
      currHist->GetXaxis()->SetTitle(DMAnalysis::getPrintVarName(varName));
    }
    else {
      currHist->GetXaxis()->SetRangeUser(2,currHist->GetNbinsX());
    }
    int nGeV = (int)((currHist->GetXaxis()->GetXmax() - 
		      currHist->GetXaxis()->GetXmin()) / currHist->GetNbinsX());
    if (((TString)currHist->GetXaxis()->GetTitle()).Contains("GeV")) {
      currHist->GetYaxis()->SetTitle(Form("Entries / %d GeV",nGeV));
    }
    else {
      currHist->GetYaxis()->SetTitle("Entries");
    }
    
    if (m_options.Contains("LogScale")) {
      if (m_options.Contains("Normalize")) {
	currHist->GetYaxis()->SetRangeUser(0.00000001, 1000);
	currHist->GetYaxis()->SetTitle("Fractional Acceptance");
      }
      else if (allEvents) {
	currHist->GetYaxis()->SetRangeUser(0.01, 100*nonDMSum_ALL);
      }
      else currHist->GetYaxis()->SetRangeUser(0.01, 100*nonDMSum_PASS);
      gPad->SetLogy();
    }
    else {
      if (m_options.Contains("Normalize")) {
	currHist->GetYaxis()->SetRangeUser(0.0, 2.0);
      }
      else if (allEvents) currHist->GetYaxis()->SetRangeUser(0.0, nonDMSum_ALL);
      else currHist->GetYaxis()->SetRangeUser(0.0, nonDMSum_PASS);
    }
    
    if (varName.Contains("njets") || varName.Contains("nleptons")) {
      currHist->GetXaxis()->SetNdivisions(currHist->GetNbinsX());
      currHist->GetXaxis()->CenterLabels();
    }
    
    // Draw the histograms:
    if (index == 0) currHist->Draw();
    else if (!m_options.Contains("StackPlot") && !currName.Contains("Data")) {
      currHist->Draw("SAME");
    }
    
    index++;
  }
  
  // Draw the stack plots:
  if (m_options.Contains("StackPlot")) hs.Draw("SAMEHIST");
    
  //Then re-draw DM signals because they might not have been drawn on top:
  for (int i_s = 0; i_s < (int)m_hists.size(); i_s++) {
    TString currName = m_names[i_s];
    if (!currName.Contains(endTag)) continue;
    else if (cateIndex >= 0 && !currName.Contains(Form("c%d",cateIndex))) {
      continue;
    }
    else if (cateIndex < 0) {
      bool badSample = false;
      for (int i_c = 0; i_c < m_config->getStr("nCategories"); i_c++) {
	if (currName.Contains(Form("c%d",i_c))) {
	  badSample = true;
	  break;
	}
      }
      if (badSample) continue;
    }
    TH1F *currHist = m_hists[i_s];
    currName = currName.ReplaceAll("_"+varName, "");
    currName = currName.ReplaceAll(endTag, "");
    currName = currName.ReplaceAll(Form("_c%d",cateIndex), "");
    if (DMAnalysis::isDMSample(m_config, currName)) {
      currHist->Draw("axisSAME");
      currHist->Draw("SAME");
    }
    else if (currName.Contains("Data")) {
      currHist->Draw("EPSAME");
    }
  }
  
  // Also draw the ATLAS text:
  TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
  l.SetTextFont(72); l.SetTextSize(0.05); l.DrawLatex(0.18,0.88,"ATLAS");
  l.SetTextFont(42); l.SetTextSize(0.05); l.DrawLatex(0.3,0.88,"Internal");
  l.DrawLatex(0.18, 0.83, Form("#scale[0.8]{#sqrt{s} = 13 TeV: #scale[0.7]{#int}Ldt = %2.1f fb^{-1}}",(m_config->getNum("analysisLuminosity")/1000.0)));
  
  // Draw the legend then print the canvas:
  leg.Draw("SAME");
  if (cateIndex < 0) {
    can->Print(Form("%s/plot_%s%s.eps", m_outputDir.Data(), varName.Data(), 
		    endTag.Data()));
  }
  else {
    can->Print(Form("%s/plot_%s_c%d%s.eps", m_outputDir.Data(), varName.Data(), 
		    cateIndex, endTag.Data()));
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Sort the stored histograms from fewest to most entries, so that the stack 
   plots look OK. The algorithm here is O(N^2), but it was easy to code with 
   vectors... Also, we don't have that many samples to sort ;).
*/
void sortHistograms() {
  std::vector<TH1F*> m_hists_sorted; m_hists_sorted.clear();
  std::vector<TString> m_names_sorted; m_names_sorted.clear();
  
  // Sort the integral values:
  std::sort(m_integrals.begin(), m_integrals.end());
  
  // Loop over the integral values:
  for (int i_i = 0; i_i < (int)m_integrals.size(); i_i++) {
    
    // Loop to find the matching histogram:
    for (int i_h = 0; i_h < (int)m_hists.size(); i_h++) {
      
      if (fabs(m_hists[i_h]->Integral() - m_integrals[i_i]) < 0.00000001 && 
	  !isMatched(m_names[i_h])) {
	m_hists_sorted.push_back(m_hists[i_h]);
	m_names_sorted.push_back(m_names[i_h]);
	m_matched.push_back(m_names[i_h]);
      }
    }
  }
  
  m_hists = m_hists_sorted;
  m_names = m_names_sorted;
}

/**
   -----------------------------------------------------------------------------
   Produces plots of important analysis variables with signal and background MC,
   and data when available.
   @param configFile - The analysis configuration file.
   @param varName - The name of the variable in the file.
   @param options - The plot options.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 4) { 
    std::cout << "\nUsage: " << argv[0] 
	      << " <configFile> <varName> <options>" << std::endl;
    exit(0);
  }
  TString configFile = argv[1];
  TString varName = argv[2];
  m_options = argv[3];
  
  // Load the config file:
  m_config = new Config(configFile);
  
  // Count the number of non-DM events for stack normalization:
  nData_ALL = 0.0;
  nData_PASS = 0.0;
  nonDMSum_ALL = 0.0;
  nonDMSum_PASS = 0.0;
  
  m_hists.clear();
  m_names.clear();
  m_integrals.clear();
  m_matched.clear();
  
  // Set the ATLAS Style:
  CommonFunc::SetAtlasStyle();

  m_outputDir = Form("%s/%s/PlotVariables", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
    
  // Loop over the SM samples, loading histograms:
  std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    loadSampleHistograms(sigSMModes[i_SM], varName);
  }

  // Loop over the Bkg samples, loading histograms:
  std::vector<TString> bkgProcesses = m_config->getStrV("BkgProcesses");
  for (int i_b = bkgProcesses.size()-1; i_b >= 0; i_b--) {
    loadSampleHistograms(bkgProcesses[i_b], varName);
  }

  // Sort the histograms (fewest to most entries) BEFORE DATA OR SIGNALS:
  sortHistograms();
  
  // Loop over the DM samples, loading histograms:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    loadSampleHistograms(sigDMModes[i_DM], varName);
  }
  
  // Load the data:
  if (!m_config->getBool("doBlind")) {
    loadSampleHistograms("Data", varName);
  }
    
  // Make the plots without cuts applied and with cuts applied.
  // Avoid this if the variable being plotted is the cutflow:
  makeCombinedPlot(true, varName, -1);
  if (!varName.Contains("cutFlowFull")) {
    makeCombinedPlot(false, varName, -1);
  
    for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
      makeCombinedPlot(false, varName, i_c);
    }
  }

  return 0;
}
