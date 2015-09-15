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
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ libraries:
#include <iostream>
#include <fstream>
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
std::map<TString,TH1F*> m_hists;
Config *m_config;
TString m_outputDir;

double nonDMSum_ALL;
double nonDMSum_PASS;

/**
   -----------------------------------------------------------------------------
   Load the histograms for a single sample and variable from file.
   @param sampleName - The name of the sample.
   @param varName - The name of the variable in the file.
*/
void loadSampleHistograms(TString sampleName, TString varName) {
  std::cout << "PlotVariables: loadSampleHistograms(" << sampleName << ", " 
	    << varName << ")..." << std::endl;
  
  TString inputDir = Form("%s/%s/DMMassPoints", 
			  (m_config->getStr("masterOutput")).Data(), 
			  (m_config->getStr("jobName")).Data());
  TFile *inputFile = new TFile(Form("%s/hists_%s.root", inputDir.Data(),sampleName.Data()));
  if (sampleName.Contains("gjet")) {
    if (!m_hists[Form("gjet_%s_ALL", varName.Data())]) {
      m_hists[Form("gjet_%s_ALL", varName.Data())]
	= (TH1F*)inputFile->Get(Form("%s_ALL",varName.Data()));
    }
    else {
      m_hists[Form("gjet_%s_ALL", varName.Data())]
	->Add((TH1F*)inputFile->Get(Form("%s_ALL",varName.Data())), 1.0);
    }
    if (!m_hists[Form("gjet_%s_PASS", varName.Data())]) {
      m_hists[Form("gjet_%s_PASS", varName.Data())]
	= (TH1F*)inputFile->Get(Form("%s_PASS",varName.Data()));
      for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
	m_hists[Form("gjet_%s_c%d_PASS",varName.Data(),i_c)]
	  = (TH1F*)inputFile->Get(Form("%s_c%d_PASS",varName.Data(),i_c));
      }
    }
    else {
      m_hists[Form("gjet_%s_PASS",varName.Data())]
	->Add((TH1F*)inputFile->Get(Form("%s_PASS",varName.Data())), 1.0);
      for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
	m_hists[Form("gjet_%s_c%d_PASS",varName.Data(),i_c)]
	  ->Add((TH1F*)inputFile
		->Get(Form("%s_c%d_PASS",varName.Data(),i_c)),1.0);
      }
    }
  }
  else if (DMAnalysis::isSMSample(m_config, sampleName)) {
    if (!m_hists[Form("SMHiggs_%s_ALL", varName.Data())]) {
      m_hists[Form("SMHiggs_%s_ALL", varName.Data())] 
	= (TH1F*)inputFile->Get(Form("%s_ALL",varName.Data()));
    }
    else {
      m_hists[Form("SMHiggs_%s_ALL", varName.Data())] 
	->Add((TH1F*)inputFile->Get(Form("%s_ALL",varName.Data())));
    }
    if (!m_hists[Form("SMHiggs_%s_PASS", varName.Data())]) {
      m_hists[Form("SMHiggs_%s_PASS", varName.Data())] 
	= (TH1F*)inputFile->Get(Form("%s_PASS",varName.Data()));
      for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
	m_hists[Form("SMHiggs_%s_c%d_PASS",varName.Data(),i_c)] 
	  = (TH1F*)inputFile->Get(Form("%s_c%d_PASS",varName.Data(),i_c));
      }
    }
    else {
      m_hists[Form("SMHiggs_%s_PASS", varName.Data())] 
	->Add((TH1F*)inputFile->Get(Form("%s_PASS",varName.Data())));
      for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
	m_hists[Form("SMHiggs_%s_c%d_PASS",varName.Data(),i_c)] 
	  ->Add((TH1F*)inputFile->Get(Form("%s_c%d_PASS",varName.Data(),i_c)));
      }
    }
  }

  else {
    m_hists[Form("%s_%s_ALL", sampleName.Data(), varName.Data())] 
      = (TH1F*)inputFile->Get(Form("%s_ALL",varName.Data()));
    m_hists[Form("%s_%s_PASS", sampleName.Data(), varName.Data())] 
      = (TH1F*)inputFile->Get(Form("%s_PASS",varName.Data()));
    for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
      m_hists[Form("%s_%s_c%d_PASS", sampleName.Data(), varName.Data(), i_c)] 
	= (TH1F*)inputFile->Get(Form("%s_c%d_PASS",varName.Data(),i_c));
    }
  }
  
  if (!DMAnalysis::isDMSample(m_config, sampleName)) {
    nonDMSum_ALL
      += ((TH1F*)inputFile->Get(Form("%s_ALL",varName.Data())))->Integral();
    nonDMSum_PASS 
      += ((TH1F*)inputFile->Get(Form("%s_PASS",varName.Data())))->Integral();
  }
  
  std::cout << "PlotVariables: Finished loading histograms from file."
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Makes a combined plot with all samples drawn for a given variable.
   @param allEvents - True iff all events are in the hist (not just selected).
   @param varName - The name of the variable in the file.
   @param options - The plot options.
*/
void makeCombinedPlot(bool allEvents, TString varName, TString options,
		      int cateIndex) {
  std::cout << "PlotVariables: makeCombinedPlot(" << allEvents << ", "
	    << varName << ", " << options << ", " << cateIndex << ")"
	    << std::endl;
  
  TString endTag = allEvents ? "_ALL" : "_PASS";
  
  // Begin plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  TLegend leg(0.6,0.67,0.9,0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.03);
  
  // For stacked histogram:
  THStack hs("hs", "stacked histogram");
  
  // Color palettes for SM Higgs, DM Higgs, and Bkg:
  Color_t colorListDM[6] = {kRed, kRed+3, kRed-2, kRed-7, kRed+4, kRed+1};
  Color_t colorListSM[6] = {kBlue,kBlue-1,kBlue+3,kBlue-3,kBlue+5,kBlue-5};
  Color_t colorListBkg[2] = {kGreen-2, kGreen+2};
  
  Color_t colorStack[8] = {kViolet+6, kTeal-9, kAzure+2, kBlue-1, 
			   kViolet-4, kTeal+2, kAzure-8, kBlue+1};
  
  // Loop over the saved histograms:
  int index = 0; int indexSM = 0; int indexDM = 0; int indexBkg = 0;
  std::map<TString,TH1F*>::iterator histIter;
  for (histIter = m_hists.begin(); histIter != m_hists.end(); histIter++) {
    TString currName = histIter->first;
    
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
    std::cout << "\tSample index " << index << ", " << histIter->first 
	      << std::endl;
    // Get pointer to the histogram and modify the name:
    TH1F *currHist = histIter->second;  
    currName = currName.ReplaceAll("_"+varName, "");
    currName = currName.ReplaceAll(endTag, "");
    currName = currName.ReplaceAll(Form("_c%d",cateIndex), "");
    std::cout << "\tThe current name is " << currName << std::endl;
    
    if (options.Contains("StackPlot") && 
	!DMAnalysis::isDMSample(m_config, currName)) {
      if (allEvents) currHist->Scale(1.0/nonDMSum_ALL);
      else currHist->Scale(1.0/nonDMSum_PASS);
      currHist->SetFillColor(colorStack[index]);
      currHist->SetLineWidth(2);
      hs.Add(currHist);
      leg.AddEntry(currHist, DMAnalysis::getPrintSampleName(m_config, currName),
		   "F");
    }
    else {
      currHist->Scale(1.0/currHist->Integral());
      currHist->SetLineWidth(2);
      if (DMAnalysis::isSMSample(m_config, currName) || 
	  currName.Contains("SMHiggs")) {
	currHist->SetLineColor(colorListSM[indexSM]);
	indexSM++;
      }
      else if (!DMAnalysis::isDMSample(m_config, currName)) {
	currHist->SetLineColor(colorListBkg[indexBkg]);
	indexBkg++;
      }
      leg.AddEntry(currHist, DMAnalysis::getPrintSampleName(m_config, currName),
		   "L");
    }
    
    if (DMAnalysis::isDMSample(m_config, currName)) {
      currHist->SetLineColor(colorListDM[indexDM]);
      currHist->SetLineStyle(1+indexDM);
      indexDM++;
    }

    // Format plot axes:
    currHist->GetXaxis()->SetTitle(DMAnalysis::getPrintVarName(varName));
    currHist->GetYaxis()->SetTitle("Normalized Entries");
    if (options.Contains("LogScale")) {
      currHist->GetYaxis()->SetRangeUser(0.001, 10.0);
      gPad->SetLogy();
    }
    else {
      currHist->GetYaxis()->SetRangeUser(0.0, 0.3);
    }
    
    // Draw the histograms:
    if (index == 0) currHist->Draw();
    else if (!options.Contains("StackPlot")) {
      currHist->Draw("SAME");
    }
    index++;
  }
  
  // Draw the stack plots:
  if (options.Contains("StackPlot")) hs.Draw("SAMEHIST");
  
  //Then re-draw DM signals because they might not have been drawn on top:
  for (histIter = m_hists.begin(); histIter != m_hists.end(); histIter++) {
    TString currName = histIter->first;
    if (!(histIter->first).Contains(endTag)) continue;
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
    TH1F *currHist = histIter->second;    
    currName = currName.ReplaceAll("_"+varName, "");
    currName = currName.ReplaceAll(endTag, "");
    currName = currName.ReplaceAll(Form("_c%d",cateIndex), "");
    if (DMAnalysis::isDMSample(m_config, currName)) {
      currHist->Draw("axisSAME");
      currHist->Draw("SAME");
    }
  }
  
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
  TString options = argv[3];
  
  // Load the config file:
  m_config = new Config(configFile);
  
  // Count the number of non-DM events for stack normalization:
  double nonDMSum_ALL = 0.0;
  double nonDMSum_PASS = 0.0;
  
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
  
  // Loop over the DM samples, loading histograms:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    loadSampleHistograms(sigDMModes[i_DM], varName);
  }

  // Loop over the Bkg samples, loading histograms:
  std::vector<TString> bkgProcesses = m_config->getStrV("BkgProcesses");
  for (int i_b = 0; i_b < (int)bkgProcesses.size(); i_b++) {
    if ((bkgProcesses[i_b]).Contains("gjet")) continue;
    loadSampleHistograms(bkgProcesses[i_b], varName);
  }
    
  // Make the plots without cuts applied and with cuts applied:
  makeCombinedPlot(true, varName, options, -1);
  makeCombinedPlot(false, varName, options, -1);
  
  for (int i_c = 0; i_c < m_config->getInt("nCategories"); i_c++) {
    makeCombinedPlot(false, varName, options, i_c);
  }

  return 0;
}
