////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMCutflowPlotter.cxx                                                //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
//  This macro plots the cutflows of several processes together.              //
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
#include "DMxAODCutflow.h"

// A map to store all of the histograms:
std::map<TString,TH1F*> m_hists;
Config *m_config;
TString m_outputDir;
double nonDMSum;

/**
   -----------------------------------------------------------------------------
   Load the histograms for a single sample and variable from file.
   @param sampleName - The name of the sample.
*/
void loadSampleHistograms(TString sampleName) {
  std::cout << "PlotVariables: loadSampleHistograms(" << sampleName << ")..."
	    << std::endl;
  
  // Tool to get the total number of events at the generator level.
  DMxAODCutflow *dmx
    = new DMxAODCutflow(DMAnalysis::nameToxAODCutFile(m_config, sampleName));
  
  if (sampleName.Contains("gjet")) {
    if (!m_hists["gjet"]) m_hists["gjet"] = dmx->getHist();
    else m_hists["gjet"]->Add(dmx->getHist());
  }
  
  else if (DMAnalysis::isSMSample(m_config, sampleName)) {
    if (!m_hists["SMHiggs"]) m_hists["SMHiggs"] = dmx->getHist();
    else m_hists["SMHiggs"]->Add(dmx->getHist());
  }
  else m_hists[sampleName] = dmx->getHist();
  
  if (!DMAnalysis::isDMSample(m_config, sampleName)) {
    nonDMSum += dmx->getHist()->GetBinContent(1);//Integral();
  }
  std::cout << "PlotVariables: Finished loading histograms from file."
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Makes a combined plot with all samples drawn for a given variable.
   @param options - The plot options.
*/
void makeCombinedPlot(TString options) {
  std::cout << "PlotVariables: makeCombinedPlot(" << options << ")" 
	    << std::endl;
  
  // Begin plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  TLegend leg(0.59,0.67,0.92,0.91);
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
    std::cout << "\tSample index " << index << ", " << histIter->first 
	      << std::endl;
    // Get pointer to the histogram and modify the name:
    TH1F *currHist = histIter->second;
    std::cout << "\tThe current name is " << currName << std::endl;
    
    if (options.Contains("StackPlot") && 
	!DMAnalysis::isDMSample(m_config, currName)) {
      currHist->Scale(100.0/nonDMSum);
      currHist->SetFillColor(colorStack[index]);
      currHist->SetLineWidth(2);
      hs.Add(currHist);
      leg.AddEntry(currHist, DMAnalysis::getPrintSampleName(m_config, currName),
		   "F");
    }
    else {
      //currHist->Scale(1.0/currHist->Integral());
      currHist->Scale(100.0/currHist->GetBinContent(1));
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
    currHist->GetYaxis()->SetTitle("Efficiency [%]");
    if (options.Contains("LogScale")) {
      currHist->GetYaxis()->SetRangeUser(0.001, 10.0);
      gPad->SetLogy();
    }
    else {
      currHist->GetYaxis()->SetRangeUser(0.0, 150.0);
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
    TH1F *currHist = histIter->second;
    if (DMAnalysis::isDMSample(m_config, currName)) {
      currHist->Draw("axisSAME");
      currHist->Draw("SAME");
    }
  }
  
  // Draw the legend then print the canvas:
  leg.Draw("SAME");
  can->Print(Form("%s/plot_cutflow.eps", m_outputDir.Data()));
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Produces comparison plots of the analysis cutflow for multiple processes.
   @param configFile - The analysis configuration file.
   @param options - The plot options.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0]
	      << " <configFile> <options>" << std::endl;
    exit(0);
  }
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Load the analysis configuration file:
  m_config = new Config(configFile);
  
  // Set the ATLAS Style:
  CommonFunc::SetAtlasStyle();

  m_outputDir = Form("%s/%s/PlotVariables", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  nonDMSum = 0;
  
  // Loop over the SM samples, loading histograms:
  std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    loadSampleHistograms(sigSMModes[i_SM]);
  }
  
  // Loop over the DM samples, loading histograms:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    loadSampleHistograms(sigDMModes[i_DM]);
  }

  // Loop over the Bkg samples, loading histograms:
  std::vector<TString> bkgProcesses = m_config->getStrV("BkgProcesses");
  for (int i_b = 0; i_b < (int)bkgProcesses.size(); i_b++) {
    if ((bkgProcesses[i_b]).Contains("gjet")) continue;
    loadSampleHistograms(bkgProcesses[i_b]);
  }
  
  // Make the plots without cuts applied and with cuts applied:
  makeCombinedPlot(options);
  return 0;
}
