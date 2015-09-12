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
//std::map<TString,TFile*> m_files;
Config *m_config;
TString m_outputDir;

/**
   -----------------------------------------------------------------------------
   Get a pretty LaTex formatted name for a variable.
   @param varName - The name of the variable in the file.
   @returns - A LaTex formatted TString.
*/
TString getPrintVarName(TString varName) {
  if (varName.EqualTo("pTyy")) return "p_{T}^{#gamma#gamma} [GeV]";
  else if (varName.EqualTo("ETMiss")) return "#slash{E}_{T} [GeV]";
  else if (varName.EqualTo("ratioETMisspTyy")) return "#slash{E}_{T}/p_{T}^{#gamma#gamma}";
  else if (varName.EqualTo("myy")) return "M_{#gamma#gamma} [GeV]";
  else return varName;
}

/**
   -----------------------------------------------------------------------------
   Get a pretty LaTex formatted name for a variable.
   @param sampleName - The name of the sample.
   @returns - A LaTex formatted TString.
*/
TString getPrintSampleName(TString sampleName) {
  if (sampleName.EqualTo("gg")) return "#gamma#gamma";
  else if (sampleName.EqualTo("gjet")) return "#gamma+jet";
  else if (sampleName.EqualTo("SMHiggs")) return "SM H#rightarrow#gamma#gamma";
  else if (DMAnalysis::isDMSample(m_config, sampleName)) {
    TString mediatorName = DMAnalysis::getMediatorName(sampleName);
    int mediatorMass = DMAnalysis::getMediatorMass(m_config, sampleName);
    TString medMassForm;
    if (mediatorName.EqualTo("shxx_gg")) {
      medMassForm = Form("m_{S}=%dGeV", mediatorMass);
    }
    else if (mediatorName.EqualTo("zphxx_gg")) {
      medMassForm = Form("m_{Z'}=%dGeV", mediatorMass);
    }
    else if (mediatorName.EqualTo("zp2hdmxx_gg")) {
      medMassForm = Form("m_{Z'}=%dGeV", mediatorMass);
    }
    
    int darkMatterMass = DMAnalysis::getDarkMatterMass(m_config, sampleName);
    TString darkMassForm;
    if (mediatorName.EqualTo("shxx_gg")) {
      darkMassForm = Form("m_{#chi}=%dGeV", darkMatterMass);
    }
    else if (mediatorName.EqualTo("zphxx_gg")) {
      darkMassForm = Form("m_{#chi}=%dGeV", darkMatterMass);
    }
    else if (mediatorName.EqualTo("zp2hdmxx_gg")) {
      darkMassForm = Form("m_{A}=%dGeV", darkMatterMass);
    }
    //(DMAnalysis::getPrintMediatorName(sampleName)).Data(),
    TString newName = Form("%s %s", medMassForm.Data(), darkMassForm.Data());
    return newName;
  }
  else return sampleName;
}


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
    }
    else {
      m_hists[Form("gjet_%s_PASS", varName.Data())]
	->Add((TH1F*)inputFile->Get(Form("%s_PASS",varName.Data())), 1.0);
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
    }
    else {
      m_hists[Form("SMHiggs_%s_PASS", varName.Data())] 
	->Add((TH1F*)inputFile->Get(Form("%s_PASS",varName.Data())));
    }
  }

  else {
    m_hists[Form("%s_%s_ALL", sampleName.Data(), varName.Data())] 
      = (TH1F*)inputFile->Get(Form("%s_ALL",varName.Data()));
    m_hists[Form("%s_%s_PASS", sampleName.Data(), varName.Data())] 
      = (TH1F*)inputFile->Get(Form("%s_PASS",varName.Data()));
  }
  std::cout << "PlotVariables: Finished loading histograms from file."
	    << std::endl;
  //inputFile.Close();
}

/**
   -----------------------------------------------------------------------------
   Makes a combined plot with all samples drawn for a given variable.
   @param allEvents - True iff all events are in the hist (not just selected).
   @param varName - The name of the variable in the file.
*/
void makeCombinedPlot(bool allEvents, TString varName) {
  std::cout << "PlotVariables: makeCombinedPlot(" << allEvents << ", "
	    << varName << ")..." << std::endl;
  
  TString endTag = allEvents ? "_ALL" : "_PASS";
  
  // Begin plotting:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  TLegend leg(0.6,0.7,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.03);
  
  // Color palettes for SM Higgs, DM Higgs, and Bkg:
  Color_t colorListSM[6] = {kRed+1, kRed-1, kRed+3, kRed-3, kRed+5, kRed-5};
  Color_t colorListDM[5] = {kBlue, kBlue+2, kBlue-2, kBlue+4, kBlue-4};
  Color_t colorListBkg[2] = {kGreen-2, kGreen+2};
  
  // Loop over the saved histograms:
  int index = 0; int indexSM = 0; int indexDM = 0; int indexBkg = 0;
  std::map<TString,TH1F*>::iterator histIter;
  for (histIter = m_hists.begin(); histIter != m_hists.end(); histIter++) {
    if (!(histIter->first).Contains(endTag)) continue;
    std::cout << "\tSample index " << index << ", " << histIter->first 
	      << std::endl;
    TH1F *currHist = histIter->second;
    currHist->Scale(1.0/currHist->Integral());
    currHist->GetXaxis()->SetTitle(getPrintVarName(varName));
    currHist->GetYaxis()->SetTitle("Normalized Entries");
    currHist->SetLineWidth(2);
    TString currName = histIter->first;
    currName = currName.ReplaceAll("_"+varName, "");
    currName = currName.ReplaceAll(endTag, "");
    std::cout << "\tThe current name is " << currName << std::endl;
    if (DMAnalysis::isSMSample(m_config, currName) || 
	currName.Contains("SMHiggs")) {
      std::cout << "\t  -> SM Higgs" << std::endl;
      currHist->SetLineColor(colorListSM[indexSM]);
      indexSM++;
    }
    
    else if (DMAnalysis::isDMSample(m_config, currName)) {
      std::cout << "\t  -> DM Higgs" << std::endl;
      currHist->SetLineColor(colorListDM[indexDM]);
      indexDM++;
    }
    else {
      std::cout << "\t  -> Bkg" << std::endl;
      currHist->SetLineColor(colorListBkg[indexBkg]);
      indexBkg++;
    }
    
    leg.AddEntry(currHist, getPrintSampleName(currName), "L");
    if (index == 0) currHist->Draw();
    else currHist->Draw("SAME");
    index++;
  }
  leg.Draw("SAME");
  can->Print(Form("%s/plot_%s%s.eps", m_outputDir.Data(), varName.Data(), 
		  endTag.Data()));
  delete can;
}

/**
   -----------------------------------------------------------------------------
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
  
  // Set the ATLAS Style:
  CommonFunc::SetAtlasStyle();

  m_outputDir = Form("%s/%s/PlotVariables", 
		     (m_config->getStr("masterOutput")).Data(), 
		     (m_config->getStr("jobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Loop over the samples, loading histograms:
  std::vector<TString> sigSMModes = m_config->getStrV("sigSMModes");
  for (int i_SM = 0; i_SM < (int)sigSMModes.size(); i_SM++) {
    loadSampleHistograms(sigSMModes[i_SM], varName);
  }
  
  // Loop over the samples, loading histograms:
  std::vector<TString> sigDMModes = m_config->getStrV("sigDMModes");
  for (int i_DM = 0; i_DM < (int)sigDMModes.size(); i_DM++) {
    loadSampleHistograms(sigDMModes[i_DM], varName);
  }

  // Loop over the samples, loading histograms:
  std::vector<TString> bkgProcesses = m_config->getStrV("BkgProcesses");
  for (int i_b = 0; i_b < (int)bkgProcesses.size(); i_b++) {
    if ((bkgProcesses[i_b]).Contains("gjet")) continue;
    loadSampleHistograms(bkgProcesses[i_b], varName);
  }
    
  // Make the plots without cuts applied and with cuts applied:
  makeCombinedPlot(true, varName);
  makeCombinedPlot(false, varName);
  
  return 0;
}
