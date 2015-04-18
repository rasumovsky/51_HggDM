////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMtoyMC.cxx                                                         //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/04/2015                                                          //
//                                                                            //
//  This program creates pseudoexperiment ensembles for the H->yy + DM search //
//  with 13 TeV data at ATLAS.                                                //
//                                                                            //
//  options:                                                                  //
//      binned, FromSnapshot, FixMu, scan, FixBkgParam, nuisplot              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMToyMC.hh"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;
using namespace DMAnalysis;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PlotNuisParams:

void PlotNuisParams(RooArgSet *nuis) {
  can->cd();
  can->SetLeftMargin(0.3);
  can->SetRightMargin(0.1);
  gPad->SetLogy(0);
  gStyle->SetHistMinimumZero();
  int index = 0;
  int number_params = nuis->getSize();
  TH1F *h_nuis = new TH1F("h_nuis","h_nuis",number_params,0,number_params);
  TIterator *iter_nuis = nuis->createIterator();
  RooRealVar* parg_nuis = NULL;
  while ((parg_nuis = (RooRealVar*)iter_nuis->Next())) {
    TString name = parg_nuis->GetName();
    TString new_name = name;
    new_name.ReplaceAll("alpha_","");
    new_name.ReplaceAll("norm_","");
    double value = parg_nuis->getVal();
    if (name.Contains("alpha")) {
      index++;
      h_nuis->SetBinContent(index, value);
      h_nuis->GetXaxis()->SetBinLabel(index, new_name);
    }
  }
  h_nuis->SetFillColor(38);
  h_nuis->SetBarWidth(0.6);
  h_nuis->SetBarOffset(0.2);
  h_nuis->SetStats(0);
  h_nuis->GetYaxis()->SetRangeUser(-2.3,2.3);
  h_nuis->GetYaxis()->SetNdivisions(6);
  h_nuis->GetYaxis()->SetTitle("Nuisance Parameter Pull (#sigma)");
  h_nuis->GetXaxis()->SetRangeUser(0, index);
  h_nuis->Draw("hbar");
  
  TLine *l0 = new TLine(0, 0, 0, index);
  l0->SetLineColor(kBlack);
  l0->SetLineWidth(1);
  l0->SetLineStyle(1);
  l0->Draw("SAME");
  TLine *lp1 = new TLine(1, 0, 1, index);
  lp1->SetLineColor(kBlack);
  lp1->SetLineWidth(1);
  lp1->SetLineStyle(2);
  lp1->Draw("SAME");
  TLine *lp2 = new TLine(2, 0, 2, index);
  lp2->SetLineColor(kBlack);
  lp2->SetLineWidth(1);
  lp2->SetLineStyle(2);
  lp2->Draw("SAME");
  TLine *ln1 = new TLine(-1, 0, -1, index);
  ln1->SetLineColor(kBlack);
  ln1->SetLineWidth(1);
  ln1->SetLineStyle(2);
  ln1->Draw("SAME");
  TLine *ln2 = new TLine(-2, 0, -2, index);
  ln2->SetLineColor(kBlack);
  ln2->SetLineWidth(1);
  ln2->SetLineStyle(2);
  ln2->Draw("SAME");
  can->Print(Form("%s/plot_nuisparams_%s.eps", outputDir.Data(),
		  signalProcess.Data()));
  can->Clear();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// generateToyData:

RooAbsData* generateToyData(RooWorkspace *w, ModelConfig *mc, int seed,
			    bool m_binned = false) {
  
  std::cout << "NPPV1_hf_toyMC::generateToyData( binned()" << std::endl;
  
  RooSimultaneous* combPdf =  (RooSimultaneous*)mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  RooArgSet* Observables = (RooArgSet*)mc->GetObservables();
  RooArgSet* originValsNP = (RooArgSet*)mc->GetNuisanceParameters()->snapshot();
  RooRealVar* firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
  
  map<string,RooDataSet*> toydatasetMap; 
  RooCategory *channellist = (RooCategory*)w->obj("channelCat");
  TIterator* iter =  combPdf->indexCat().typeIterator();
  RooCatType* tt = NULL;
  RooDataSet* datatmp[20];
  
  // loop over all channels
  int idx = 0;
  
  // WARNING! THIS WAS PREVIOUSLY COMMENTED, WHILE SIMILAR LINE BELOW WAS UNCOMMENTED.
  statistics::randomizeSet( combPdf, globalObservables, seed ); 
  statistics::constSet( globalObservables, true );
  
  cout << "Printing the globalObservables #1!" << endl;
  globalObservables->Print("v");
  cout << "Printing the combined PDF." << endl;
  combPdf->Print();
  
  nToyEvents.clear();
  
  while( (tt=(RooCatType*) iter->Next()) )
  {
    cout << idx << " " << tt->GetName() << endl;
    RooAbsPdf* pdftmp = combPdf->getPdf(tt->GetName());
    RooArgSet* obstmp = pdftmp->getObservables( Observables );
    // this is very important. this is fucking mandatory
    RooArgSet* glotmp = pdftmp->getObservables( globalObservables );
    RooRealVar *t = (RooRealVar*)obstmp->first();
    
    cout << "    Here's exactly what's being fitted: " << endl;
    cout << "    all observables: " << endl;
    obstmp->Print("v");
    cout << "    single observable: " << endl;
    t->Print("v");
    cout << "    pdftmp: " << endl;
    pdftmp->Print("v");
    
    //statistics::randomizeSet( pdftmp, glotmp, -1 );
    //statistics::constSet( glotmp, true );
    
    
    if( m_binned )
    {
      pdftmp->setAttribute("PleaseGenerateBinned");
      TIterator *iter_obs = obstmp->createIterator();
      RooRealVar* parg_obs = NULL;
      // Was at more bins...
      while( (parg_obs=(RooRealVar*)iter_obs->Next()) ) { parg_obs->setBins(NbinsT); cout << "obstmp names special : " << parg_obs->GetName() << endl; }
      datatmp[idx] = (RooDataSet*)pdftmp->generate( *obstmp, AutoBinned(true), Extended(pdftmp->canBeExtended()), GenBinned("PleaseGenerateBinned") );
    }
    else datatmp[idx] = (RooDataSet*)pdftmp->generate( *obstmp, Extended(true) );
    toydatasetMap[(string)tt->GetName()] = datatmp[idx];
    nToyEvents.push_back((double)datatmp[idx]->sumEntries());
    idx++;
    
    /*
     // Finally, print the data and PDF for this category...
    if( isFirstToy )
    {
      can->cd();
      can->Clear();
      RooPlot* frame = t->frame(6);
      (toydatasetMap[(string)tt->GetName()])->plotOn( frame, LineColor(kBlack), MarkerColor(kBlack) );
      pdftmp->plotOn( frame, LineColor(kBlue), LineWidth(2) );
      frame->GetYaxis()->SetTitle("Events / bin");
      frame->GetXaxis()->SetTitle("t_{#gamma} bin");
      frame->Draw();
      TLatex cs; cs.SetNDC(); cs.SetTextColor(1); cs.SetTextFont(42); cs.DrawLatex(0.55,0.85,Form("|z_{DCA}| category %i",idx));
      TString catename = (TString)tt->GetName();
      can->Print(Form("%s/plots/plot_toy_mu%i_%s_%i.eps",outputDir.Data(),inputMuDM,catename.Data(),seed));
    }
    */
    
  }
  
  cout << "Printing the globalObservables #2!" << endl;
  globalObservables->Print("v");
  
  RooArgSet* args = new RooArgSet();
  RooRealVar *wt = w->var("weightVar");
  args->add(*Observables);
  args->add(*wt);
  RooDataSet* toyData = new RooDataSet( "toyData", "toyData", *args, Index(*channellist), Import(toydatasetMap), WeightVar(*wt) );

  toyData->Print();
  // release nuisance parameters:
  statistics::constSet( nuisanceParameters, false );
  return toyData;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// calcLHCCTestStat:

map<string,double> calcLHCTestStat( ModelConfig *mc, RooAbsData* data, int& status_mu1, int& status_mu0, int& status_MuFree )
{
  cout << "NPPV1_hf_toyMC::calcLHCTestStat()" << endl;
  nameNP.clear();
  valueNPMu1.clear();
  valueNPMu0.clear();
  valueNPMuFree.clear();
  
  nameGlobs.clear();
  valueGlobsMu1.clear();
  valueGlobsMu0.clear();
  valueGlobsMuFree.clear();
  
  cout << "Entering function calcLHCTestStat()" << endl;
  statistics::setDefaultPrintLevel(0);
  RooAbsPdf* combPdf = mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  RooArgSet* Observables = (RooArgSet*)mc->GetObservables();
  RooArgSet* originValsNP = (RooArgSet*)mc->GetNuisanceParameters()->snapshot();
  RooArgSet* poi = (RooArgSet*)mc->GetParametersOfInterest();
  RooRealVar* firstpoi = (RooRealVar*)poi->first();
  
  // First, print the "original values" of the NP:
  cout << "Printing original values of the NP." << endl;
  originValsNP->Print("v");
  
  cout << "Here's the toy data under suspicion of being fake:" << endl;
  data->Print("v");
  cout << "data->sumEntries() = " << data->sumEntries() << endl;
  
  //----------------------------------------//
  // first fit at mu = 1.0
  cout << "Doing mu = 1:" << endl;
  
  // Prepare to fit toy data:
  statistics::constSet( nuisanceParameters, false );
  statistics::constSet( globalObservables, true );
  firstpoi->setVal(1.0);
  firstpoi->setConstant(true);
  
  RooNLLVar* nll_mu1 = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize( status_mu1, nll_mu1, "", NULL, false );
  cout << "  Step 1: Completed mu=1 fit" << endl;
  
  // iterate over nuisanceParameters:
  TIterator *iter_nuifix_mu1 = nuisanceParameters->createIterator();
  RooRealVar* parg_nuifix_mu1 = NULL;
  nuisanceParameters->Print("v");
  while( (parg_nuifix_mu1 = (RooRealVar*)iter_nuifix_mu1->Next()) )
  {
    nameNP.push_back(parg_nuifix_mu1->GetName());
    valueNPMu1.push_back(parg_nuifix_mu1->getVal());
  }
  delete iter_nuifix_mu1;
  
  // iterate over globalObservables:
  cout << "Iterate over globalObservables." << endl;
  TIterator *iter_globfix_mu1 = globalObservables->createIterator();
  RooRealVar* parg_globfix_mu1 = NULL;
  globalObservables->Print("v");
  while( (parg_globfix_mu1 = (RooRealVar*)iter_globfix_mu1->Next()) )
  {
    nameGlobs.push_back(parg_globfix_mu1->GetName());
    valueGlobsMu1.push_back(parg_globfix_mu1->getVal());
  }
  delete iter_globfix_mu1;
  
  // Save the NLL value:
  double NLL_mu1 = nll_mu1->getVal();
  delete nll_mu1;
 
  //----------------------------------------//
  // second fit mu = 0.0:
  cout << "Doing mu = 0:" << endl;  
  
  // release nuisance parameters and set to default values after fit:
  statistics::constSet( nuisanceParameters, false, originValsNP );
  statistics::constSet( globalObservables, true );
  firstpoi->setVal(0.0);
  firstpoi->setConstant(true);
  
  RooNLLVar* nll_mu0 = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize( status_mu0, nll_mu0, "", NULL, false );
  cout << "  Step 2: Completed mu=0 fit" << endl;
  
  // iterate over nuisanceParameters:
  TIterator *iter_nuifix_mu0 = nuisanceParameters->createIterator();
  RooRealVar* parg_nuifix_mu0 = NULL;
  while( (parg_nuifix_mu0=(RooRealVar*)iter_nuifix_mu0->Next()) )
    valueNPMu0.push_back(parg_nuifix_mu0->getVal());
  delete iter_nuifix_mu0;
  
  // iterate over globalObservables:
  cout << "Iterate over globalObservables." << endl;
  TIterator *iter_globfix_mu0 = globalObservables->createIterator();
  RooRealVar* parg_globfix_mu0 = NULL;
  globalObservables->Print("v");
  while( (parg_globfix_mu0 = (RooRealVar*)iter_globfix_mu0->Next()) )
    valueGlobsMu0.push_back(parg_globfix_mu0->getVal());
  delete iter_globfix_mu0;
  
  // Save the NLL value:
  double NLL_mu0 = nll_mu0->getVal();
  delete nll_mu0;
 
  //----------------------------------------//
  // third fit mu free:
  cout << "Doing mu free:" << endl;  
  
  // release nuisance parameters and set to default values after fit:
  statistics::constSet( nuisanceParameters, false, originValsNP );
  statistics::constSet( globalObservables, true );
  // BAD! Why was the line below present in the first place?!
  //statistics::constSet( globalObservables, false );// release global observables after fit
  firstpoi->setVal(0.0);
  firstpoi->setConstant(false);
  
  cout << "TEST OF NP BEFORE MUFREE FIT" << endl;
  nuisanceParameters->Print("v");
  cout << "TEST OF GLOBS BEFORE MUFREE FIT" << endl;
  globalObservables->Print("v");
  

  RooNLLVar* nll_MuFree = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  statistics::minimize( status_MuFree, nll_MuFree, "", NULL, false );
  cout << "  Step 3: Completed mu=free fit" << endl;
  double mu_profiled = firstpoi->getVal();
  
  cout << "TEST OF NP AFTER MUFREE FIT" << endl;
  nuisanceParameters->Print("v");
  cout << "TEST OF GLOBS AFTER MUFREE FIT" << endl;
  globalObservables->Print("v");
  
  
  // iterate over nuisanceParameters:
  TIterator *iter_nuifix_MuFree = nuisanceParameters->createIterator();
  RooRealVar* parg_nuifix_MuFree = NULL;
  while( (parg_nuifix_MuFree=(RooRealVar*)iter_nuifix_MuFree->Next()) )
    valueNPMuFree.push_back(parg_nuifix_MuFree->getVal());
  delete iter_nuifix_MuFree;
  
  // iterate over globalObservables:
  cout << "Iterate over globalObservables." << endl;
  TIterator *iter_globfix_MuFree = globalObservables->createIterator();
  RooRealVar* parg_globfix_MuFree = NULL;
  globalObservables->Print("v");
  while( (parg_globfix_MuFree = (RooRealVar*)iter_globfix_MuFree->Next()) )
    valueGlobsMuFree.push_back(parg_globfix_MuFree->getVal());
  delete iter_globfix_MuFree;

  // Save the NLL value:
  double NLL_MuFree = nll_MuFree->getVal();
  delete nll_MuFree;
  
  if( option.Contains("nuisplot") && firstPlot )
  {
    PlotNuisParams( nuisanceParameters );
    firstPlot = false;
  }
  
  // release nuisance parameters and set to default values after fit:
  statistics::constSet( nuisanceParameters, false, originValsNP );
  // WHY WAS THIS EVER ADDED: ?!
  statistics::constSet( globalObservables, false );// release global observables after fit
  
  map<string,double> result;
  result["profiled_mu"] = mu_profiled;
  result["nll_mu1"] = NLL_mu1; 
  result["nll_mu0"] = NLL_mu0;
  result["nll_MuFree"] = NLL_MuFree;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// ProfileToData:

void ProfileToData( ModelConfig *mc, RooAbsData* data )
{
  cout << "  NPPV1_hf_toyMC::ProfileToData() " << endl;
  RooAbsPdf* combPdf  =  (RooAbsPdf*)mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)mc->GetGlobalObservables();
  RooArgSet* Observables = (RooArgSet*)mc->GetObservables();
  RooRealVar* firstPOI=(RooRealVar*)mc->GetParametersOfInterest()->first();
  firstPOI->setConstant(false);

  TString pdfname=combPdf->GetName();
  
  // Profile to data
  // Fix nuisance parameters before randomize global observables.
  statistics::constSet(nuisanceParameters,false);
  // Release global observables before randomize global observables.
  statistics::constSet(globalObservables,true);  
  
  // set the mu = 0 if doing background fix:
  if( option.Contains("FixBkgParam") )
  {
    firstPOI->setVal(0.0);
    firstPOI->setConstant(true);
  }
  
  RooNLLVar* nll = (RooNLLVar*)combPdf->createNLL(*data, Constrain(*nuisanceParameters), Extended(combPdf->canBeExtended()));
  RooFitResult *res = statistics::minimize(nll);
  
  delete res;
  // Release nuisance parameters:
  statistics::constSet(nuisanceParameters,false);
  // Release global observables before randomize global observables.
  statistics::constSet(globalObservables,false); 
  firstPOI->setConstant(false);
}





/**
   The main method. 
   @param jobName - the name of the analysis job.
   @param option - the options (see header note).
   @param seed - the random seed for pseudoexperiment creation.
   @param toysPerJob - the number of pseudoexperiments to create per job.
   @param muDMVal - the value of the DM signal strength to use.
   @param signalProcess - the DM signal production process of interest. 
*/
int main(int argc , char **argv) {

  // Root macros:
  SetAtlasStyle();
  
  if (argc < 7) {
    std::cout << "Usage: " << argv[0] << "./bin/NPPV1_hf_toyMC <jobName> <option> <seed> <toysperjob> <mu> <lambda> <lifetime>" << std::endl;
    return 0;
  }
  
  // Collect the inputs:
  jobName = argv[1];
  option = argv[2];
  seed = atoi(argv[3]);
  nToysPerJob = atoi(argv[4]);
  inputMuDM = atoi(argv[5]);
  signalProcess = argv[6];
  
  // Copy the input file locally:
  TString inputLocation = Form("%s/%s/workspaces/rootfiles/workspace_%s.root",
			       masterOutput.Data(), jobName.Data(),
			       signalProcess.Data());
  TString nameWS = "workspace_NPP.root";
  system(Form("cp %s %s",inputLocation.Data(),nameWS.Data()));
  
  // Create output TTree:
  outputDir = Form("%s/%s/DMToyMC_%s/", masterOutput.Data(), jobName.Data(),
		   signalProcess.Data());
  TString tempOutputFileName = Form("%s/single_files/toy_mu%i_%i.root",
				  outputDir.Data(),inputMuDM,seed);
  
  // Output file:
  TFile fOutputFile(tempOutputFileName,"recreate");
  TTree fOutputTree("toy","toy");
  
  // Variables to store in the TTree:
  int statusMu1, statusMu0, statusMuFree;
  double muDMVal, nllMu0, nllMu1, nllMuFree, llrL1L0, llrL1Lfree, llrL0Lfree;
  double nevt;
  
  nameNP.clear();
  valueNPMu1.clear();
  valueNPMu0.clear();
  valueNPMuFree.clear();
  nToyEvents.clear();
  nameGlobs.clear();
  valueGlobsMu1.clear();
  valueGlobsMu0.clear();
  valueGlobsMuFree.clear();
  
  fOutputTree.Branch("statusMu1",&statusMu1,"statusMu1/I");
  fOutputTree.Branch("statusMu0",&statusMu0,"statusMu0/I");
  fOutputTree.Branch("statusMuFree",&statusMuFree,"statusMuFree/I");
  fOutputTree.Branch("muDMVal", &muDMVal, "muDMVal/D");
  fOutputTree.Branch("nllMu0", &nllMu0, "nllMu0/D");
  fOutputTree.Branch("nllMu1", &nllMu1, "nllMu1/D");
  fOutputTree.Branch("nllMuFree", &nllMuFree, "nllMuFree/D");
  fOutputTree.Branch("llrL1L0", &llrL1L0, "llrL1L0/D");
  fOutputTree.Branch("llrL1Lfree", &llrL1Lfree, "llrL1Lfree/D");
  fOutputTree.Branch("llrL0Lfree", &llrL0Lfree, "llrL0Lfree/D");
  
  fOutputTree.Branch("seed", &seed, "seed/I");
  fOutputTree.Branch("nevt", &nevt, "nevt/D");
  
  fOutputTree.Branch("nameNP", &nameNP);
  fOutputTree.Branch("valueNPMu1", &valueNPMu1);
  fOutputTree.Branch("valueNPMu0", &valueNPMu0);
  fOutputTree.Branch("valueNPMuFree", &valueNPMuFree);
  fOutputTree.Branch("nToyEvents", &nToyEvents);
  
  fOutputTree.Branch("nameGlobs", &nameGlobs);
  fOutputTree.Branch("valueGlobsMu1", &valueGlobsMu1);
  fOutputTree.Branch("valueGlobsMu0", &valueGlobsMu0);
  fOutputTree.Branch("valueGlobsMuFree", &valueGlobsMuFree);
  
  // initialize canvas:
  can = new TCanvas("can","can",800,800);
  
  // Load model, data, etc. from workspace:
  TFile f(Form("%s",nameWS.Data()), "read");
  RooWorkspace *w = (RooWorkspace*)f.Get("combined");    
  ModelConfig* mc = (ModelConfig*)w->obj("ModelConfig");
  TString dname = "obsData";
  if (option.Contains("FixBkgParam")) dname = "AsimovMu0";
  RooAbsData *obsData = w->data(dname);
  
  // Get the POI (mu):
  RooRealVar* firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();
  statistics::setDefaultPrintLevel(0);
  statistics::setDefaultMinimizer("Minuit");
  
  // Load snapshot of profiling or profile data. Then save as "Profile".
  TString snapshotName = "";
  if (option.Contains("FromSnapshot")) snapshotName = "paramsOrigin";
  else ProfileToData(mc, obsData);
  if (snapshotName != "") cout << w->loadSnapshot(snapshotName) << endl;
  RooArgSet* NuisAndpoi=new RooArgSet();
  NuisAndpoi->add(*(RooArgSet*)mc->GetNuisanceParameters());
  NuisAndpoi->add(*(RooArgSet*)mc->GetParametersOfInterest());
  w->saveSnapshot("Profile",*NuisAndpoi);
  cout << "Snapshot saved." << endl;
  
  double fixedMu = (double)inputMuDM;
  
  //----------------------------------------//
  // loop to generate toys:
  cout << "Generating " << nToysPerJob << " toys with mu = " << fixedMu << endl;
  
  for( int itoy = 0; itoy < nToysPerJob; itoy++ )
  {
    // load the snapshot chosen before the toy loop.
    w->loadSnapshot("Profile");
    if( option.Contains("FixMu") )
    {
      firstPOI->setVal(fixedMu);
      firstPOI->setConstant(true);
    }
    bool m_binned = option.Contains("binned") ? true : false;
    RooAbsData *toyData = generateToyData( w, mc, seed, m_binned );
    nevt = toyData->sumEntries();
    cout << "Printing the toy dataset with " << nevt << " events: " << endl;
    toyData->Print("v");
    
    firstPOI->setConstant(false);
    
    //----------------------------------------//
    // Function to fit toy data:
    map<string, double> results = calcLHCTestStat(mc, toyData, statusMu1, 
						  statusMu0, statusMuFree);
    muDMVal = results["profiled_mu"];
    nllMu0 = results["nllMu0"];
    nllMu1 = results["nllMu1"];
    nllMuFree = results["nllMuFree"];
    
    llrL1L0 = nllMu1 - nllMu0;
    llrL1Lfree = muDMVal > 1.0 ? 0.0 : nllMu1 - nllMuFree;
    llrL0Lfree = muDMVal < 0.0 ? 0.0 : nllMu0 - nllMuFree;
    
    // Fill the tree:
    fOutputTree.Fill();
    
    // Count the toys:
    seed++;
    SafeDelete(toyData);
    fOutputTree.AutoSave("SaveSelf");
    isFirstToy = false;
  }
  f.Close();
  fOutputFile.cd();
  fOutputTree.Write();
  fOutputFile.Close();
  system(Form("rm %s",nameWS.Data()));
  return 0;
}
