////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMWorkspace.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard, Hongtao Yang, Haichen Wang                          //
//  Email: ahard@cern.ch                                                      //
//  Date: 02/04/2015                                                          //
//                                                                            //
//  This class builds the workspace for the dark matter analysis fits.        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMWorkspace.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;
using namespace DMAnalysis;

/**
   Instantiate the class.
   @param newJobName - The name of the job
   @param newDMSignal - The Dark Matter signal to incorporate in the model.
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile"), etc.
   @returns void
*/
DMWorkspace::DMWorkspace(TString newJobName, TString newDMSignal,
			 TString newCateScheme, TString newOptions) {
  jobName = newJobName;
  DMSignal = newDMSignal;
  cateScheme = newCateScheme;
  options = newOptions;
  allGoodFits = true;
  
  combinedWS = NULL;
  mConfig = NULL;
  
  std::cout << "\nDMWorkspace: Initializing..."
	    << "\n\tjobName = " << jobName
	    << "\n\tsignal = " << DMSignal
	    << "\n\tcateScheme = " << cateScheme 
	    << "\n\toptions = " << options << std::endl;
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMWorkspace",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/Plots/",outputDir.Data()));
  system(Form("mkdir -vp %s/rootfiles/",outputDir.Data()));
  system(Form("mkdir -vp %s/mu/",outputDir.Data()));
  
  // Set style for plots:
  CommonFunc::SetAtlasStyle();
    
  muNominalSM = 1;
  dataToPlot = (DMAnalysis::doBlind) ? "asimovDataMu1" : "obsData";

  // Make new or load old workspace:
  if (options.Contains("FromFile")) loadWSFromFile();
  else createNewWS();
  
  if (fitsAllConverged()) {
    std::cout << "DMWorkspace: Successfully initialized!" << std::endl;
  }
  else {
    std::cout << "DMWorkspace: Fit failure during initialization." << std::endl;
  }
  return;
}

/**
   Checks whether all of the fits converged.
   @returns - true iff the fits all converged.
*/
bool DMWorkspace::fitsAllConverged() {
  return allGoodFits;
}

/**
   Retrieves the workspace created by this program.
*/
RooWorkspace* DMWorkspace::getCombinedWorkspace() {
  return combinedWS;
}

/**
   Retrieves a pointer to the model config.
*/
ModelConfig* DMWorkspace::getModelConfig() {
  return mConfig;
}

/**
   Load a previously created workspace.
*/
void DMWorkspace::loadWSFromFile() {
  //Check to see if the workspace has actually been made.
  TFile inputFile(Form("%s/rootfiles/workspaceDM_%s.root", outputDir.Data(),
		       DMSignal.Data()),"read");
  if (inputFile.IsOpen()) {
    std::cout << "Loading workspace from file"<< std::endl;
    combinedWS = (RooWorkspace*)inputFile.Get("combinedWS");
    mConfig = (ModelConfig*)combinedWS->obj("modelConfig");
  }
  else {
    std::cout << "WARNING! Cannot locate requested workspace!"<< std::endl;
    createNewWS();
  }
}

/**
   Create a workspace from scratch. 
*/
void DMWorkspace::createNewWS() {
  
  std::cout << "DMWorkspace: Create a new workspace from scratch." << std::endl;
  std::cout << "\n........................................" << std::endl;
  std::cout << "Workspace parameters:" << std::endl;
  
  // Define and name analysis categories:
  nCategories = DMAnalysis::getNumCategories(cateScheme);
  std::cout << "  Number of categories = " << nCategories << std::endl;
  
  vector<TString> cateNames; cateNames.clear();
  vector<string> cateNamesS; cateNamesS.clear();
  for (int i_c = 0; i_c < nCategories; i_c++) {
    currCateName = Form("%s_%d",cateScheme.Data(),i_c);
    cateNames.push_back(currCateName);
    cateNamesS.push_back((string)currCateName);
    std::cout << "  \t" << currCateName << std::endl;
  }
  std::cout << "Luminosity at 13 TeV: " << analysisLuminosity 
	    << " pb-1." << std::endl;
  std::cout << "........................................" << std::endl;
  
  // Read tables of PES and PER and store values:
  pes = new PESReader(fileNamePESValues, nCategories);
  per = new PERReader(fileNamePERValues, nCategories);
  
  // Instantiate the signal parameterization class using the observable:
  currSigParam = new DMSigParam(jobName, cateScheme, "FromFile", NULL);
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // Everything for simultaneous fit:
  RooWorkspace* cateWS[nCategories];
  RooCategory* categories = new RooCategory(Form("categories_%s",
						 cateScheme.Data()),
					    Form("categories_%s",
						 cateScheme.Data()));
  combinedWS = new RooWorkspace("combinedWS");
  combinedWS->importClassCode();
  
  // Define the combined PDF:
  RooSimultaneous *combinedPdf = new RooSimultaneous("combinedPdf",
						     "combinedPdf",
						     *categories);
  
  // Parameter sets:
  RooArgSet* nuisanceParameters = new RooArgSet();
  RooArgSet* muSMConstants = new RooArgSet();
  RooArgSet* globalObservables = new RooArgSet();
  RooArgSet* observables = new RooArgSet();
  RooArgSet* constraints = new RooArgSet();
  
  // maps for datasets:
  map<string,RooDataSet*> dm;
  map<string,RooDataSet*> dmBinned;
  map<string,RooDataSet*> dmAsimovMu0;
  map<string,RooDataSet*> dmAsimovMu1;
  
  //--------------------------------------//
  // Loop over channels:
  std::cout << "DMWorkspace: Loop over categories to define WS" << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    
    currCateIndex = i_c;
    currCateName = cateNames[i_c];

    // Create the workspace for a single category:
    cateWS[i_c] = createNewCategoryWS();
    categories->defineType(cateNames[i_c]);
    
    // Add category PDF to combined PDF:
    TString namePdf = Form("model_%s",cateNames[i_c].Data());
    combinedPdf->addPdf(*cateWS[i_c]->pdf(namePdf), cateNames[i_c]);
    
    // Add parameter sets to combined workspace:
    TString nameNP = Form("nuisanceParameters_%s",cateNames[i_c].Data());
    TString nameGlob = Form("globalObservables_%s",cateNames[i_c].Data());
    TString nameMuC = Form("muConstants_%s",cateNames[i_c].Data());
    TString nameObs = Form("observables_%s",cateNames[i_c].Data());
    nuisanceParameters->add(*cateWS[i_c]->set(nameNP));
    globalObservables->add(*cateWS[i_c]->set(nameGlob));
    muSMConstants->add(*cateWS[i_c]->set(nameMuC));
    nuisanceParameters->add(*cateWS[i_c]->set(nameMuC));//Added to NP as well.
    observables->add(*cateWS[i_c]->set(nameObs));
    
    // Add category datasets to combined workspace and combined datasets:
    TString nameOD = Form("obsData_%s",cateNames[i_c].Data());
    TString nameODB = Form("obsDataBinned_%s",cateNames[i_c].Data());
    TString nameAD0 = Form("asimovDataMu0_%s",cateNames[i_c].Data());
    TString nameAD1 = Form("asimovDataMu1_%s",cateNames[i_c].Data());
    combinedWS->import(*(RooDataSet*)cateWS[i_c]->data(nameOD));
    combinedWS->import(*(RooDataSet*)cateWS[i_c]->data(nameODB));
    combinedWS->import(*(RooDataSet*)cateWS[i_c]->data(nameAD0));
    combinedWS->import(*(RooDataSet*)cateWS[i_c]->data(nameAD1));
    dm[cateNamesS[i_c]] = (RooDataSet*)combinedWS->data(nameOD);
    dmBinned[cateNamesS[i_c]] =(RooDataSet*)combinedWS->data(nameODB);
    dmAsimovMu0[cateNamesS[i_c]] = (RooDataSet*)combinedWS->data(nameAD0);
    dmAsimovMu1[cateNamesS[i_c]] = (RooDataSet*)combinedWS->data(nameAD1);
  }
  std::cout << "DMWorkspace: Beginning to combine all categories." << std::endl;
  
  // Define the combined datasets:
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*observables);
  args->add(wt);
  RooDataSet* obsData = new RooDataSet("obsData","obsData", *args,
				       Index(*categories), Import(dm), 
				       WeightVar(wt));
  RooDataSet* obsDataBinned = new RooDataSet("obsDataBinned", "obsDataBinned",
					     *args, Index(*categories),
					     Import(dmBinned), WeightVar(wt));
  RooDataSet* asimovDataMu0 = new RooDataSet("asimovDataMu0", "asimovDataMu0",
					     *args, Index(*categories), 
					     Import(dmAsimovMu0),WeightVar(wt));
  RooDataSet* asimovDataMu1 = new RooDataSet("asimovDataMu1", "asimovDataMu1",
					     *args, Index(*categories), 
					     Import(dmAsimovMu1),WeightVar(wt));
  
  // Import PDFs, parameters, and dataset into workspace:
  combinedWS->import(*combinedPdf);
  combinedWS->defineSet("nuisanceParameters", *nuisanceParameters);
  combinedWS->defineSet("observables", *observables);
  combinedWS->defineSet("globalObservables", *globalObservables);
  combinedWS->defineSet("poi", RooArgSet(*combinedWS->var("mu_DM")));   
  combinedWS->defineSet("muSMConstants", *muSMConstants);
  combinedWS->import(*obsData);
  combinedWS->import(*obsDataBinned);
  combinedWS->import(*asimovDataMu0);
  combinedWS->import(*asimovDataMu1);
  
  // Define the ModelConfig:
  mConfig = new ModelConfig("modelConfig",combinedWS);
  mConfig->SetPdf((*combinedWS->pdf("combinedPdf")));
  mConfig->SetObservables((*combinedWS->set("observables")));
  mConfig->SetParametersOfInterest((*combinedWS->set("poi")));
  mConfig->SetNuisanceParameters((*combinedWS->set("nuisanceParameters")));
  mConfig->SetGlobalObservables((*combinedWS->set("globalObservables")));
  combinedWS->import(*mConfig);
  
  std::cout << "DMWorkspace: Printing the combined workspace." << std::endl;
  combinedWS->Print("v");
  
  // Start profiling the data:
  std::cout << "DMWorkspace: Start profiling data" << std::endl;
  
  // Profile and save snapshots of the data:
  double nllMu0, nllMu1, nllMuFree, profiledMuValue;
  profileAndSnapshot("0", nllMu0, profiledMuValue);
  profileAndSnapshot("1", nllMu1, profiledMuValue);
  profileAndSnapshot("Free", nllMuFree, profiledMuValue);
  
  // Print summary of the fits:
  std::cout.precision(10);
  std::cout << "\nPrinting likelihood results: " << std::endl;
  std::cout << "\tnll(muDM = 1):  " << nllMu1 << std::endl;
  std::cout << "\tnll(muDM = 0):  " << nllMu0 << std::endl;
  std::cout << "\tnll(muDM free): " << nllMuFree << std::endl;
  std::cout << " " << endl;
  std::cout << "\tnll(S+B)/nll(B) " << nllMu1 - nllMu0 << std::endl;
  std::cout << "\tnll(muDM=1)/nll(muhat) = " << nllMu1 - nllMuFree << std::endl;
  std::cout << "\tnll(muDM=0)/nll(muhat) = " << nllMu0 - nllMuFree << std::endl;
  if (allGoodFits) std::cout << "allGoodFits = TRUE" << std::endl;
  else std::cout << "allGoodFits = FALSE" << std::endl;
  std::cout << "Profiled muDM value : " << profiledMuValue << std::endl;
  
  // Write the profiled mu value to file:
  ofstream fileMuProf;
  fileMuProf.open(Form("%s/mu/mu_%s.txt", outputDir.Data(), DMSignal.Data()));
  fileMuProf << profiledMuValue << endl;
  fileMuProf.close();
  
  // Write workspace to file:
  combinedWS->writeToFile(Form("%s/rootfiles/workspaceDM_%s.root",
			       outputDir.Data(), DMSignal.Data()));
}

/**
   Create the workspace for a single analysis category.
   @param currCategory
*/
RooWorkspace* DMWorkspace::createNewCategoryWS() {
  
  // The bools that control the systematic uncertainties:
  bool inclusive = currCateName == "inclusive";
  bool channel_constraints_attached = (currCateIndex == 0);
  bool m_norm = !options.Contains("nonorm");
  bool m_pes = !options.Contains("nopes");
  bool m_per = !options.Contains("noper");
  bool m_ss  = !options.Contains("noss");
  bool m_bgm = !options.Contains("nobgm");
  bool m_mig = !options.Contains("nomig");
  bool m_nosys = options.Contains("nosys");
  if (m_nosys) {
    std::cout << "\tDMWorkspace: ALL systematics OFF" << endl;
    m_norm = false;   m_pes = false;   m_per = false;
    m_ss = false;     m_bgm = false;   m_mig = false;
  }
  std::cout << "\tNormalization systematics = " << m_norm << std::endl;
  std::cout << "\tEnergy scale systematics  = " << m_pes  << std::endl;
  std::cout << "\tResolution systematics    = " << m_per  << std::endl;
  std::cout << "\tShape systematics         = " << m_ss   << std::endl;
  std::cout << "\tBackground systematics    = " << m_bgm  << std::endl;
  std::cout << "\tMigration systematics     = " << m_mig  << std::endl;
  
  //--------------------------------------//
  // Create the individual channel workspace:
  RooWorkspace *tempWS = new RooWorkspace(Form("tmpWS_%s",currCateName.Data()));
  
  // nuispara:
  RooArgSet *nuisParams = new RooArgSet();
  RooArgSet *nuisParamsBkg = new RooArgSet();
  RooArgSet *nuisParamsUncorrelated = new RooArgSet();
  // constraints:
  RooArgSet *constraints = new RooArgSet();
  RooArgSet *constraintsBias = new RooArgSet();
  // globobs:
  RooArgSet *globalObs = new RooArgSet();
  RooArgSet *globalObsProc = new RooArgSet();
  // expected:
  RooArgSet *expectedShape = new RooArgSet();
  RooArgSet *expectedBias = new RooArgSet();
  RooArgSet *expected = new RooArgSet();
  //RooArgSet *expectedSM = new RooArgSet();
  RooArgSet *expectedDM = new RooArgSet();
  /*
  RooArgSet *expectedProc_ggH = new RooArgSet();
  RooArgSet *expectedProc_VBF = new RooArgSet();
  RooArgSet *expectedProc_WH = new RooArgSet();
  RooArgSet *expectedProc_ZH = new RooArgSet();
  RooArgSet *expectedProc_bbH = new RooArgSet();
  RooArgSet *expectedProc_ttH = new RooArgSet();
  */
  
  // array setup[5] is used to configure a nuisance parameter
  // [0]    [1]       [2]   [3]     
  // sigma, sigmalow, beta, nominal,
  
  //--------------------------------------//
  // Normalization systematics:
  if (m_norm) {
    double setupLumi[4] = {0.036, 0, 1, 1};
    makeNP("Luminosity", setupLumi, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupTrigger[4] = {0.005, 0, 1, 1};
    makeNP("Trigger", setupTrigger, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupIsEM[4] = {0.0526, 0, 1, 1};
    makeNP("PhotonID", setupIsEM, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupIso[4] = {0.004, 0, 1, 1};
    makeNP("Isolation", setupIso, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
    double setupESCALE[4] = {0.003, 0, 1, 1};
    makeNP("ESCALE", setupESCALE, *&nuisParams, *&constraints, *&globalObs,
	   *&expected);
  }
  
  //--------------------------------------//
  // Migration systematics:
  /*
  if (m_mig) {
    int number_SS_sources = ss_tool->GetNumberOfSources(energy);
    // loop over ss sources.
    for( int i_s = 0; i_s < number_SS_sources; i_s++ )
    {
      TString current_SS_source_name = ss_tool->GetNameOfSource( i_s, energy );
      TString ss_np_name = Form("shape_%s",current_SS_source_name.Data());
      double current_ss_value = ss_tool->GetValue( current_SS_source_name, currCateIndex, energy );
      int current_ss_sign = ss_tool->GetSign( current_SS_source_name, currCateIndex, energy );
      
      // Asymmetric migration uncertainties:
      double setup_ss_current[4] = {current_ss_value, 0, current_ss_sign, 1};
      if( current_SS_source_name.Contains("_up") )
      {
	TString current_SS_source_name_down = ss_tool->GetNameOfSource( i_s+1, energy );
	double current_ss_value_down = ss_tool->GetValue( current_SS_source_name_down, currCateIndex, energy );
	setup_ss_current[1] = current_ss_value_down;
	ss_np_name.ReplaceAll("_up","");
      }
      // down values must follow the up case in the list, and are included in the step above
      if( current_SS_source_name.Contains("_down") ) continue;
      
      // spin0 shape systematics:
      makeNP(ss_np_name, setup_ss_current, nuispara, constraints, globobs, expected_spin0p);
      
    }
  }
  */
  //--------------------------------------//
  // SYSTEMATICS: Spurious signal
  if (m_bgm) {
    double ssEvents = spuriousSignal();
    double setupBias[4] = {ssEvents, -999, 1, 0}; //Gaussian constraint
    makeNP("bias", setupBias, *&nuisParamsUncorrelated, *&constraintsBias,
	   *&globalObs, *&expectedBias);
    // two lines below were added.
    RooProduct sigBias("sigBias","sigBias",*expectedBias);
    tempWS->import(sigBias);
  }
  else tempWS->factory("sigBias[0]");//expectedBias
  
  //--------------------------------------//
  // SYSTEMATICS: Resolution:
  vector<TString> perList; perList.clear();
  if (m_per) {
    double setupPER[4] = {0.0, 0, 1, 1};
    // Loop over sources of resolution systematic uncertainty:
    for (int i_s = 0; i_s < per->getNumberOfSources(); i_s++) {
      TString currPERSource = per->getNameOfSource(i_s);
      TString currPERName = Form("EM_%s",currPERSource.Data());
      perList.push_back(currPERName);
      setupPER[0] = per->getValue(currPERSource, currCateIndex);
      setupPER[2] = per->getSign(currPERSource, currCateIndex);
      
      // resolution on the inclusive shape:
      makeShapeNP(currPERName, "DM", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      //makeShapeNP(currPERName, "SM", setupPER, *&nuisParams, *&constraints,
      //	  *&globalObs, *&expectedShape);
      /*
      makeShapeNP(currPERName, "ggH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "VBF", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "WH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "ZH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "bbH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
      makeShapeNP(currPERName, "ttH", setupPER, *&nuisParams, *&constraints,
		  *&globalObs, *&expectedShape);
		  */   
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Energy-scale
  vector<TString> pesList; pesList.clear();
  if (m_pes) {
    double setupPES[4] = {0.0, 0, 1, 1};
    // loop over sources of energy scale systematic uncertainty:
    for (int i_s = 0; i_s < pes->getNumberOfSources(); i_s++) {
      TString currPESSource = pes->getNameOfSource(i_s);
      TString currPESName = Form("EM_%s",currPESSource.Data());
      pesList.push_back(currPESName);
      setupPES[0] = pes->getValue(currPESSource, currCateIndex);
      setupPES[2] = pes->getSign(currPESSource, currCateIndex);
      makeNP(currPESName, setupPES, *&nuisParams, *&constraints, *&globalObs,
	     *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // Parameters of interest (POIs):
  double muMin = 0; double muMax = 100;
  RooRealVar *mu_DM = new RooRealVar("mu_DM", "mu_DM", 1, muMin, muMax);
  //RooRealVar *mu_SM = new RooRealVar("mu_SM", "mu_SM", 1, muMin, muMax);
  /*
  RooRealVar *mu_ggH = new RooRealVar("mu_ggH", "mu_ggH", 1, muMin, muMax);
  RooRealVar *mu_VBF = new RooRealVar("mu_VBF", "mu_VBF", 1, muMin, muMax);
  RooRealVar *mu_WH = new RooRealVar("mu_WH", "mu_WH", 1, muMin, muMax);
  RooRealVar *mu_ZH = new RooRealVar("mu_ZH", "mu_ZH", 1, muMin, muMax);
  RooRealVar *mu_bbH = new RooRealVar("mu_bbH", "mu_bbH", 1, muMin, muMax);
  RooRealVar *mu_ttH = new RooRealVar("mu_ttH", "mu_ttH", 1, muMin, muMin);
  */
  expectedDM->add(RooArgSet(*mu_DM));
  //expectedSM->add(RooArgSet(*mu_SM));
  /*
  expectedProc_ggH->add(RooArgSet(*mu_ggH));
  expectedProc_VBF->add(RooArgSet(*mu_VBF));
  expectedProc_WH->add(RooArgSet(*mu_WH));
  expectedProc_ZH->add(RooArgSet(*mu_ZH));
  expectedProc_bbH->add(RooArgSet(*mu_bbH));
  expectedProc_ttH->add(RooArgSet(*mu_ttH));
  */
  
  // Expectation values:
  RooProduct expectationDM("expectationDM","expectationDM", *expectedDM);
  //RooProduct expectationSM("expectationSM","expectationSM", *expectedSM);
  //RooProduct expectationCommon("expectationCommon","expectationCommon",
  //			       *expected);
  /*
  RooProduct expectationProc_ggH("expectationProc_ggH","expectationProc_ggH",
				 *expectedProc_ggH);
  RooProduct expectationProc_VBF("expectationProc_VBF","expectationProc_VBF",
				 *expectedProc_VBF);
  RooProduct expectationProc_WH("expectationProc_WH","expectationProc_WH",
				*expectedProc_WH);
  RooProduct expectationProc_ZH("expectationProc_ZH","expectationProc_ZH",
				*expectedProc_ZH);
  RooProduct expectationProc_bbH("expectationProc_bbH","expectationProc_bbH",
				 *expectedProc_bbH);
  RooProduct expectationProc_ttH("expectationProc_ttH","expectationProc_ttH",
				 *expectedProc_ttH);
				 */
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  //tempWS->import(expectationCommon);
  //tempWS->import(expectationSM);
  tempWS->import(expectationDM);
  /*
  tempWS->import(expectationProc_ggH);
  tempWS->import(expectationProc_VBF);
  tempWS->import(expectationProc_WH);
  tempWS->import(expectationProc_ZH);
  tempWS->import(expectationProc_bbH);
  tempWS->import(expectationProc_ttH);
  */
  tempWS->import(*expectedShape);
  tempWS->import(*expectedBias);
  
  // Declare the observable m_yy, and the observables set:
  tempWS->factory(Form("m_yy[%f,%f]",DMMyyRangeLo,DMMyyRangeHi));
  tempWS->defineSet("obsprelim","m_yy");
  
  // Construct the signal PDFs:
  std::cout << "DMWorkspace: Adding signal parameterizations." << std::endl;
  currSigParam->addSigToCateWS(tempWS, pesList, perList, DMSignal,
			       currCateIndex);
  //currSigParam->addSigToCateWS(tempWS, pesList, perList, "SM", currCateIndex);
  /*
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "ggH", currCateIndex);
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "VBF", currCateIndex);
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "WH", currCateIndex);
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "ZH", currCateIndex);
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "bbH", currCateIndex);
  currSigParam->addSigToCateWS(tempWS, pesList, perList, "ttH", currCateIndex);
  */
  // Construct the background PDF:
  DMBkgModel *currBkgModel = new DMBkgModel(jobName, cateScheme, "FromFile", 
  					    tempWS->var("m_yy"));
  currBkgModel->addBkgToCateWS(tempWS, nuisParamsBkg, currCateIndex);
  
  // Add background parameters to uncorrelated collection:
  nuisParamsUncorrelated->add(*nuisParamsBkg);
    
  // Normalization for each process follows such pattern:
  // Definition of expectationCommon = mu*isEM*lumi*migr
  //tempWS->factory(Form("prod::nSigSM(nSM[%f],expectationCommon,expectationSM)",
  //		       currSigParam->getCateSigYield(currCateIndex,"SM")));
  
  //tempWS->factory(Form("prod::nSigDM(nDM[%f],expectationCommon,expectationDM)",
  //		       currSigParam->getCateSigYield(currCateIndex,DMSignal)));
  
  tempWS->factory(Form("prod::nSigDM(nDM[%f],expectationDM)",
		       currSigParam->getCateSigYield(currCateIndex,DMSignal)));
  
  /*
  tempWS->factory(Form("prod::nSigggH(nggH[%f],expectationCommon,expectationProc_ggH)", currSigParam->getCateSigYield(currCateIndex,"ggH")));
  tempWS->factory(Form("prod::nSigVBF(nVBF[%f],expectationCommon,expectationProc_VBF)", currSigParam->getCateSigYield(currCateIndex,"VBF")));
  tempWS->factory(Form("prod::nSigWH(nWH[%f],expectationCommon,expectationProc_WH)", currSigParam->getCateSigYield(currCateIndex,"WH")));
  tempWS->factory(Form("prod::nSigZH(nZH[%f],expectationCommon,expectationProc_ZH)", currSigParam->getCateSigYield(currCateIndex,"ZH")));
  tempWS->factory(Form("prod::nSigbbH(nbbH[%f],expectationCommon,expectationProc_bbH)", currSigParam->getCateSigYield(currCateIndex,"bbH")));
  tempWS->factory(Form("prod::nSigttH(nttH[%f],expectationCommon,expectationProc_ttH)", currSigParam->getCateSigYield(currCateIndex,"ttH")));
  */
  
  // Model with combined SM production modes:
  //tempWS->factory("SUM::modelSB(nSigSM*sigPdfSM,nSigDM*sigPdfDM,sigBias*sigPdfDM,nBkg*bkgPdf)");
  tempWS->factory("SUM::model(nSigDM*sigPdfDM,nBkg*bkgPdf)");
  // Model with separated SM production modes:
  //tempWS->factory("SUM::modelProdSB(nSigggH*sigPdfggH,nSigVBF*sigPdfVBF,nSigWH*sigPdfWH,nSigZH*sigPdfZH,nSigbbH*sigPdfbbH,nSigttH*sigPdfttH,nSigDM*sigPdfDM,expectedBias*sigPdfDM,nBkg*bkgPdf)");
  
  // Only attach constraint term to first category. If constraint terms were
  // attached to each category, constraints would effectively be multiplied.
  if (currCateIndex == 0) {
    constraints->add(*constraintsBias);
    RooProdPdf constraint("constraint", "constraint", *constraints);
    tempWS->import(constraint);
    //tempWS->factory("PROD::model(modelSB,constraint)");
    //tempWS->factory("PROD::modelProd(modelProdSB,constraint)");
  }
  // Except in the case where the constraints are uncorrelated between
  // categories, as with the spurious signal:
  else {
    RooProdPdf constraint("constraint", "constraint", *constraintsBias);
    tempWS->import(constraint);
    //tempWS->factory("PROD::model(modelSB,constraint)");
    //tempWS->factory("PROD::modelProd(modelProdSB,constraint)");
  }
  
  /*
    Specify the group of nuisance parameters that are correlated between
    categories. Technically, this is done by sharing the same name for nuisance
    parameter between sub-channels. Their respective global observables should
    also share the same name. nuisParams should contain all correlated nuisance
    parameters. All uncorrelated nuisance parameters should be included in
    nuisParamsUncorrelated.
  */
  //TString corrNPNames = "mu_DM,mu_SM,mu_ggH,mu_VBF,mu_WH,mu_ZH,mu_bbH,mu_ttH";
  //TString corrNPNames = "mu_DM,mu_SM";
  TString corrNPNames = "mu_DM";
  
  // Iterate over nuisance parameters:
  TIterator *iterNuis = nuisParams->createIterator();
  RooRealVar* currNuis;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    std::cout << "\t" << currNuis->GetName() << std::endl;
    corrNPNames += Form(",nuisPar_%s,globOb_%s",currNuis->GetName(),
			currNuis->GetName());
  }
  std::cout << "For category " << currCateName << ", correlate variables: "
	    << corrNPNames << std::endl;
  
  /*
    Sub-channel labeling
    Import the workspace tempWS to another workspace and add currCateName as a 
    suffix to all nodes and variables of w. the correlated nuisance parameters
    and their respective global observables will not be renamed.
  */
  RooWorkspace* categoryWS = new RooWorkspace("workspace_"+currCateName);
  categoryWS->import((*tempWS->pdf("model")), RenameAllNodes(currCateName),
		     RenameAllVariablesExcept(currCateName,corrNPNames),
		     Silence());
  /*
  categoryWS->import((*tempWS->pdf("modelProd")), RenameAllNodes(currCateName),
		     RenameAllVariablesExcept(currCateName,corrNPNames),
		     Silence());
  */
  // Adding correlated nuisance parameters to nuisanceParameters:
  RooArgSet* nuisCateWS = new RooArgSet();
  iterNuis->Reset();
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(currNuis->GetName()));
  }
  
  // Adding uncorrelated nuisance parameters to nuisanceParameters:
  TIterator *iterNuisUncorrelated = nuisParamsUncorrelated->createIterator();
  RooRealVar* currNuisUncorrelated;
  while ((currNuisUncorrelated = (RooRealVar*)iterNuisUncorrelated->Next())) {
    TString nuisName = (currNuisUncorrelated->GetName() 
			+ (TString)"_" + currCateName);
    nuisCateWS->add(*(RooRealVar*)categoryWS->obj(nuisName));
  }
  
  // Adding unconstrained NPs from the background pdf:
  RooArgSet* nuisBkgCateWS = new RooArgSet();
  TIterator *iterNuisBkg = nuisParamsBkg->createIterator();
  RooRealVar* currNuisBkg;
  while ((currNuisBkg = (RooRealVar*)iterNuisBkg->Next())) {
    TString parName = currNuisBkg->GetName()+(TString)"_"+currCateName;
    nuisBkgCateWS->add(*categoryWS->var(parName));
  }
  
  /*
    Global observables:
    Global observables only appear in the constraint terms. All constraint terms
    of correlated nuisance parameters are attached to the pdf of the first
    subchannel. For those global observables, their names should be the same as
    those in the w. For other subchannels, only the bias constraint term is
    attached.
  */  
  RooArgSet *globsCateWS = new RooArgSet();
  TIterator *iterGlobs = globalObs->createIterator();
  RooRealVar *currGlobs;
  while ((currGlobs = (RooRealVar*)iterGlobs->Next())) {
    TString globName = currGlobs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(globName)) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(globName));
      categoryWS->var(globName)->setConstant();
    }
    else if ((bool)categoryWS->obj(currGlobs->GetName())) {
      globsCateWS->add(*(RooRealVar*)categoryWS->obj(currGlobs->GetName()));
      categoryWS->var(currGlobs->GetName())->setConstant();
    }
  }
  
  RooArgSet *obsCateWS = new RooArgSet();
  TIterator *iterObs = tempWS->set("obsprelim")->createIterator();
  RooRealVar *currObs;
  while ((currObs = (RooRealVar*)iterObs->Next())) {
    TString obsName = currObs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(obsName)) {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(obsName));
    }
    else {
      obsCateWS->add(*(RooRealVar*)categoryWS->obj(currObs->GetName()));
    }
  }
  
  // Set some of the mu values constant:
  std::cout << "DMWorkspace: Setting SM signals constant." << std::endl;
  RooArgSet* muConstCateWS = new RooArgSet();
  //muConstCateWS->add(*categoryWS->var("mu_SM"));
  /*
  muConstCateWS->add(*categoryWS->var("mu_ggH"));
  muConstCateWS->add(*categoryWS->var("mu_VBF"));
  muConstCateWS->add(*categoryWS->var("mu_WH"));
  muConstCateWS->add(*categoryWS->var("mu_ZH"));
  muConstCateWS->add(*categoryWS->var("mu_bbH"));
  muConstCateWS->add(*categoryWS->var("mu_ttH"));
  */ 
  TIterator *iterMuConst = muConstCateWS->createIterator();
  RooRealVar *currMuConst;
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    currMuConst->setVal(1.0);
    currMuConst->setConstant(true);
  }
  
  categoryWS->defineSet(Form("muConstants_%s",currCateName.Data()),
			*muConstCateWS);
  categoryWS->defineSet(Form("observables_%s",currCateName.Data()),
			*obsCateWS);
  categoryWS->defineSet(Form("nuisanceParameters_%s",currCateName.Data()),
			*nuisCateWS);
  categoryWS->defineSet(Form("globalObservables_%s",currCateName.Data()),
			*globsCateWS);
  
  // Import the observed data set:
  DMMassPoints *currMassPoints = NULL;
  if (DMAnalysis::doBlind) {
    currMassPoints = new DMMassPoints(jobName,"gg_gjet",cateScheme,"FromFile",
				      categoryWS->var("m_yy_"+currCateName));
  }
  else {
    currMassPoints = new DMMassPoints(jobName, "data", cateScheme, "FromFile",
				      categoryWS->var("m_yy_"+currCateName));
  }
  
  RooDataSet *obsData = currMassPoints->getCateDataSet(currCateIndex);
  TString obsDataName = Form("obsData_%s",currCateName.Data());
  obsData->SetNameTitle(obsDataName, obsDataName);
  
  // Set the background normalization parameter:
  std::cout << "DMWorkspace: Fitting bkg. in " << currCateName << std::endl;
  (*categoryWS->var("nBkg_"+currCateName)).setVal(obsData->numEntries());
  (*categoryWS->pdf("bkgPdf_"+currCateName)).fitTo(*obsData, Minos(RooArgSet(*nuisBkgCateWS)), SumW2Error(kTRUE));//should be false?
  (*categoryWS->var("nBkg_"+currCateName)).setVal(obsData->numEntries());
  
  categoryWS->import(*obsData);
    
  // Create Asimov mu_DM = 0,1 data:
  RooRealVar wt("wt","wt",1);
  createAsimovData(categoryWS, obsData, wt, 0, muNominalSM);
  createAsimovData(categoryWS, obsData, wt, 1, muNominalSM);
  
  // Create a binned observed data set:
  RooArgSet* obsPlusWt = new RooArgSet();
  obsPlusWt->add(wt);
  obsPlusWt->add(*categoryWS->var("m_yy_"+currCateName));
  
  // Create a histogram to store binned data:
  TH1F* h_data = new TH1F("h_data", "h_data", 240, DMMyyRangeLo, DMMyyRangeHi);
  RooArgSet* obsArgSet = (RooArgSet*)obsData->get();
  RooRealVar* massVar = (RooRealVar*)obsArgSet->find("m_yy_"+currCateName);
  for (int i_e = 0; i_e < obsData->numEntries(); i_e++) {
    obsData->get(i_e);
    h_data->Fill(massVar->getVal());
  }
  
  // Fill obsdatabinned dataset with binned data:
  TString obsDataBinName = Form("obsDataBinned_%s", currCateName.Data());
  RooDataSet *obsDataBinned = new RooDataSet(obsDataBinName, obsDataBinName,
					     *obsPlusWt, WeightVar(wt));
  int nBin = h_data->GetNbinsX();
  for (int i_b = 1; i_b < nBin; i_b++) {
    // 240 bins -> 0.25 GeV per bin
    double massVal = h_data->GetBinCenter(i_b);
    double weightVal = h_data->GetBinContent(i_b);
    categoryWS->var("m_yy_"+currCateName)->setVal(massVal);
    wt.setVal(weightVal);
    obsDataBinned->add(RooArgSet(*categoryWS->var("m_yy_"+currCateName),wt),
		       weightVal);
  }
  categoryWS->import(*obsDataBinned);
    
  // Plot the single-channel fit:
  plotSingleCateFit(categoryWS, Form("obsData_%s",currCateName.Data()));
  plotSingleCateFit(categoryWS, Form("asimovDataMu1_%s",currCateName.Data()));
  
  std::cout << "DMWorkspace: Printing workspace for category:" << currCateName
	    << std::endl;
  categoryWS->Print("v");
  
  delete h_data;
  return categoryWS;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// spuriousSignal:

double DMWorkspace::spuriousSignal() {
  double spurious[10] = { 1.0, 1.0, 1.0, 1.0, 1.0 };// per fb-1
  return spurious[currCateIndex] * analysisLuminosity;
}

/**
   A private method for constructing a normalization systematic uncertainty.
   @param varName - the name of the systematic uncertainty.
   @param setup - the systematic uncertainty configuration parameters.
   @param nuisParams - the set of nuisance parameters to which this will add.
   @param constraints - the set of constraints to which this will add a term. 
   @param globalObs - the set of global observables, to which this will add.
   @param expected - the set of expected terms. 
   @returns - void. 
*/
void DMWorkspace::makeNP(TString varName, double setup[4],
			 RooArgSet *&nuisParams, RooArgSet *&constraints,
			 RooArgSet *&globalObs, RooArgSet *&expected) {
  
  // Settings for the systematic:
  double sigma    = setup[0];
  double sigmaLow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  
  RooWorkspace* workspace = new RooWorkspace(varName);
  std::cout << "Creating a nuisance parameter " << varName << std::endl;

  // Create nuisance parameter with asymmetric uncertainty:
  if (sigmaLow > 0) {
    std::cout << "  parameter has an asymmetric uncertainty" << std::endl;
    
    RooRealVar* var = new RooRealVar(Form("nuisPar_%s",varName.Data()),
				     Form("nuisPar_%s",varName.Data()),
				     0,-5,5);
    RooRealVar* varBeta = new RooRealVar(Form("beta_%s",varName.Data()), 
					  Form("beta_%s",varName.Data()),
					  beta);
    RooProduct* varXBeta = new RooProduct(Form("%s_times_beta",varName.Data()),
					  Form("%s_times_beta",varName.Data()),
					  RooArgSet(*var,*varBeta));
    vector<double> sVarHi; sVarHi.clear(); sVarHi.push_back(1+sigma);
    vector<double> sVarLo; sVarLo.clear(); sVarLo.push_back(1-sigmaLow);
    RooArgList nuisList(*varXBeta);
    
    TString expName = Form("expected_%s",varName.Data());
    RooStats::HistFactory::FlexibleInterpVar expVar(expName, expName, nuisList,
						    nominal, sVarLo, sVarHi);
    workspace->import(expVar);
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Create nuisance parameter with gaussian (symmetric) uncertainty:
  else if (sigmaLow == -999) {
    std::cout << "  parameter has a Gaussian constraint term" << std::endl;
    
    workspace->factory(Form("sum::expected_%s(nominal_%s[%f], prod::uncer_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), sigma_%s[%f]))", varName.Data(), varName.Data(), nominal, varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data(), sigma));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Create nuisance parameter with bifurcated Gaussian constraint term:
  else if (sigmaLow<0 && sigmaLow != -999) {
    
    TString valLogKappa = Form("%f",sqrt(log(1+pow(sigma,2))));
    TString valA = Form("%f",fabs(sigma/sigmaLow)); 
    
    std::cout << "  parameter has a Bif. Gauss constraint term" << std::endl;
    std::cout << "  asymmetric factor is " << valA << std::endl;
    
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), valLogKappa_%s)", varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s, nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooBifurGauss::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1,%s)", varName.Data(), varName.Data(), varName.Data(), valA.Data()));
  }
  
  // Create a nuisance parameter with log-normal constraint term:
  else {
    std::cout << "  parameter has logNormal constraint term" << std::endl;
    TString valLogKappa = Form("%f",sqrt(log(1+pow(sigma,2))));
    
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), valLogKappa_%s)", varName.Data(), varName.Data(), varName.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s,nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varName.Data(), varName.Data(), varName.Data()));
  }
  
  // Add parameters and PDFs to relevant sets:
  nuisParams->add(*workspace->var(Form("nuisPar_%s",varName.Data())));
  constraints->add(*workspace->pdf(Form("constrPdf_%s",varName.Data())));
  globalObs->add(*workspace->var(Form("globOb_%s",varName.Data())));
  expected->add(*workspace->function(Form("expected_%s",varName.Data())));
}

/**
   A private method for constructing a shape systematic uncertainty. The shape
   NP maker will give variables used in parameterization process dependent name,
   but keep the same name for nuisance parameter, and global observables.
   @param varNameNP - the name of the systematic uncertainty.
   @param process - the corresponding signal process for the systematic.
   @param setup - the systematic uncertainty configuration parameters.
   @param nuisParams - the set of nuisance parameters to which this will add.
   @param constraints - the set of constraints to which this will add a term. 
   @param globalObs - the set of global observables, to which this will add.
   @param expected - the set of expected terms. 
   @returns - void. 
*/
void DMWorkspace::makeShapeNP(TString varNameNP, TString process,
			      double setup[4], RooArgSet *&nuisParams,
			      RooArgSet *&constraints, RooArgSet *&globalObs,
			      RooArgSet *&expected) {
  
  // Settings for the systematic:
  double sigma    = setup[0];
  double sigmaLow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  TString varName = varNameNP + process;
  
  RooWorkspace* workspace = new RooWorkspace(varName);
  
  // Create nuisance parameter with asymmetric uncertainty:
  if (sigmaLow > 0) {
    std::cout << "  parameter for an asymmetric uncertainty" << std::endl;
    
    RooRealVar* var = new RooRealVar(Form("nuisPar_%s",varNameNP.Data()),
				     Form("nuisPar_%s",varNameNP.Data()),
				     0, -5, 5);
    RooRealVar* varBeta = new RooRealVar(Form("beta_%s",varName.Data()),
					  Form("beta_%s",varName.Data()),
					  beta);
    RooProduct* varXBeta = new RooProduct(Form("%s_times_beta",varName.Data()),
					  Form("%s_times_beta",varName.Data()),
					  RooArgSet(*var,*varBeta));
    
    vector<double> sVarHi; sVarHi.clear(); sVarHi.push_back(1+sigma);
    vector<double> sVarLo; sVarLo.clear(); sVarLo.push_back(1-sigmaLow);
    RooArgList nuisList(*varXBeta);
    RooStats::HistFactory::FlexibleInterpVar expVar(Form("expected_%s",varName.Data()), Form("expected_%s",varName.Data()), nuisList, nominal, sVarLo, sVarHi);
    
    workspace->import(expVar);
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Create a nuisance parameter with Gaussian constraint term:
  else if (sigmaLow == -999) {
    std::cout << "  parameter with a Gaussian constraint term" << std::endl;
    
    workspace->factory(Form("sum::expected_%s(nominal_%s[%f], prod::uncer_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]), sigma_%s[%f]))", varName.Data(), varName.Data(), nominal, varName.Data(), varName.Data(), varNameNP.Data(), varName.Data(), beta, varName.Data(), sigma));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Create a nuisance parameter with log-normal constraint term:
  else {
    TString valLogKappa = Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    workspace->factory(Form("valLogKappa_%s[%s]", varName.Data(), valLogKappa.Data()));
    workspace->factory(Form("RooExponential::expTerm_%s(prod::%s_times_beta(nuisPar_%s[0,-5,5], beta_%s[%f]),valLogKappa_%s)", varName.Data(), varName.Data(), varNameNP.Data(), varName.Data(), beta, varName.Data()));
    workspace->factory(Form("prod::expected_%s(expTerm_%s,nominal_%s[%f])", varName.Data(), varName.Data(), varName.Data(), nominal));
    workspace->factory(Form("RooGaussian::constrPdf_%s(globOb_%s[0,-5,5],nuisPar_%s,1)", varNameNP.Data(), varNameNP.Data(), varNameNP.Data()));
  }
  
  // Add parameters and PDFs to relevant sets:
  nuisParams->add(*workspace->var(Form("nuisPar_%s",varNameNP.Data())));
  constraints->add(*workspace->pdf(Form("constrPdf_%s",varNameNP.Data())));
  globalObs->add(*workspace->var(Form("globOb_%s",varNameNP.Data())));
  expected->add(*workspace->function(Form("expected_%s",varName.Data())));
}

/**
   Create Asimov data for the statistical model, using a fit to observed data
   for the shape and normalizaiton of the background.
   @param cateWS - the current category workspace.
   @param obsData - the observed dataset for the category.
   @param wt - the weight variable for the dataset.
   @param valMuDM - the value of dark matter signal strength to use.
   @returns - void.
*/
void DMWorkspace::createAsimovData(RooWorkspace* cateWS, RooDataSet *obsData,
				   RooRealVar wt, int valMuDM, int valMuSM) {
  
  std::cout << "DMWorkspace: Creating Asimov data, mu=" << valMuDM << std::endl;
  
  int nPointsAsimov = 275;
  // This is the dataset to be returned:
  TString asimovName = Form("asimovDataMu%d_%s", valMuDM, currCateName.Data());
  RooDataSet *asimovData = new RooDataSet(asimovName, asimovName, RooArgSet(*cateWS->var("m_yy_"+currCateName),wt), WeightVar(wt));
  // Load the PDF from the workspace:
  RooAbsPdf *currPdf = (RooAbsPdf*)(cateWS->pdf("model_"+currCateName));
  double initialMuDM = (cateWS->var("mu_DM"))->getVal();
  //double initialMuSM = (cateWS->var("mu_SM"))->getVal();
  cateWS->var("mu_DM")->setVal((double)valMuDM);
  //cateWS->var("mu_SM")->setVal((double)valMuSM);
  cateWS->var("mu_DM")->setConstant(true);
  //cateWS->var("mu_SM")->setConstant(true);
  
  // Use fit result to get the estimate of the background:
  double totalEvents = (obsData->sumEntries() +
			currSigParam->getCateSigYield(currCateIndex,DMSignal));
  //currSigParam->getCateSigYield(currCateIndex,"SM") + 
  

  double width = (DMMyyRangeHi - DMMyyRangeLo) / ((double)nPointsAsimov);
  cateWS->var("m_yy_"+currCateName)->setRange("nominalRange", DMMyyRangeLo,
					      DMMyyRangeHi);
  
  RooAbsReal *totIntegral = (RooAbsReal*)currPdf->createIntegral(RooArgSet(*cateWS->var("m_yy_"+currCateName)), NormSet(*cateWS->var("m_yy_"+currCateName)), Range("nominalRange"));

  // Loop over the number of Asimov points:
  double countAsimov = 0.0;
  for ( int i_p = 0; i_p < nPointsAsimov; i_p++) {
  
    double massVal = DMMyyRangeLo + (0.5*width) + (width*(double)i_p);
    cateWS->var("m_yy_"+currCateName)->setRange(Form("rangeIntegral%d",i_p),
						massVal-(0.5*width),
						massVal+(0.5*width));
    
    RooAbsReal *currIntegral = (RooAbsReal*)currPdf->createIntegral(RooArgSet(*cateWS->var("m_yy_"+currCateName)), NormSet(*cateWS->var("m_yy_"+currCateName)), Range(Form("rangeIntegral%d",i_p)));
   
    double weightVal = totalEvents * (currIntegral->getVal() / 
				      totIntegral->getVal());
    countAsimov += weightVal;
    cateWS->var("m_yy_"+currCateName)->setVal(massVal);
    wt.setVal(weightVal);
    asimovData->add(RooArgSet(*cateWS->var("m_yy_"+currCateName),wt),weightVal);
  }
  /*
  if (fabs((countAsimov - obsData->sumEntries()) / countAsimov) > 0.04) {
    std::cout << "Bad Asimov Data: Data = " << obsData->sumEntries()
	      << ", Asimov = " << countAsimov << std::endl;
    exit(0);
  }
  */
  cateWS->var("m_yy_"+currCateName)->setRange("nominalRange", DMMyyRangeLo,
					      DMMyyRangeHi);
  cateWS->import(*asimovData);
  cateWS->var("mu_DM")->setVal(initialMuDM);
  //cateWS->var("mu_SM")->setVal(initialMuSM);
  cateWS->var("mu_DM")->setConstant(false);
  //cateWS->var("mu_SM")->setConstant(true);
  std::cout << "DMWorkspace: Asimov data has " << asimovData->sumEntries() 
	    << " entries" << std::endl;
}

/**
   Plot the fits produced by the specified model.
   @param plotOptions - options for what fits to plot etc.
   @returns void
*/
void DMWorkspace::plotSingleCateFit(RooWorkspace *cateWS, TString dataset) {
  std::cout << "DMWorkspace: Plot single category fit for "
	    << currCateName << std::endl;
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  RooPlot* frame =  (*cateWS->var("m_yy_"+currCateName)).frame(55);
  cateWS->data(dataset)->plotOn(frame);
  (*cateWS->pdf("model_"+currCateName)).plotOn(frame, LineColor(2));
  (*cateWS->pdf("model_"+currCateName)).plotOn(frame, Components((*cateWS->pdf("bkgPdf_"+currCateName))), LineColor(4));
  (*cateWS->pdf("model_"+currCateName)).plotOn(frame, Components((*cateWS->pdf("sigPdfDM_"+currCateName))), LineColor(3));
  //(*cateWS->pdf("model_"+currCateName)).plotOn(frame, Components((*cateWS->pdf("sigPdfSM_"+currCateName))), LineColor(5));
  
  //double chi2 = frame->chiSquare();
  frame->SetYTitle("Events / GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  frame->Draw();
  
  TLatex text; text.SetNDC(); text.SetTextColor(1);
  text.DrawLatex(0.2, 0.81, Form("Category %d", currCateIndex));
  text.DrawLatex(0.2, 0.87, Form("Signal %s", DMSignal.Data()));
  //TH1F *histSM = new TH1F("histSM", "histSM", 1, 0, 1);
  TH1F *histDM = new TH1F("histDM", "histDM", 1, 0, 1);
  TH1F *histNR = new TH1F("histNR", "histNR", 1, 0, 1);
  TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
  histDM->SetLineColor(3);
  //histSM->SetLineColor(5);
  histNR->SetLineColor(4);
  histSig->SetLineColor(2);
  TLegend leg(0.61, 0.6, 0.89, 0.77);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.SetBorderSize(0);
  leg.AddEntry(histSig, "Sig. + bkg.", "l");
  leg.AddEntry(histDM, "Dark matter", "l");
  //leg.AddEntry(histSM, "SM Higgs", "l");
  leg.AddEntry(histNR, "Non-resonant", "l");
  leg.Draw("SAME");
  can->Print(Form("%s/Plots/cateFit_%s_%s_%s.eps", outputDir.Data(),
		  dataset.Data(), DMSignal.Data(), currCateName.Data()));
  delete can;
}

/**
   Plot the fits produced by the specified model.
   @param plotOptions - options for what fits to plot etc.
   @returns void
*/
void DMWorkspace::plotFinalFits(RooWorkspace *combWS, TString fitType) {
  std::cout << "DMWorkspace: Plot final fits for " << fitType << std::endl;
  
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  
  // loop over categories:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    currCateName = Form("%s_%d", cateScheme.Data(), i_c);
    currCateIndex = i_c;
    RooPlot* frame =  (*combWS->var("m_yy_"+currCateName)).frame(55);
    //dataMap[(string)currCateName]->plotOn(frame);
    combWS->data(Form("%s_%s", dataToPlot.Data(), currCateName.Data()))->plotOn(frame);
    (*combWS->pdf("model_"+currCateName)).plotOn(frame, Components((*combWS->pdf("sigPdfDM_"+currCateName))), LineColor(3));
    //(*combWS->pdf("model_"+currCateName)).plotOn(frame, Components((*combWS->pdf("sigPdfSM_"+currCateName))), LineColor(5));
    (*combWS->pdf("model_"+currCateName)).plotOn(frame, Components((*combWS->pdf("bkgPdf_"+currCateName))), LineColor(4));
    (*combWS->pdf("model_"+currCateName)).plotOn(frame, LineColor(2));
    
    //double chi2 = frame->chiSquare();
    frame->SetYTitle("Events / GeV");
    frame->SetXTitle("M_{#gamma#gamma} [GeV]");
    frame->Draw();
    
    TLatex text; text.SetNDC(); text.SetTextColor(1);
    text.DrawLatex(0.2, 0.81, Form("Category %d", currCateIndex));
    text.DrawLatex(0.2, 0.87, Form("Signal %s", DMSignal.Data()));
    TH1F *histDM = new TH1F("histDM", "histDM", 1, 0, 1);
    //TH1F *histSM = new TH1F("histSM", "histSM", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDM->SetLineColor(3);
    //histSM->SetLineColor(5);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    TLegend leg(0.61, 0.63, 0.89, 0.77);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.SetBorderSize(0);
    leg.AddEntry(histDM, "Dark matter", "l");
    //leg.AddEntry(histSM, "SM Higgs", "l");
    leg.AddEntry(histBkg, "Non-resonant", "l");
    leg.AddEntry(histSig, "Sig. + bkg.", "l");
    leg.Draw("SAME");
    can->Print(Form("%s/Plots/combFit_%s_%s_%s.eps", outputDir.Data(),
		    DMSignal.Data(), fitType.Data(), currCateName.Data()));
  }
  delete can;
}

/**
   Plot the pulls of the nuisance parameters related to systematic uncertainties
   in each of the fits.
   @param nuisParams - the set of nuisance parameters
   @param type - the type of fit
*/
void DMWorkspace::plotNuisParams(RooArgSet nuisParams, TString type) {
  // TBD
}

// "1", "0", "Free"
void DMWorkspace::profileAndSnapshot(TString muDMValue, double &nllValue,
				     double &profiledMu) {
  std::cout << "\nDMWorkspace: Profile mu_DM = " << muDMValue << std::endl;
  
  RooRealVar *poi = (RooRealVar*)mConfig->GetParametersOfInterest()->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*mConfig->GetNuisanceParameters());
  poiAndNuis->add(*poi);
  combinedWS->saveSnapshot("paramsOrigin",*poiAndNuis);
  RooArgSet* globs = (RooArgSet*)mConfig->GetGlobalObservables();
  RooArgSet* nuis = (RooArgSet*)mConfig->GetNuisanceParameters();
  
  RooAbsPdf *pdf = mConfig->GetPdf();
  //RooAbsPdf *pdf = (RooAbsPdf*)combinedPdf;
  //RooAbsPdf *pdf = (RooAbsPdf*)combinedWS->pdf("combinedPdf");
    
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  
  // Choose the constant mu_DM value:
  if (muDMValue.EqualTo("1")) poi->setVal(1.0);
  else poi->setVal(0.0);
  
  // Choose free or fixed:
  if (muDMValue.EqualTo("Free")) poi->setConstant(false);
  else poi->setConstant(true);
  
  // Set the mu_SM fixed and to the nominal value.
  //combinedWS->var("mu_SM")->setVal(muNominalSM);
  //combinedWS->var("mu_SM")->setConstant(true);
  
  // Perform the fit:
  RooFitResult* resMu = pdf->fitTo(*combinedWS->data(dataToPlot), 
				   Constrain(*nuis),
				   Extended(pdf->canBeExtended()),
				   PrintLevel(0), Save(true));
  // Track whether all fits converge, save results:
  if (resMu->status() != 0) allGoodFits = false;
  nllValue = resMu->minNll();
  profiledMu = poi->getVal();
  
  combinedWS->saveSnapshot(Form("paramsProfileMu%s",muDMValue.Data()),
			   *poiAndNuis);
  
  // Plots of invariant mass and nuisance parameters:
  if (!options.Contains("noplot")) {
    plotFinalFits(combinedWS, Form("mu%s",muDMValue.Data()));
    if (!options.Contains("nosys")) {
      plotNuisParams(*combinedWS->set("nuisanceParameters"),
		     Form("mu%s",muDMValue.Data()));
    }
  }
  // Then revert to original nuisance parameter values and settings:
  combinedWS->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  //combinedWS->var("mu_SM")->setConstant(true);
}
