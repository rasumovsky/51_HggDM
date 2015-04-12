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
//  First: build signal and background models.                                //
//  Second: add asimov data function.                                         //
//  Third: make plots a la spin analysis or better yet NPP.                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                      //
//  Options:                            //
//    noess - no energy scale sys.      //
//    nores - no resolution sys.        //
//    noss  - no shape sys.             //
//    nosys - no sys. period            //
//    decorrmu                          //
//    fitasimov0p                       //
//    fitasimov2p                       //
//                                      //
//////////////////////////////////////////

#include "DMWorkspace.h"

/**
   Instantiate the class.
   @param newJobName - The name of the job
   @param newCateScheme - The name of the event categorization
   @param newOptions - The job options ("New", "FromFile"), etc.
   @returns void
*/
DMWorkspace::DMWorkspace(TString newJobName, TString newCateScheme,
			 TString newOptions) {
  jobName = newJobName;
  cateScheme = newCateScheme;
  options = newOptions;
  
  // Assign output directory, and make sure it exists:
  outputDir = Form("%s/%s/DMWorkspace",masterOutput.Data(),jobName.Data());
  system(Form("mkdir -vp %s",outputDir.Data()));
  system(Form("mkdir -vp %s/figures/",outputDir.Data()));
  
  // Set style for plots:
  SetAtlasStyle();
  
  // Instantiate important classes for global info.
  selector = new DMEvtSelect();// for cate info
  
  // Make new or load old workspace:
  if (options.Contains("FromFile")) loadWSFromFile();
  else createNewWS();
  return;
}

/**
   Create a workspace from scratch. 
*/
void DMWorkspace::createNewWS() {
  
  std::cout << "DMWorkspace: Create a new workspace from scratch." << std::endl;
  std::cout << "\n........................................" << std::endl;
  std::cout << "Workspace parameters:" << std::endl;
  
  // Define and name analysis categories:
  int nCategories = selector->getNCategories(cateScheme);
  std::cout << "  Number of categories = " << nCategories << std::endl;
  
  vector<TString> cateNames; cateNames.clear();
  vector<string> cateNamesS; cateNamesS.clear();
  for (int i_c = 0; i_c < nCategories; i_c++) {
    currCateName = Form("%s_%d",newCateScheme.Data(),i_c);
    cateNames.push_back(currCateName);
    cateNamesS.push_back((string)currCateName);
    std::cout << "  \t" << currCateName << std::endl;
  }
  std::cout << "Luminosity at 13 TeV:" << analysisLuminosity << std::endl;
  std::cout << "........................................" << std::endl;
  
  // Read tables of ESS and Res and store values:
  //ess_tool = new ESSReader( file_name_ESS_values, nCategories);
  //res_tool = new ResReader( file_name_Res_values, nCategories);
  //ss_tool  = new SigShapeReader( file_name_SS_values, nCategories);
  
  // Instantiate the signal parameterization class using the observable:
  currSigParam = new DMSigParam(jobName, cateScheme, "FromFile");
  // Instantiate the background parameterization class using the observable:
  currBkgModel = new DMBkgModel(jobName, cateScheme, "FromFile");
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // Everything for simultaneous fit:
  RooWorkspace* cateWS[nCategories];
  RooCategory* categories = new RooCategory(Form("categories_%s",
						 newCateScheme.Data()),
					    Form("categories_%s",
						 newCateScheme.Data()));
  
  combinedWS = new RooWorkspace("combinedWS");
  combinedWS->importClassCode();
  RooSimultaneous combinedPdf("combinedPdf","",*categories);
  
  // Parameter sets:
  RooArgSet* nuisanceParameters = new RooArgSet();
  RooArgSet* globalObservables = new RooArgSet();
  RooArgSet* observables = new RooArgSet();
  RooArgSet* constraints = new RooArgSet();

  // maps for datasets:
  map<string,RooDataSet*> dataMap;
  map<string,RooDataSet*> dataMapBinned;
  map<string,RooDataSet*> dataMapAsimov;
  
  //--------------------------------------//
  // Loop over channels:
  std::cout << "DMWorkspace: Loop over categories to define WS" << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    
    currCateIdex = i_c;
    currCateName = cateNames[i_c];

    // Create the workspace for a single category:
    cateWS[i_c] = createNewCategoryWS();
    categories->defineType(cateNames[i_c]);
    
    // Add PDF and parameters from category to global collections:
    combinedPdf.addPdf(*cateWS[i_c]->pdf("model_"+CN[i_c]),CN[i_c]);
    nuisanceParameters->add(*cateWS[i_c]->set("nuisanceParameters"));
    globalObservables->add(*cateWS[i_c]->set("globalObservables"));
    observables->add(*cateWS[i_c]->set("observables"));
    
    // Retrieve the datasets produced for each category:
    dataMap[cateNamesS[i_c]] = (RooDataSet*)w[i_c]->data("obsData");
    dataMapBinned[cateNamesS[i_c]] = (RooDataSet*)w[i_c]->data("obsDataBinned");
    dataMapAsimov[cateNamesS[i_c]] = (RooDataSet*)w[i_c]->data("asimovBinned");
  }
  
  // Import PDFs and parameters to combined workspace:
  combinedWS->import(combinedPdf);
  combinedWS->defineSet("nuisanceParameters",*nuisanceParameters);
  combinedWS->defineSet("observables",*observables);
  combinedWS->defineSet("globalObservables",*globalObservables);
  combinedWS->defineSet("poi",RooArgSet(*combinedWS->var("mu")));   
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*observables);
  args->add(wt);
  RooDataSet* obsData = new RooDataSet("obsData","combined data",*args,
				       Index(*categories), Import(dataMap), 
				       WeightVar(wt));
  RooDataSet* obsDataBinned = new RooDataSet("obsDatabinned",
					     "combined data binned",*args,
					     Index(*categories),
					     Import(dataMapBinned),
					     WeightVar(wt));
  RooDataSet* asimovData = new RooDataSet("asimovData","combined asimov data",
					  *args, Index(*categories), 
					  Import(dataMapAsimov), WeightVar(wt));
 
  combinedWS->import(*obsData);
  combinedWS->import(*obsDataBinned);
  combinedWS->import(*asimovData);
  
  // Define the ModelConfig:
  mconfig = new ModelConfig("mconfig",combinedWS);
  mconfig->SetPdf(*combinedWS->pdf("combinedPdf"));
  mconfig->SetObservables(*combinedWS->set("observables"));
  mconfig->SetParametersOfInterest((*combinedWS->set("poi")));
  mconfig->SetNuisanceParameters((*combinedWS->set("nuisanceParameters")));
  mconfig->SetGlobalObservables((*combinedWS->set("globalObservables")));
  combinedWS->import(*mconfig);
  
  // Start profiling the data:
  cout << "Start profiling data" << endl;
  RooRealVar *poi = (RooRealVar*)mconfig->GetParametersOfInterest()->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*mconfig->GetNuisanceParameters());
  RooArgSet* globs = (RooArgSet*)mconfig->GetGlobalObservables();
  poiAndNuis->add(*poi);
  poiAndNuis->Print();
  combinedWS->saveSnapshot("paramsOrigin",*poiAndNuis);
  RooAbsPdf *pdf = mconfig->GetPdf();
  
  can = new TCanvas("can","can",800,800);
  
  //----------------------------------------//
  // Profile and save snapshot for mu=1 fixed:
  std::cout << "\nProfile mu = 1:" << std::endl;
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  poi->setVal(1.0);
  poi->setConstant(true);
  RooFitResult* resMu1;
  if (options.Contains("fitasimov")) {
    resMu1 = pdf->fitTo(*combinedWS->data("asimovMu1"), PrintLevel(0),
			Save(true));
  }
  else {
    resMu1 = pdf->fitTo(*combinedWS->data("obsData"), PrintLevel(0),
			Save(true));
  }
  
  // Track whether all fits converge:
  if (resMu1->status() != 0) AllGoodFits = false;
  
  double nllMu1 = resMu1->minNll();
  combinedWS->saveSnapshot("paramsProfileMu1",*poiAndNuis);
  
  // Plots of invariant mass and nuisance parameters:
  if (!option.Contains("noplot")) {
    if (options.Contains("fitasimov")) {
      PlotFinalFits(combinedWS,"asimovMu1","mu1");
    }
    else {
      PlotFinalFits(combinedWS,"obsData","mu1");
    }
    if (!options.Contains("nosys")) {
      PlotNuisParams(*combinedWS->set("ModelConfig_NuisParams"), "mu1");
    }
  }
  
  //----------------------------------------//
  // Profile and save snapshot for mu=0 fixed:
  std::cout << "\nProfile mu = 0:" << std::endl;
  combinedWS->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  poi->setVal(0.0);
  poi->setConstant(true);
  RooFitResult* resMu0;
  if (options.Contains("fitasimov")) {
    resMu0 = pdf->fitTo(*combinedWS->data("asimovMu1"), PrintLevel(0),
			Save(true));
  }
  else {
    resMu0 = pdf->fitTo(*combinedWS->data("obsData"), PrintLevel(0),
			Save(true));
  }
  
  // Track whether all fits converge:
  if (resMu0->status() != 0) AllGoodFits = false;
  
  double nllMu0 = resMu0->minNll();
  combinedWS->saveSnapshot("paramsProfileMu0",*poiAndNuis);

  // Plots of invariant mass and nuisance parameters:
  if (!options.Contains("noplot")) {
    if (options.Contains("fitasimov")) {
      PlotFinalFits(combinedWS, "asimovMu1", "mu0");
    }
    else {
      PlotFinalFits(combinedWS, "obsData", "mu0");
    }
    if (!options.Contains("nosys")) {
      PlotNuisParams(*combinedWS->set("ModelConfig_NuisParams"), "mu0");
    }
  }
  
  //----------------------------------------//
  // Profile and save snapshot for mu floating:
  std::cout << "\nProfile mu floating:" << std::endl;
  combinedWS->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  poi->setVal(0.0);
  poi->setConstant(false);
  RooFitResult* resMuFree;
  if (options.Contains("fitasimov")) {
    resMuFree = pdf->fitTo(*combinedWS->data("asimovMu1"), PrintLevel(0),
			    Save(true));
  }
  else {
    resMuFree = pdf->fitTo(*combinedWS->data("obsData"), PrintLevel(0), 
			    Save(true));
  }
  
  // Track whether all fits converge:
  if (resMsuFree->status() != 0) AllGoodFits = false;
  
  double nllMuFree = resMuFree->minNll();
  double profiledMuValue = poi->getVal();
  combinedWS->saveSnapshot("paramsProfile_muFree",*poiAndNuis);
 
  // Plots of invariant mass and nuisance parameters:
  if (!option.Contains("noplot")) {
    if (option.Contains("fitasimov")) {
      PlotFinalFits(combinedWS, "asimovMu1", "mufree");
    }
    else {
      PlotFinalFits(combinedWS, "obsData", "mufree");
    }
    if (!option.Contains("nosys")) {
      PlotNuisParams(*combinedWS->set("ModelConfig_NuisParams"), "mufree");
    }
  }
  
  combinedWS->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  
  // Print summary of the fits:
  std::cout.precision(10);
  std::cout << "\nPrinting likelihood results: " << std::endl;
  std::cout << "\tnll(mu = 1):  " << nllMu1 << std::endl;
  std::cout << "\tnll(mu = 0):  " << nllMu0 << std::endl;
  std::cout << "\tnll(mu free): " << nllMuFree << std::endl;
  std::cout << " " << endl;
  std::cout << "\tnll(S+B)/nll(B) " << nllMu1 - nllMu0 << std::endl;
  std::cout << "\tnll(mu=1)/nll(muhat) = " << nllMu1 - nllMuFree << std::endl;
  std::cout << "\tnll(mu=0)/nll(muhat) = " << nllMu0 - nllMuFree << std::endl;
  if (AllGoodFits) std::cout << "AllGoodFits = TRUE" << std::endl;
  else std::cout << "AllGoodFits = FALSE" << std::endl;
  std::cout << "Profiled mu value : " << profiledMuValue << std::endl;
  
  //ofstream file_muprof;
  //file_muprof.open(Form("%s/profiled_mu_%iTeV_%ips.txt",output_directory_ws.Data(),lambda,lifetime));
  //file_muprof << profiled_mu_value << endl;
  //file_muprof.close();
  
  // Write workspace to file:
  combinedWS->writeToFile(Form("%s/workspaceDM.root",outputDir.Data()));
}

/**
   Create the workspace for a single analysis category.
   @param currCategory
*/
RooWorkspace* DMWorkspace::createNewCategoryWS() {
  
  // The bools that control the systematic uncertainties:
  bool inclusive = currCateName == "inclusive";
  bool channel_constraints_attached = (currCateName == leadingchannel);
  bool m_norm = !options.Contains("nonorm");
  bool m_ess = !options.Contains("noess");
  bool m_res = !options.Contains("nores");
  bool m_ss  = !options.Contains("noss");
  bool m_bgm = !options.Contains("nobgm");
  bool m_mig = !options.Contains("nomig");
  bool m_nosys = options.Contains("nosys");
  if (m_nosys) {
    std::cout << "\tALL systematics = OFF" << endl;
    m_norm = false;
    m_ess = false;
    m_res = false;
    m_ss  = false;
    m_bgm = false;
    m_mig = false;
  }
  std::cout << "\tNormalization systematics = " << m_norm << std::endl;
  std::cout << "\tEnergy scale systematics  = " << m_ess << std::endl;
  std::cout << "\tResolution systematics    = " << m_res << std::endl;
  std::cout << "\tShape systematics         = " << m_ss<< std::endl;
  std::cout << "\tBackground systematics    = " << m_bgm << std::endl;
  std::cout << "\tMigration systematics     = " << m_mig << std::endl;
  
  //--------------------------------------//
  // Create the individual channel workspace:
  currWS = new RooWorkspace(currCateName);
  
  // nuispara:
  RooArgSet *nuisParams = new RooArgSet();
  RooArgSet *nuisParamsNull = new RooArgSet();
  RooArgSet *nuisParamsBkg = new RooArgSet();
  RooArgSet *nuisParamsUncorrelated = new RooArgSet();
  RooArgSet *nuisParamsProc = new RooArgSet();
  // constraints:
  RooArgSet *constraints = new RooArgSet();
  RooArgSet *constraintsNull = new RooArgSet();
  RooArgSet *constraintsBias = new RooArgSet();
  RooArgSet *constraintsProc = new RooArgSet();
  // globobs:
  RooArgSet *globalObs = new RooArgSet();
  RooArgSet *globalObsNull = new RooArgSet();
  RooArgSet *globalObsProc = new RooArgSet();
  // expected:
  RooArgSet *expected = new RooArgSet();
  RooArgSet *expectedSM = new RooArgSet();
  RooArgSet *expectedDM = new RooArgSet();
  RooArgSet *expectedShape = new RooArgSet();
  RooArgSet *expectedBias = new RooArgSet();
  RooArgSet *expectedProc_ggF = new RooArgSet();
  RooArgSet *expectedProc_VBF = new RooArgSet();
  RooArgSet *expectedProc_WH = new RooArgSet();
  RooArgSet *expectedProc_ZH = new RooArgSet();
  RooArgSet *expectedProc_ttH = new RooArgSet();
    
  // array setup[5] is used to configure a nuisance parameter
  // [0]    [1]       [2]   [3]      [4]
  // sigma, sigmalow, beta, nominal, nonATLAS
  
  //--------------------------------------//
  // Normalization systematics:
  if (m_norm) {
    double setupLumi[5] = {0.036, 0, 1, 1, 1};
    NPmaker("Luminosity", setupLumi, *&nuisParams, *&constraints, *&globalObs,
	    *&expected);
    double setupTrigger[5] = {0.005, 0, 1, 1, 1};
    NPmaker("Trigger", setupTrigger, *&nuisParams, *&constraints, *&globalObs,
	    *&expected);
    double setupIsEM[5] = {0.0526, 0, 1, 1, 1};
    NPmaker("PhotonID", setupIsEM, *&nuisParams, *&constraints, *&globalObs,
	    *&expected);
    double setupIso[5] = {0.004, 0, 1, 1, 1};
    NPmaker("Isolation", setupIso, *&nuisParams, *&constraints, *&globalObs,
	    *&expected);
    double setupESCALE[5] = {0.003, 0, 1, 1, 1};
    NPmaker("ESCALE", setupESCALE, *&nuisParams, *&constraints, *&globalObs,
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
      double setup_ss_current[5] = {current_ss_value, 0, current_ss_sign, 1, 1};
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
      NPmaker(ss_np_name, setup_ss_current, nuispara, constraints, globobs, expected_spin0p);
      
    }
  }
  */
  //--------------------------------------//
  // SYSTEMATICS: Spurious signal
  if (m_bgm) {
    double ss_events = spurious_signal(currCateName);
    double setup_bias[5] = {ss_events, -999, 1, 0, 1}; //Gaussian constraint
    NPmaker("bias", setup_bias, *&nuisParamsUncorrelated, *&constraintsBias,
	    *&globalObs, *&expectedBias);
  }
  else currWS->factory("expectedBias[0]");
  
  //--------------------------------------//
  // SYSTEMATICS: Resolution:
  vector<TString> resList; resList.clear();
  if (m_res) {
    double setup_AllRes[5] = {0.0, 0, 1, 1, 1};
    
    // Loop over sources of resolution systematic uncertainty:
    for (int i_s = 0; i_s < res_tool->GetNumberOfSources(); i_s++) {
      TString currResSource = res_tool->GetNameOfSource(i_s);
      TString currResName = Form("EM_%s",currResSource.Data());
      resList.push_back(currResName);
      setupRes[0] = res_tool->GetValue(currResSource, currCateIndex);
      setupRes[2] = res_tool->GetSign(currResSource, currCateIndex);
      
      // resolution on the inclusive shape:
      shapeNPmaker(currResName, "_inc", setupRes, *&nuisParams, *&constraints,
		   *&globalObs, *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Energy-scale
  vector<TString> essList; essList.clear();
  if( !m_noess )
  {
    double setup_AllESS[5] = {0.0, 0, 1, 1, 1};
    
    // loop over sources of energy scale systematic uncertainty:
    for (int i_s = 0; i_s < ess_tool->GetNumberOfSources(); i_s++) {
      TString currESSSource = ess_tool->GetNameOfSource(i_s);
      TString currESSName = Form("EM_%s",currESSSource.Data());
      essList.push_back(currESSName);
      setupESS[0] = ess_tool->GetValue(currESSSource, currCateIndex);
      setupESS[2] = ess_tool->GetSign(currESSSource, currCateIndex);
      NPmaker(currESSName, setupESS, *&nuisParams, *&constraints, *&globalObs,
	      *&expectedShape);
    }
  }
  
  //--------------------------------------//
  // Parameters of interest (POIs):
  //RooRealVar *mu_DM = new RooRealVar("mu_DM","mu_DM",1,-100,100);
  //RooRealVar *mu_SM = new RooRealVar("mu_SM","mu_SM",1,-100,100);
  RooRealVar *mu_BR_gg = new RooRealVar("mu_BR_gg","mu_BR_gg",1,-100,100);
  RooRealVar *mu_ggF = new RooRealVar("mu_ggF","mu_ggF",1,-100,100);
  RooRealVar *mu_VBF = new RooRealVar("mu_VBF","mu_VBF",1,-100,100);
  RooRealVar *mu_WH = new RooRealVar("mu_WH","mu_WH",1,-100,100);
  RooRealVar *mu_ZH = new RooRealVar("mu_ZH","mu_ZH",1,-100,100);
  RooRealVar *mu_ttH = new RooRealVar("mu_ttH","mu_ttH",1,-100,100);
  //expected->add(RooArgSet(*mu, *mu_BR_gg));NO
  expectedSM->add(RooArgSet(*mu_SM, *mu_BR_gg));
  expectedDM->add(RooArgSet(*mu_DM, *mu_BR_gg));
  expectedProc_ggF->add(RooArgSet(*mu_ggF));
  expectedProc_VBF->add(RooArgSet(*mu_VBF));
  expectedProc_WH->add(RooArgSet(*mu_WH));
  expectedProc_ZH->add(RooArgSet(*mu_ZH));
  expectedProc_ttH->add(RooArgSet(*mu_ttH));
  
  // Expectation values:
  RooProduct expectationSM("expectationSM","expectationSM", *expectedSM);
  RooProduct expectationDM("expectationDM","expectationDM", *expectedDM);
  RooProduct expectationCommon("expectationCommon","expectationCommon",
			       *expected);
  RooProduct expectationProc_ggF("expectationProc_ggF","expectationProc_ggF",
				 *expectedProc_ggF);
  RooProduct expectationProc_VBF("expectationProc_VBF","expectationProc_VBF",
				 *expectedProc_VBF);
  RooProduct expectationProc_WH("expectationProc_WH","expectationProc_WH",
				*expectedProc_WH);
  RooProduct expectationProc_ZH("expectationProc_ZH","expectationProc_ZH",
				*expectedProc_ZH);
  RooProduct expectationProc_ttH("expectationProc_ttH","expectationProc_ttH",
				 *expectedProc_ttH);
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  currWS->import(expectationSM);
  currWS->import(expectationDM);
  currWS->import(expectationCommon);
  currWS->import(expectationProc_ggF);
  currWS->import(expectationProc_VBF);
  currWS->import(expectationProc_WH);
  currWS->import(expectationProc_ZH);
  currWS->import(expectationProc_ttH);
  currWS->import(*expectedShape);
  currWS->import(*expectedBias);
  
  // Declare the observable m_yy, and the observables set:
  currWS->factory(Form("m_yy[%f,%f]",DMMyyRangeLo,DMMyyRangeHi));
  currWS->defineSet("observables","m_yy");
  
  // Construct the signal PDFs:
  currSigParam->addSigToCateWS(currWS, essList, resList, "DM", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "SM", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ggH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "VBF", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "WH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ZH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ttH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "bbH", currCateIndex);
  
  // Construct the background PDF:
  currBkgModel->addBkgToCateWS(currWS, nuisParamsBkg, currCateName);
  
  // Add background parameters to uncorrelated collection:
  nuisParamsUncorrelated->add(*nuisParamsBkg);
  
  // build the signal normalization
  TString nSM = Form("%f", currSigParam->getCateSigYield(currCateIndex,"SM"));
  TString nDM = Form("%f", currSigParam->getCateSigYield(currCateIndex,"DM"));
  std::cout << "\tnSM for " << currCateName << " = " << nSM << std::endl;
  std::cout << "\tnDM for " << currCateName << " = " << nDM << std::endl;
  
  // Normalization for each process follows such pattern:
  // mu*isEM*lumi*migr => expectationCommon
  w->factory((TString)"prod::nSigSM(nSM["+nDM+(TString)"],expectationCommon,expectationSM)");
  w->factory((TString)"prod::nSigDM(nDM["+nSM_2p+(TString)"],expectationCommon,expectationDM)");
  w->factory("SUM::modelSB(nSigSM*sigPdfSM,nSigDM*sigPdfDM,expectedBias*sigPdfInc,nBkg*bkgPdf)");
  w->Print();
  
  // Only attach constraint term to first category. If constraint terms were
  // attached to each category, constraints would effectively be multiplied.
  if (currCateIndex == 0) {
    constraints->add(*constraints_bias);
    RooProdPdf constraint("constraint", "constraint", *constraints);
    w->import(constraint);
    w->factory("PROD::model(modelSB,constraint)");
  }
  // Except in the case where the constraints are uncorrelated between
  // categories, as with the spurious signal:
  else {
    RooProdPdf constraint("constraint","constraint",*constraintsBias);
    w->import(constraint);
    w->factory("PROD::model(modelSB,constraint)");
  }
  
  /*
    Specify the group of nuisance parameters that are correlated between
    categories. Technically, this is done by sharing the same name for nuisance
    parameter between sub-channels. Their respective global observables should
    also share the same name. nuisParams should contain all correlated nuisance
    parameters. All uncorrelated nuisance parameters should be included in
    nuisParamsUncorrelated.
  */
  TString corrNPNames = "mu_DM,mu_SM,mu_BR_gg";
  
  // Iterate over nuisance parameters:
  TIterator *iterNuis = nuisParams->createIterator();
  RooRealVar* currNuis;
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    std::cout << "\t" << currNuis->GetName() << std::endl;
    corrNPNames += ("," + currNuis->GetName() + ",R_" + currNuis->GetName());
  }
  std::cout << "For category " << currCateName
	    << " the following variables will be correlated: "
	    << corrNPNames << std::endl;
  /*
    Sub-channel labeling
    Import the workspace currWS to another workspace and add currCateName as a 
    suffix to all nodes and variables of w. the correlated nuisance parameters
    and their respective global observables will not be renamed.
  */
  RooWorkspace* categoryWS = new RooWorkspace("workspace"+currCateName);
  categoryWS->import( (*currWS->pdf("model")), RenameAllNodes(currCateName),
		      RenameAllVariablesExcept(currCateName,corrNPNames),
		      Silence());
  
  // Adding correlated nuisance parameters to nuisanceParameters:
  RooArgSet* nuisance_categoryWS = new RooArgSet();
  iterNuis->Reset();
  while ((currNuis = (RooRealVar*)iterNuis->Next())) {
    nuisance_categoryWS->add(*(RooRealVar*)categoryWS->obj(currNuis->GetName()));
  }
  
  // Adding uncorrelated nuisance parameters to nuisanceParameters:
  TIterator *iterNuisUncorrelated = nuisParamsUncorrelated->createIterator();
  RooRealVar* currNuisUncorrelated;
  while ((currNuisUncorrelated = (RooRealVar*)iterNuisUncorrelated->Next())) {
    TString nuisName = (currNuisUncorrelated->GetName() 
			+ (TString)"_" + currCateName);
    nuisance_categoryWS->add(*(RooRealVar*)categoryWS->obj(currNuisName));
  }
  
  // Adding unconstrained NPs from the background pdf:
  RooArgSet* nuispara_bkg_categoryWS = new RooArgSet();
  TIterator *iterNuisBkg = nuisParamsBkg->createIterator();
  RooRealVar* currNuisBkg;
  while ((currNuisBkg = (RooRealVar*)iterNuisBkg->Next())) {
    TString parName = currNuisBkg->GetName()+(TString)"_"+currCateName;
    nuispara_bkg_categoryWS->add(*categoryWS->var(parName));
  }
  
  /*
    Global observables:
    Global observables only appear in the constraint terms. All constraint terms
    of correlated nuisance parameters are attached to the pdf of the first
    subchannel. For those global observables, their names should be the same as
    those in the w. For other subchannels, only the bias constraint term is
    attached.
  */  
  RooArgSet *global_categoryWS = new RooArgSet();
  TIterator *iterGlobs = globalObs->createIterator();
  RooRealVar *currGlobs;
  while ((currGlobs = (RooRealVar*)iterGlobs->Next())) {
    
    TString globName = currGlobs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(globName)) {
      global_categoryWS->add(*(RooRealVar*)categoryWS->obj(globName));
      categoryWS->var(globName)->setConstant();
    }
    else if ((bool)categoryWS->obj(currGlobs->GetName())) {
      global_categoryWS->add(*(RooRealVar*)categoryWS->obj(currGlobs->GetName()));
      categoryWS->var(currGlobs->GetName())->setConstant();
    }
  }
  
  RooArgSet *observable_categoryWS = new RooArgSet();
  TIterator *iterObs = currWS->set("observables")->createIterator();
  RooRealVar *currObs;
  while (currObs = (RooRealVar*)iterObs->Next()) {
    TString obsName = currObs->GetName()+(TString)"_"+currCateName;
    if ((bool)categoryWS->obj(obsName)) {
      observable_categoryWS->add(*(RooRealVar*)categoryWS->obj(obsName));
    }
    else {
      observable_categoryWS->add(*(RooRealVar*)categoryWS->obj(currObs->GetName()));
    }
  }
  
  RooArgSet* muConstants_categoryWS = new RooArgSet();
  muConstants_categoryWS->add(*categoryWS->var("mu_ggF"));
  muConstants_categoryWS->add(*categoryWS->var("mu_VBF"));
  muConstants_categoryWS->add(*categoryWS->var("mu_WH"));
  muConstants_categoryWS->add(*categoryWS->var("mu_ZH"));
  muConstants_categoryWS->add(*categoryWS->var("mu_ttH"));
  muConstants_categoryWS->add(*categoryWS->var("mu_BR_gg"));
  
  TIterator *iterMuConst = muConstants_categoryWS->createIterator();
  RooRealVar *currMuConst;
  while ((currMuConst = (RooRealVar*)iterMuConst->Next())) {
    currMuConst->setConstant();
  }
  
  categoryWS->defineSet("muConstants", *muConstants_categoryWS);
  categoryWS->defineSet("observables", *observable_categoryWS);
  categoryWS->defineSet("nuisanceParameters", *nuisance_categoryWS);
  categoryWS->defineSet("globalObservables", *global_categoryWS);
  
  // Import the observed data set:
  currMassPoints = new DMMassPoints(jobName, "obsData", cateScheme, "FromFile",
				currWS->var("m_yy_"+currCateName));
  RooDataSet *obsData = currMassPoints->getCateDataSet(currCateIndex);
  obsData->SetNameTitle("obsData","obsData");
  
  // Set the background normalization parameter:
  (*categoryWS->var("nBkg_"+currCateName) ).setVal(obsdata->numEntries());
  (*categoryWS->pdf("bkgPdf_"+currCateName)).fitTo(*obsData, Minos(RooArgSet(*nuispara_bkg_categoryWS)));
  (*categoryWS->var("nBkg_"+currCateName)).setVal(obsdata->numEntries());
  nuispara_bkg_categoryWS->Print("v");
  categoryWS->import(*obsData);
  
  // Create a binned observed data set:
  RooRealVar wt("wt","wt",1);
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
  RooDataSet *obsdatabinned = new RooDataSet("obsDataBinned", "obsDataBinned",
					     *obsPlusWt, WeightVar(wt));
  int nBin = h_data->GetNbinsX();
  for (int i_b = 1; i_b < nBin; i_b++) {
    // 240 bins -> 0.25 GeV per bin
    double massVal = h_data->GetBinCenter(i_b);
    double weightVal = h_data->GetBinContent(i_b);
    categoryWS->var("atlas_invMass_"+currCateName)->setVal(massVal);
    wt.setVal(weightVal);
    obsDataBinned->add(RooArgSet(*categoryWS->var("m_yy_"+currCateName),wt),
		       weightVal);
  }
  categoryWS->import(*obsDataBinned);
  
  // Create Asimov mu DM = 1 data:
  createAsimovData(currCateName, categoryWS, obsData, wt, 1);
  // Create Asimov mu_DM = 0 data:
  createAsimovData(currCateName, categoryWS, obsData, wt, 0);
  
  //--------------------------------------//
  // Plot the single-channel fit:
  plotFit(currCateName, categoryWS);
  
  delete h_data;
  return categoryWS;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// spurious_signal:

double DMWorkspace::spurious_signal();
{
  double spurious[10] = { 1.0. 1.0, 1.0, 1.0, 1.0 };// per fb-1
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
void DMWorkspace::NPmaker(const char* varName, double setup[4],
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
    
    RooRealVar* var = new RooRealVar(Form("nuisPar_%s",varNameNP.Data()),
				     Form("nuisPar_%s",varNameNP.Data()),
				     varName,0,-5,5);
    RooRealVar* beta_var = new RooRealVar(Form("beta_%s",varName.Data()), 
					  Form("beta_%s",varName.Data()),
					  beta);
    RooProduct* varXBeta = new RooProduct(Form("%s_times_beta",varName.Data()),
					  Form("%s_times_beta",varName.Data()),
					  RooArgSet(*var,*beta_var));
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
    
    std::cout << "  parameter has a Bif. Gauss constraint term" << std::endl;
    std::cout << "  asymmetric factor is " << valA << std::endl;
    
    TString valLogKappa = Form("%f",sqrt(log(1+pow(sigma,2))));
    TString valA = Form("%f",fabs(sigma/sigmaLow)); 
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
  nuisParams->add(*w->var(Form("nuisPar_%s",varName.Data())));
  constraints->add(*w->pdf(Form("constrPdf_%s",varName.Data())));
  globalObs->add(*w->var(Form("globOb_%s",varName.Data())));
  expected->add(*w->function(Form("expected_%s",varName.Data())));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// shapeNPmaker:

void DMWorkspace::shapeNPmaker( const char* varNameNP, const char* proc, double setup[5], RooArgSet *&nuisParams, RooArgSet *&constraints, RooArgSet *&globalObs, RooArgSet *&expected )
{
  // The "shape" NP maker will give variables used in parameterization process dependent name, but keep the same name for the nuisance parameter, and the global observables.
  
  double sigma    = setup[0];
  double sigmaLow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  double nonATLAS = setup[4];
  TString varName = (TString)varNameNP + (TString)proc;
  
  RooWorkspace* w = new RooWorkspace(varName);
  //----------------------------------------//
  // Asymmetric uncertainty:
  if( sigmaLow > 0 )
  {
    cout << " Set up nuisance parameter for an asymmetric uncertainty " << endl;
    RooRealVar* var = new RooRealVar(varNameNP,varNameNP,0,-5,5);
    if( nonATLAS != 0 )
    { 
      TString atlasNPname = (TString)"nuisPar_"+varNameNP;
      var->SetName(atlasNPname); 
      var->SetTitle(atlasNPname); 
    }
    
    RooRealVar* beta_var = new RooRealVar((TString)"beta_"+varName,(TString)"beta_"+varName,beta);
    RooProduct* varXBeta = new RooProduct(varName+(TString)"_times_beta",varName+(TString)"_times_beta",RooArgSet(*var,*beta_var));
    vector<double> sVarHi, sVarLo;
    sVarHi.push_back( 1+sigma );
    sVarLo.push_back( 1-sigmaLow );
    RooArgList nuisList(*varXBeta);
    RooStats::HistFactory::FlexibleInterpVar expVar("expected_"+(TString)varName,"expected_"+(TString)varName,nuisList,nominal,sVarLo,sVarHi);
    w->import(expVar);
    
    if( nonATLAS == 0 )
    {
      cout << " Nuisance parameter is shared between ATLAS and CMS " << endl;
      w->factory((TString)"RooGaussian::constrPdf_"+(TString)varNameNP+(TString)"(R_"+(TString)varNameNP+(TString)"[0,-5,5],"+(TString)varNameNP+(TString)",1)");
    }
    else
    {
      w->factory((TString)"RooGaussian::constrPdf_"+(TString)varNameNP+(TString)"(globOb_"+(TString)varNameNP+(TString)"[0,-5,5],nuisPar_"+(TString)varNameNP+(TString)",1)");
    }
  }
  //----------------------------------------//
  // Gaussian uncertainty:
  else if( sigmaLow == -999 )
  {
    cout << " Set up nuisance parameter with a Gaussian constraint term " << endl;
    TString valSigma=Form("%f", sigma);
    TString valBeta=Form("%f", beta);
    TString valNom=Form("%f", nominal );  
    w->factory((TString)"sum::expected_"+(TString)varName+(TString)"(nominal_"+(TString)varName+"["+valNom+(TString)"]  , prod::uncer_"+(TString)varName+(TString)"( prod::"+varName+(TString)"_times_beta(nuisPar_"+(TString)varNameNP+(TString)"[ 0 ,-5 , 5 ] ,beta_"+varName+(TString)"["+valBeta+(TString)"]), sigma_"+(TString)varName+(TString)"["+valSigma+(TString)" ]))");
    w->factory("RooGaussian::constrPdf_"+(TString)varNameNP+(TString)"(globOb_"+(TString)varNameNP+(TString)"[0,-5,5],nuisPar_"+(TString)varNameNP+(TString)",1)");
  }
  //----------------------------------------//
  // Other case?
  else
  {
    TString valBeta=Form("%f", beta);
    TString valLogKappa=Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    TString valNom=Form("%f", nominal );
    w->factory((TString)"atlas_valLogKappa_"+(TString)varName+"["+(TString)valLogKappa+(TString)"]") ;
    w->factory("RooExponential::atlas_expTerm_"+(TString)varName+"(prod::"+varName+(TString)"_times_beta(nuisPar_"+(TString)varNameNP+(TString)"[ 0 , -5 , 5 ], beta_"+varName+(TString)"["+valBeta+(TString)"]),atlas_valLogKappa_"+(TString)varName+")");
    w->factory((TString)"prod::expected_"+(TString)varName+"(atlas_expTerm_"+(TString)varName+",nominal_"+(TString)varName+"["+(TString)valNom+(TString)"])");
    
    if( nonATLAS != 0 )
      w->factory("RooGaussian::constrPdf_"+(TString)varNameNP+"(globOb_"+(TString)varNameNP+"[0,-5,5],nuisPar_"+(TString)varNameNP+",1)");
    else
    {
      cout << " Set up constraint term for " << varNameNP << endl;
      w->factory("RooGaussian::constrPdf_"+(TString)varNameNP+"(R_"+(TString)varNameNP+"[0,-5,5],"+(TString)varNameNP+",1)");
    }
  }
  
  // declare the NP and constraint term only when it's not declared in the workspace to avoid duplication.   
  if( nonATLAS == 0 ) nuisParams->add(*w->var(varNameNP));
  else nuisParams->add(*w->var("nuisPar_"+(TString)varNameNP));
  constraints->add(*w->pdf("constrPdf_"+(TString)varNameNP));
  if( nonATLAS == 0 ) globalObs->add(*w->var("R_"+(TString)varNameNP));
  else globalObs->add(*w->var("globOb_"+(TString)varNameNP));
  expected->add(*w->function("expected_"+(TString)varName));
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CreateAsimovData:

void DMWorkspace::createAsimovData( TString currCateName, RooWorkspace* categoryWS, RooDataSet *obsdata, RooRealVar wt, double xmin, double DMMyyRangeHi, int epsilon, TString option )
{
  TString spin = ( epsilon == 1 ) ? "0p" : "2p";
  cout << "CreateAsimovData( " << spin << " )" << endl;
  
  int npoints_Asimov = 275;
  
  // This is the dataset to be returned:
  RooDataSet *AsimovData = new RooDataSet( Form("asimovdatabinned%s",spin.Data()), Form("asimovdatabinned%s",spin.Data()), RooArgSet(*categoryWS->var("atlas_invMass_"+currCateName),wt), WeightVar(wt) );
    
  // Load the PDF from the workspace:
  //categoryWS->Print("v");
  RooAbsPdf *current_pdf = (RooAbsPdf*)(categoryWS->pdf("modelSB_"+currCateName));
  double initial_epsilon = (categoryWS->var("epsilon"))->getVal();
  (categoryWS->var("epsilon"))->setVal(epsilon);
  double initial_mu;
  if( option.Contains("decorrmu") )
  {
    initial_mu = (categoryWS->var(Form("mu_%s",currCateName.Data())))->getVal();
    (categoryWS->var(Form("mu_%s",currCateName.Data())))->setVal(1.0);
  }
  else
  {
    if( currCateName.Contains("7TeV") )
    {
      initial_mu = (categoryWS->var("mu_7TeV"))->getVal();
      (categoryWS->var("mu_7TeV"))->setVal(1.0);
    }
    else if( currCateName.Contains("8TeV") )
    {
      initial_mu = (categoryWS->var("mu_8TeV"))->getVal();
      (categoryWS->var("mu_8TeV"))->setVal(1.0);
    }
  }
  
  // use fit result or integrals and sidebands to get the estimate of the background:
  double total_BkgEvents = obsdata->sumEntries();
  double width = ( DMMyyRangeHi - xmin ) / ((double)npoints_Asimov);
  
  // loop over the number of asimov points:
  double count_Asimov = 0.0;
  for( int i_p = 0; i_p < npoints_Asimov; i_p++ )
  {
    double mass_value = DMMyyRangeLo + ( 0.5 * width ) + ( width * (double)i_p );
    (categoryWS->var("atlas_invMass_"+currCateName))->setRange("range_Integral", mass_value-(0.5*width), mass_value+(0.5*width));
    RooAbsReal *integral = (RooAbsReal*)current_pdf->createIntegral(RooArgSet(*categoryWS->var("atlas_invMass_"+currCateName)), NormSet(*categoryWS->var("atlas_invMass_"+currCateName)), Range("range_Integral"));
    double weight_value = total_BkgEvents * integral->getVal();
    count_Asimov += weight_value;
    (categoryWS->var("atlas_invMass_"+currCateName))->setVal(mass_value);
    wt.setVal(weight_value);
    AsimovData->add( RooArgSet( *categoryWS->var("atlas_invMass_"+currCateName), wt ), weight_value );
  }
  if( fabs((count_Asimov-obsdata->sumEntries())/count_Asimov) > 0.04 ){ cout << "Bad Asimov Data: D=" << obsdata->sumEntries() << " A=" << count_Asimov << endl; exit(0); }
  categoryWS->import(*AsimovData);
  (categoryWS->var("epsilon"))->setVal(initial_epsilon);

  if( option.Contains("decorrmu") ) (categoryWS->var(Form("mu_%s",currCateName.Data())))->setVal(initial_mu);
  else
  {
    if( currCateName.Contains("7TeV") )(categoryWS->var("mu_7TeV"))->setVal(initial_mu);
    else if( currCateName.Contains("8TeV") )(categoryWS->var("mu_8TeV"))->setVal(initial_mu);
  }
}

/**
   Plot the fits produced by the specified model.
   @param plotOptions - options for what fits to plot etc.
   @returns void
*/
void DMWorkspace::plotFit(TString plotOptions)
{
  cout << "plotBackgroundOnlyFit( " << currCateName << " )" << endl;
  TCanvas *c = new TCanvas();
  RooPlot* frame =  (*categoryWS->var("atlas_invMass_"+currCateName)).frame(55);
  categoryWS->data("obsdata")->plotOn(frame);
  (*categoryWS->pdf("model_"+currCateName)).plotOn(frame, LineColor(2));
  (*categoryWS->pdf("model_"+currCateName)).plotOn(frame,Components( (*categoryWS->pdf("bkgPdf_"+currCateName)) ) , LineColor(4));
  double chi2 = frame->chiSquare() ;
  frame->SetYTitle("Events / GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  frame->Draw();
  
  TLatex lresult3;
  lresult3.SetNDC();
  lresult3.SetTextColor(1);
  lresult3.DrawLatex(0.5,0.78, currCateName);
  
  system(Form("mkdir -vp %s/figures/",outputDir.Data()));
  PrintCanvas(c, Form("%s/figures/data_fit_%s",outputDir.Data(),currCateName.Data()));
  delete c;
}

/**
   Plots the values of the nuisance parameters in a fit.
   @param plotOptions - options for what fit parameters to plot etc.
   @returns void
*/
void DMWorkspace::plotNuisParams(TString plotOptions)
{
  TCanvas *can = new TCanvas("can","can",2500,1800);
  can->cd();
  can->SetBottomMargin(0.5);
  int index = 0;
  int number_params = nuis.getSize();
  TH1F *h_nuis = new TH1F("h_nuis","h_nuis",number_params,0,number_params);
  TIterator *iter_nuis = nuis.createIterator();
  RooRealVar* parg_nuis = NULL;
  while( (parg_nuis = (RooRealVar*)iter_nuis->Next()) )
  {
    TString name = parg_nuis->GetName();
    double value = parg_nuis->getVal();
    double error = parg_nuis->getError();
    if( !name.Contains("nBkg") && !name.Contains("p0") && !name.Contains("p1") && !name.Contains("p2") )
    {
      index++;
      h_nuis->SetBinContent( index, value );
      h_nuis->SetBinError( index, error );
      h_nuis->GetXaxis()->SetBinLabel( index, name );
    }
  }
  h_nuis->SetLineColor(kBlack);
  h_nuis->GetYaxis()->SetTitle("Nuisance Parameter Pull (#sigma)");
  h_nuis->GetYaxis()->SetTitleOffset(0.3);
  h_nuis->GetXaxis()->SetRangeUser(-1, index+1);
  h_nuis->GetYaxis()->SetRangeUser(-2, 2);
  h_nuis->Draw();
  TBox *b = new TBox( 0, -1, index, 1 );
  b->SetFillColor(kGreen);
  b->SetLineColor(kGreen+2);
  b->Draw("SAME");
  TLine *l = new TLine( 0, 0, index, 0 );
  l->SetLineColor(kRed);
  l->SetLineWidth(2);
  l->SetLineStyle(2);
  l->Draw("SAME");
  h_nuis->Draw("same");
  can->Print( Form("%s/figures/nuisparams_%s.eps",outputDir.Data(),signal_type.Data()) );
  can->Print( Form("%s/figures/nuisparams_%s.png",outputDir.Data(),signal_type.Data()) );
  can->Print( Form("%s/figures/nuisparams_%s.C",outputDir.Data(),signal_type.Data()) );
  delete can;
}
