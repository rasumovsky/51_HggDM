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
  for (int i_c = 0; i_c < nCategories; i_c++) {
    TString currCateName = Form("%s_%d",newCateScheme.Data(),i_c);
    cateNames.push_back(currCateName);
    std::cout << "  \t" << currCateName << std::endl;
  }
  std::cout << "Luminosity at 13 TeV:" << analysisLuminosity << std::endl;
  std::cout << "........................................" << std::endl;
  
  //--------------------------------------//
  // Read tables of ESS and Res and store values:
  //ess_tool = new ESSReader( file_name_ESS_values, nCategories);
  //res_tool = new ResReader( file_name_Res_values, nCategories);
  //ss_tool  = new SigShapeReader( file_name_SS_values, nCategories);
  sigParam = new DMSigParam(jobName, cateScheme, "FromFile");
  massPoints = new DMMassPoints(jobName, "data", cateScheme, "FromFile");
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // Everything for simultaneous fit:
  vector<string> catName;
  RooWorkspace* cateWS[nCategories];
  RooCategory* categories = new RooCategory(Form("categories_%s",
						 newCateScheme.Data()),
					    Form("categories_%s",
						 newCateScheme.Data()));
  
  combination = new RooWorkspace("combination");
  combination->importClassCode();
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
    
    cateWS[i_c] = newCategoryWS(cateNames[i_c]);
    
    categories->defineType(cateNames[i_c]);
    
    // Add PDF and parameters from category to global collections:
    combinedPdf.addPdf(*cateWS[i_c]->pdf("model_"+CN[i_c]),CN[i_c]);
    nuisanceParameters->add(*cateWS[i_c]->set("nuisanceParameters"));
    globalObservables->add(*cateWS[i_c]->set("globalObservables"));
    observables->add(*cateWS[i_c]->set("observables"));
    
    // Retrieve the datasets produced for each category:
    dataMap[catName[i_c]] = (RooDataSet*)w[i_c]->data("obsData");
    dataMapBinned[catName[i_c]] = (RooDataSet*)w[i_c]->data("obsDataBinned");
    dataMapAsimov[catName[i_c]] = (RooDataSet*)w[i_c]->data("asimovBinned");
  }
  
  // Import PDFs and parameters to combined workspace:
  combination->import(combinedPdf);
  combination->defineSet("nuisanceParameters",*nuisanceParameters);
  combination->defineSet("observables",*observables);
  combination->defineSet("globalObservables",*globalObservables);
  combination->defineSet("poi",RooArgSet(*combination->var("mu")));   
  
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
 
  combination->import(*obsData);
  combination->import(*obsDataBinned);
  combination->import(*asimovData);
  
  // Define the ModelConfig:
  mconfig = new ModelConfig("mconfig",combination);
  mconfig->SetPdf(*combination->pdf("combinedPdf"));
  mconfig->SetObservables(*combination->set("observables"));
  mconfig->SetParametersOfInterest((*combination->set("poi")));
  mconfig->SetNuisanceParameters((*combination->set("nuisanceParameters")));
  mconfig->SetGlobalObservables((*combination->set("globalObservables")));
  combination->import(*mconfig);
  
  // Start profiling the data:
  cout << "Start profiling data" << endl;
  RooRealVar *poi = (RooRealVar*)mconfig->GetParametersOfInterest()->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*mconfig->GetNuisanceParameters());
  RooArgSet* globs = (RooArgSet*)mconfig->GetGlobalObservables();
  poiAndNuis->add(*poi);
  poiAndNuis->Print();
  combination->saveSnapshot("paramsOrigin",*poiAndNuis);
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
    resMu1 = pdf->fitTo(*combination->data("asimovMu1"), PrintLevel(0),
			Save(true));
  }
  else {
    resMu1 = pdf->fitTo(*combination->data("obsData"), PrintLevel(0),
			Save(true));
  }
  
  // Track whether all fits converge:
  if (resMu1->status() != 0) AllGoodFits = false;
  
  double nllMu1 = resMu1->minNll();
  combination->saveSnapshot("paramsProfileMu1",*poiAndNuis);
  
  // Plots of invariant mass and nuisance parameters:
  if (!option.Contains("noplot")) {
    if (options.Contains("fitasimov")) {
      PlotFinalFits(combination,"asimovMu1","mu1");
    }
    else {
      PlotFinalFits(combination,"obsData","mu1");
    }
    if (!options.Contains("nosys")) {
      PlotNuisParams(*combination->set("ModelConfig_NuisParams"), "mu1");
    }
  }
  
  //----------------------------------------//
  // Profile and save snapshot for mu=0 fixed:
  std::cout << "\nProfile mu = 0:" << std::endl;
  combination->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  poi->setVal(0.0);
  poi->setConstant(true);
  RooFitResult* resMu0;
  if (options.Contains("fitasimov")) {
    resMu0 = pdf->fitTo(*combination->data("asimovMu1"), PrintLevel(0),
			Save(true));
  }
  else {
    resMu0 = pdf->fitTo(*combination->data("obsData"), PrintLevel(0),
			Save(true));
  }
  
  // Track whether all fits converge:
  if (resMu0->status() != 0) AllGoodFits = false;
  
  double nllMu0 = resMu0->minNll();
  combination->saveSnapshot("paramsProfileMu0",*poiAndNuis);

  // Plots of invariant mass and nuisance parameters:
  if (!options.Contains("noplot")) {
    if (options.Contains("fitasimov")) {
      PlotFinalFits(combination, "asimovMu1", "mu0");
    }
    else {
      PlotFinalFits(combination, "obsData", "mu0");
    }
    if (!options.Contains("nosys")) {
      PlotNuisParams(*combination->set("ModelConfig_NuisParams"), "mu0");
    }
  }
  
  //----------------------------------------//
  // Profile and save snapshot for mu floating:
  std::cout << "\nProfile mu floating:" << std::endl;
  combination->loadSnapshot("paramsOrigin");
  statistics::constSet(poiAndNuis, false);
  statistics::constSet(globs, true);
  poi->setVal(0.0);
  poi->setConstant(false);
  RooFitResult* resMuFree;
  if (options.Contains("fitasimov")) {
    resMuFree = pdf->fitTo(*combination->data("asimovMu1"), PrintLevel(0),
			    Save(true));
  }
  else {
    resMuFree = pdf->fitTo(*combination->data("obsData"), PrintLevel(0), 
			    Save(true));
  }
  
  // Track whether all fits converge:
  if (resMsuFree->status() != 0) AllGoodFits = false;
  
  double nllMuFree = resMuFree->minNll();
  double profiledMuValue = poi->getVal();
  combination->saveSnapshot("paramsProfile_muFree",*poiAndNuis);
 
  // Plots of invariant mass and nuisance parameters:
  if (!option.Contains("noplot")) {
    if (option.Contains("fitasimov")) {
      PlotFinalFits(combination, "asimovMu1", "mufree");
    }
    else {
      PlotFinalFits(combination, "obsData", "mufree");
    }
    if (!option.Contains("nosys")) {
      PlotNuisParams(*combination->set("ModelConfig_NuisParams"), "mufree");
    }
  }
  
  combination->loadSnapshot("paramsOrigin");
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
  combination->writeToFile(Form("%s/workspaceDM.root",outputDir.Data()));
}

/**
   Create the workspace for a single analysis category.
   @param currCategory
*/
RooWorkspace* DMWorkspace::newChannelWS(TString currCategory) {
  
  // The bools that control the systematic uncertainties:
  bool inclusive = currCategory == "inclusive";
  bool channel_constraints_attached = (currCategory == leadingchannel);
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
  RooWorkspace* currWS = new RooWorkspace(currCategory);
  
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
  RooArgSet *expected_spin0p = new RooArgSet();
  RooArgSet *expected_spin2p = new RooArgSet();
  RooArgSet *expected_shape = new RooArgSet();
  RooArgSet *expected_bias = new RooArgSet();
  RooArgSet *expected_proc_ggF = new RooArgSet();
  RooArgSet *expected_non_ggF = new RooArgSet();
  RooArgSet *expected_proc_VHttH = new RooArgSet();
  RooArgSet *expected_proc_VBF = new RooArgSet();
  RooArgSet *expected_proc_WH = new RooArgSet();
  RooArgSet *expected_proc_ZH = new RooArgSet();
  RooArgSet *expected_proc_ttH = new RooArgSet();
  
  // an energy-dependent index, 0 inclusive, 1-11 categories:
  int cateIndex = CatNameToBinNumber(currCategory);
  
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
    int bin_index = CatNameToBinNumber( currCategory );
    int number_SS_sources = ss_tool->GetNumberOfSources(energy);
    // loop over ss sources.
    for( int i_s = 0; i_s < number_SS_sources; i_s++ )
    {
      TString current_SS_source_name = ss_tool->GetNameOfSource( i_s, energy );
      TString ss_np_name = Form("shape_%s",current_SS_source_name.Data());
      double current_ss_value = ss_tool->GetValue( current_SS_source_name, bin_index, energy );
      int current_ss_sign = ss_tool->GetSign( current_SS_source_name, bin_index, energy );
      
      // Asymmetric migration uncertainties:
      double setup_ss_current[5] = {current_ss_value, 0, current_ss_sign, 1, 1};
      if( current_SS_source_name.Contains("_up") )
      {
	TString current_SS_source_name_down = ss_tool->GetNameOfSource( i_s+1, energy );
	double current_ss_value_down = ss_tool->GetValue( current_SS_source_name_down, bin_index, energy );
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
    double ss_events = spurious_signal(currCategory);
    double setup_bias[5] = {ss_events, -999, 1, 0, 1}; //Gaussian constraint
    NPmaker("bias", setup_bias, *&nuisParamsUncorrelated, *&constraints_bias,
	    *&globalObs, *&expected_bias);
  }
  else currWS->factory("atlas_expected_bias[0]");
  
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
      setupRes[0] = res_tool->GetValue(currResSource, cateIndex);
      setupRes[2] = res_tool->GetSign(currResSource, cateIndex);
      
      // resolution on the inclusive shape:
      shapeNPmaker(currResName, "_inc", setupRes, *&nuisParams, *&constraints,
		   *&globalObs, *&expected_shape);
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
      setupESS[0] = ess_tool->GetValue(currESSSource, cateIndex);
      setupESS[2] = ess_tool->GetSign(currESSSource, cateIndex);
      NPmaker(currESSName, setupESS, *&nuisParams, *&constraints, *&globalObs,
	      *&expected_shape);
    }
  }
  
  //--------------------------------------//
  // Parameters of interest (POIs):
  RooRealVar *mu = new RooRealVar("mu","mu",1,-100,100);
  RooRealVar *mu_BR_gg = new RooRealVar("mu_BR_gg","mu_BR_gg",1,-100,100);
  RooRealVar *mu_ggF = new RooRealVar("mu_ggF","mu_ggF",1,-100,100);
  RooRealVar *mu_VBF = new RooRealVar("mu_VBF","mu_VBF",1,-100,100);
  RooRealVar *mu_WH = new RooRealVar("mu_WH","mu_WH",1,-100,100);
  RooRealVar *mu_ZH = new RooRealVar("mu_ZH","mu_ZH",1,-100,100);
  RooRealVar *mu_VH = new RooRealVar("mu_VH","mu_VH",1,-100,100);
  RooRealVar *mu_tH = new RooRealVar("mu_tH","mu_tH",1,-100,100);
  RooRealVar *mu_ttH = new RooRealVar("mu_ttH","mu_ttH",1,-100,100);
  RooRealVar *mu_VBFVH = new RooRealVar("mu_VBFVH","mu_VBFVH",1,-100,100);
  expected->add(RooArgSet(*mu, *mu_BR_gg));
  expected_proc_ggF->add(RooArgSet(*mu_ggF, *mu_tH));
  expected_proc_VBF->add(RooArgSet(*mu_VBF, *mu_VBFVH));
  expected_proc_WH->add(RooArgSet(*mu_WH, *mu_VH, *mu_VBFVH));
  expected_proc_ZH->add(RooArgSet(*mu_ZH, *mu_VH, *mu_VBFVH));
  expected_proc_ttH->add(RooArgSet(*mu_ttH, *mu_tH));
  
  // Expectation values:
  RooProduct expectation_proc_VHttH("expectation_proc_VHttH",
				    "expectation_proc_VHttH",
				    *expected_proc_VHttH);
  RooProduct expectation_common("expectation_common","expectation_common",
				*expected);
  RooProduct expectation_spin0p("expectation_spin0p","expectation_spin0p",
				*expected_spin0p);
  RooProduct expectation_spin2p("expectation_spin2p","expectation_spin2p",
				*expected_spin2p);
  RooProduct expectation_proc_ggF("expectation_proc_ggF","expectation_proc_ggF",
				  *expected_proc_ggF);
  RooProduct expectation_proc_VBF("expectation_proc_VBF","expectation_proc_VBF",
				  *expected_proc_VBF);
  RooProduct expectation_proc_WH("expectation_proc_WH","expectation_proc_WH",
				 *expected_proc_WH);
  RooProduct expectation_proc_ZH("expectation_proc_ZH","expectation_proc_ZH",
				 *expected_proc_ZH);
  RooProduct expectation_proc_ttH("expectation_proc_ttH","expectation_proc_ttH",
				  *expected_proc_ttH);
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  currWS->import(expectation_proc_VHttH);
  currWS->import(expectation_common);
  //currWS->import(expectation_spin0p);
  //currWS->import(expectation_spin2p);
  currWS->import(expectation_proc_ggF);
  currWS->import(expectation_proc_VBF);
  currWS->import(expectation_proc_WH);
  currWS->import(expectation_proc_ZH);
  currWS->import(expectation_proc_ttH);
  
  currWS->import(*expected_shape);
  currWS->import(*expected_bias);
  
  // Declare the observable m_yy, and the observables set:
  currWS->factory(Form("m_yy[%f,%f]",DMMyyRangeLo,DMMyyRangeHi));
  currWS->defineSet("observables","m_yy");
  
  // Set the observables for the various input classes:
  sigParam->setMassObservable(currWS->("m_yy"));
  massPoints->setMassObservable(currWS->var("m_yy"));
  bkgModel->setMassObservable(currWS->var("m_yy"));
  
  ////////
  ////////
  // Construct the signal and background PDFs:
  currWS->factory("epsilon[0.5,0,1]");
  currWS->factory("sum::epsilon_min_1(plusone[1.], prod::mineps(minusone[-1.],epsilon))");
  signalPdfBuilder(*&currWS, essList, resList, "_inc");
  backgroundPdfBuilder(*&currWS, *&nuisParamsBkg, currCategory);
  
  // Add background parameters to uncorrelated collection:
  nuisParamsUncorrelated->add(*nuisParamsBkg);
  
  // build the signal normalization
  TString nSM_0p = Form("%f", value_0p[1]);
  TString nSM_2p = Form("%f", value_2p[1]);
  cout << "  nSM_0p for " << currCategory << " = " << nSM_0p << endl;
  cout << "  nSM_2p for " << currCategory << " = " << nSM_2p << endl;
  
  
  // Normalization for each process follows such pattern:
  // mu*isEM*lumi*migr => expectation_common
  w->factory((TString)"prod::atlas_nsig_spin0p(atlas_nSM_spin0p["+nSM_0p+(TString)"],expectation_common,expectation_spin0p,epsilon)");
  w->factory((TString)"prod::atlas_nsig_spin2p(atlas_nSM_spin2p["+nSM_2p+(TString)"],expectation_common,expectation_spin2p,epsilon_min_1)");
  w->factory("SUM::modelSB(atlas_nsig_spin0p*signalPdf_spin0p,atlas_nsig_spin2p*signalPdf_spin2p,atlas_expected_bias*signalPdf_inc,atlas_nbkg*bkgPdf)");
  w->Print();
  
    
  if( currCategory == leadingchannel )
  {
    constraints->add(*constraints_bias);
    RooProdPdf constraint( "constraint", "constraint", *constraints );
    w->import(constraint);
    w->factory("PROD::model(modelSB,constraint)");
  }
  else
  {
    RooProdPdf constraint( "constraint", "constraint", *constraints_bias );
    w->import(constraint);
    w->factory("PROD::model(modelSB,constraint)");
  }
  
  // Specify the group of nuisance parameters that are correlated between sub-channels.
  // Technically, this is done by sharing the same name for nuisance parameter between sub-channels.
  // Their respective global observables should also share the same name.
  // nuispara should contain all correlated nuisance parameters.
  // all uncorrelated nuisance parameters should be included in nuisParamsUncorrelated.
  TString correlated;
  
  if( m_decorr_mu ) correlated = "epsilon";
  else
  {
    if( currCategory.Contains("8TeV") ) correlated = "mu_8TeV,epsilon,mu_BR_gg_8TeV";
    else correlated = "mu_7TeV,epsilon,mu_BR_gg_7TeV";
  }
  
  // Iterate over nuisance parameters:
  TIterator *iter_nui = nuispara->createIterator();
  RooRealVar* parg_nui = NULL;
  while( (parg_nui=(RooRealVar*)iter_nui->Next()) )
  {
    cout << parg_nui->GetName() << endl;
    correlated = correlated +","+parg_nui->GetName()+",R_"+parg_nui->GetName();
  }
  cout << " For channel " << currCategory << " the following variables will not be renamed : " << correlated << endl;
  
  // sub-channel labeling
  // import the workspace w to another workspace and add currCategory as a suffix to all nodes and variables of w.
  // the correlated nuisance parameters and their respective global observables will not be renamed.
  RooWorkspace* wchannel = new RooWorkspace("wchannel"+currCategory);
  wchannel->import( (*w->pdf("model")), RenameAllNodes(currCategory), RenameAllVariablesExcept(currCategory,correlated), Silence() );
  
  // Adding correlated nuisance parameters to nuisanceParameters:
  //     From nuispara
  RooArgSet* nuisance_wchannel = new RooArgSet();
  iter_nui->Reset();
  cout << " Adding correlated nuisance parameters to nuisanceParameters RooArgSet"<< endl;
  while( (parg_nui=(RooRealVar*)iter_nui->Next()) )
  {
    cout << " Adding variable : " << parg_nui->GetName() << endl;
    cout << (bool)wchannel->obj(parg_nui->GetName()) << endl;
    nuisance_wchannel->add( *(RooRealVar*)wchannel->obj(parg_nui->GetName()) );
  }
  
  // Adding uncorrelated nuisance parameters to nuisanceParameters:
  //   From nuisParamsUncorrelated:
  cout << " Adding uncorrelated nuisance parameters to nuisanceParameters RooArgSet" << endl;
  TIterator *iter_nui_uncorrelated = nuisParamsUncorrelated->createIterator();
  RooRealVar* parg_nui_uncorrelated = NULL;
  while( (parg_nui_uncorrelated = (RooRealVar*)iter_nui_uncorrelated->Next()) )
  {
    TString name_of_nuisance = parg_nui_uncorrelated->GetName()+(TString)"_"+currCategory;
    nuisance_wchannel->add( *(RooRealVar*)wchannel->obj(name_of_nuisance) );
  }
  
  // The following are nps from background pdf, which don't have constraints: 
  //   From nuispara_bkg:
  RooArgSet* nuispara_bkg_wchannel = new RooArgSet();
  TIterator *iter_nui_bkg = nuispara_bkg->createIterator();
  RooRealVar* parg_nui_bkg = NULL;
  while( (parg_nui_bkg = (RooRealVar*)iter_nui_bkg->Next()) )
  {
    TString name_of_parameter = parg_nui_bkg->GetName()+(TString)"_"+currCategory;
    nuispara_bkg_wchannel->add( *wchannel->var(name_of_parameter) );
  }
  
  // Global observables:
  // Global observables only appear in the constraint terms.
  // All constraint terms of correlated nuisance parameters are attached to the pdf of the first subchannel.
  // For those global observables, their names should be the same as those in the w.
  // For other subchannels, only the bias constraint term is attached.
  
  RooArgSet *global_wchannel = new RooArgSet();
  TIterator *iter_global = globobs->createIterator();
  RooRealVar *parg_global;
  while( (parg_global = (RooRealVar*)iter_global->Next()) )//&& (currCategory==channel_constraints_attached) )
  {
    TString name_of_global = parg_global->GetName()+(TString)"_"+currCategory;
    cout << " Channel Name " << currCategory << " getting global observable " << parg_global->GetName() << endl;
    
    if( (bool)wchannel->obj(name_of_global) == true )
    {
      global_wchannel->add( *(RooRealVar*)wchannel->obj(name_of_global) );
      wchannel->var(name_of_global)->setConstant();
    }
    else if( (bool)wchannel->obj(parg_global->GetName()) == true )
    {
      global_wchannel->add( *(RooRealVar*)wchannel->obj(parg_global->GetName()) );
      wchannel->var(parg_global->GetName())->setConstant();
    }
  }
  
  RooArgSet *observable_wchannel = new RooArgSet();
  TIterator *iter_observable = w->set("observables")->createIterator();
  RooRealVar *parg_observable;
  while( (parg_observable = (RooRealVar*)iter_observable->Next()) )
  {
    TString name_of_observable = parg_observable->GetName()+(TString)"_"+currCategory;
    if( (bool)wchannel->obj(name_of_observable) == true )
      observable_wchannel->add( *(RooRealVar*)wchannel->obj(name_of_observable) );
    else
      observable_wchannel->add( *(RooRealVar*)wchannel->obj(parg_observable->GetName()) );
  }
  
  RooArgSet* muconstants_wchannel = new RooArgSet();
  if( currCategory.Contains("7TeV") )
  {
    muconstants_wchannel->add(*wchannel->var("mu_ggF_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VBF_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_WH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_ZH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_ttH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_tH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VBFVH_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_muo_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_ele_7TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_BR_gg_7TeV"));
  }
  else if( currCategory.Contains("8TeV") )
  {
    muconstants_wchannel->add(*wchannel->var("mu_ggF_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VBF_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_WH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_ZH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_ttH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_tH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VBFVH_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_muo_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_VH_ele_8TeV"));
    muconstants_wchannel->add(*wchannel->var("mu_BR_gg_8TeV"));
  }
  
  TIterator *iter_muconst = muconstants_wchannel->createIterator();
  RooRealVar* parg_muconst;
  while( (parg_muconst=(RooRealVar*)iter_muconst->Next()) )
  {
    parg_muconst->setConstant();
  }
  
  wchannel->defineSet("muConstants", *muconstants_wchannel);
  wchannel->defineSet("observables", *observable_wchannel);
  wchannel->defineSet("nuisanceParameters", *nuisance_wchannel);
  wchannel->defineSet("globalObservables", *global_wchannel);
  
  //--------------------------------------//
  // Import data set and set up the background related nuisance parameter values
  dataInputDir = Form("%s/%s/mass_points_%iTeV/",master_output.Data(),jobname.Data(),energy);
  RooDataSet *obsdata = RooDataSet::read(dataInputDir+(TString)"mass_"+currCategory+(TString)".txt",RooArgList(*wchannel->var("atlas_invMass_"+currCategory)));
  obsdata->SetNameTitle("obsdata","obsdata");
  (*wchannel->var("atlas_nbkg_"+currCategory) ).setVal(obsdata->numEntries() );
  (*wchannel->pdf("bkgPdf_"+currCategory)).fitTo( *obsdata, Minos(RooArgSet(*nuispara_bkg_wchannel ) ) );
  (*wchannel->var("atlas_nbkg_"+currCategory) ).setVal(obsdata->numEntries() );
  nuispara_bkg_wchannel->Print("v");
  wchannel->import(*obsdata);
  
  //--------------------------------------//
  // Create a binned data set:
  RooRealVar wt("wt","wt",1);
  RooArgSet* obs_plus_wt = new RooArgSet();
  obs_plus_wt->add(wt);
  obs_plus_wt->add(*wchannel->var("atlas_invMass_"+currCategory));
  
  // histogram to store binned data:
  TH1F* h_data = new TH1F("h_data", "", 240, DMMyyRangeLo, DMMyyRangeHi );
  RooArgSet* obs = (RooArgSet*)obsdata->get();
  RooRealVar* xdata = (RooRealVar*)obs->find("atlas_invMass_"+currCategory);
  for( int i = 0; i < obsdata->numEntries(); i++ )
  {
    obsdata->get(i);
    h_data->Fill( xdata->getVal() );
  }
  // fill obsdatabinned dataset with binned data:
  RooDataSet *obsdatabinned = new RooDataSet( "obsdatabinned", "obsdatabinned", *obs_plus_wt, WeightVar(wt) );
  int nbin = h_data->GetNbinsX();
  for( int ibin = 1; ibin < nbin; ibin++ )
  {
    // 240 bins -> 0.25 GeV per bin
    double mass_val = h_data->GetBinCenter(ibin);
    wchannel->var("atlas_invMass_"+currCategory)->setVal( mass_val );
    double weight = h_data->GetBinContent(ibin);
    wt.setVal(weight);
    obsdatabinned->add( RooArgSet(*wchannel->var("atlas_invMass_"+currCategory), wt), weight );
    mass_val += 0.25;
  }
  wchannel->import(*obsdatabinned);
  
  //--------------------------------------//
  // Create a binned Asimov dataset:
  // Create Asimov spin 0+ data:
  CreateAsimovData( currCategory, wchannel, obsdata, wt, DMMyyRangeLo, DMMyyRangeHi, 1, option );
  // Create Asimov spin 2+ data:
  CreateAsimovData( currCategory, wchannel, obsdata, wt, DMMyyRangeLo, DMMyyRangeHi, 0, option );
  
  //--------------------------------------//
  // Plot the single-channel fit:
  plotBackgroundOnlyFit( wchannel, currCategory );
  
  delete h_data;
  return wchannel;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// readinput:

vector<double> DMWorkspace::readinput( TString currCategory, TString signal_name )
{
  paramInputDir = Form("%s/%s/parameterization_%iTeV/FinalSignal/",master_output.Data(),jobName.Data(),CatNameToEnergy(currCategory));
  TString SignalFileName = Form("%s/fitpars_%s_%s.txt",paramInputDir.Data(),currCategory.Data(),signal_name.Data());
  cout << "Reading file " << SignalFileName.Data() << endl;
  ifstream file_to_read(SignalFileName.Data(),ios::in);
  assert(file_to_read);
  double value[9];
  while( !file_to_read.eof() )
    file_to_read >> value[0] >> value[1] >> value[2] >> value[3] >> value[4] >> value[5] >> value[6] >> value[7] >> value[8];
  file_to_read.close();
  
  // The signal yield we use now correspond to 1 fb-1. Scaling to current luminosity.
  if( currCategory.Contains("7TeV") ) value[1] *= luminosity_7TeV;
  else value[1] *= luminosity_8TeV;
  
  vector<double> result;
  for( int i = 0; i < 9; i++ ) result.push_back(value[i]);
  return result;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// signalPdfBuilder:

void DMWorkspace::signalPdfBuilder( RooWorkspace *&w, vector<double> value, vector<TString> ess_parnames, vector<TString> res_parnames, TString procname )
{
  //----------------------------------------//
  // Create list of ess to multiply:
  TString list_of_ess = "";
  for( int i_e = 0; i_e < (int)ess_parnames.size(); i_e++ )
  {
    TString atlas_exp_name_ess = Form("atlas_expected_%s",ess_parnames[i_e].Data());
    if( (bool)w->obj(atlas_exp_name_ess) != true ) w->factory(Form("%s[1]",atlas_exp_name_ess.Data()));
  
    if( i_e < ((int)ess_parnames.size()-1) ) list_of_ess.Append(Form("%s,",atlas_exp_name_ess.Data()));
    else list_of_ess.Append(Form("%s",atlas_exp_name_ess.Data()));
  }
  //----------------------------------------//
  // Create list of res to multiply:
  // one important difference from ESS: it is atlas_expected_mRes+procname, where procname = _inc,...
  TString list_of_res = "";
  for( int i_r = 0; i_r < (int)res_parnames.size(); i_r++ )
  {
    TString atlas_exp_name_res = Form("atlas_expected_%s",res_parnames[i_r].Data());
    // fix this here:
    if( (bool)w->obj(atlas_exp_name_res) != true ) w->factory(Form("%s%s[1]",atlas_exp_name_res.Data(),procname.Data()));
    
    if( i_r < ((int)res_parnames.size()-1) ) list_of_res.Append(Form("%s%s,",atlas_exp_name_res.Data(),procname.Data()));
    else list_of_res.Append(Form("%s%s",atlas_exp_name_res.Data(),procname.Data()));
  }
  
  cout << "Building a signal pdf " << endl;
  TString mHiggs = Form("%f", value[2]);
  TString mResVal = Form("%f", value[3]);
  TString tailAlpha = Form("%f", value[4]);
  TString mTail = Form("%f", value[6]);
  TString sigTail = Form("%f", value[7]/value[3]);
  TString frac = Form("%f", value[8]);

  // Previous code before modifying resolution systematics:
  //w->factory((TString)"RooCBShape::peakPdf"+procname+(TString)"(atlas_invMass , prod::mHiggs"+procname+(TString)"(mHiggs0"+procname+(TString)"["+mHiggs+(TString)"],"+list_of_ess+(TString)") , atlas_expected_mRes"+procname+(TString)", tailAlpha"+procname+(TString)"["+tailAlpha+(TString)"] , 10)");
  //w->factory((TString)"RooGaussian::tailPdf"+procname+(TString)"(atlas_invMass, prod::mTail"+procname+(TString)"(mTail0"+procname+(TString)"["+mTail+(TString)"],"+list_of_ess+(TString)+"), prod::sigTail"+procname+(TString)"(atlas_expected_mRes"+procname+(TString)","+sigTail+"))");
  
  w->factory((TString)"RooCBShape::peakPdf"+procname+(TString)"(atlas_invMass, prod::mHiggs"+procname+(TString)"(mHiggs0"+procname+(TString)"["+mHiggs+(TString)"],"+list_of_ess+(TString)"), prod::mRes"+procname+(TString)"(mRes0"+procname+(TString)"["+mResVal+(TString)"],"+list_of_res+(TString)"), tailAlpha"+procname+(TString)"["+tailAlpha+(TString)"] , 10)");
  w->factory((TString)"RooGaussian::tailPdf"+procname+(TString)"(atlas_invMass, prod::mTail"+procname+(TString)"(mTail0"+procname+(TString)"["+mTail+(TString)"],"+list_of_ess+(TString)+"), prod::sigTail"+procname+(TString)"(mRes0"+procname+(TString)"["+mResVal+(TString)"],"+list_of_res+(TString)","+sigTail+"))");
  // the implementation of sigTail above scales the resolution of the CB component to that of the GA component.
  w->factory((TString)"SUM::signalPdf"+procname+(TString)"(frac"+procname+(TString)"["+frac+(TString)"]*peakPdf"+procname+(TString)",tailPdf"+procname+(TString)")");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// backgroundPdfBuilder:

void DMWorkspace::backgroundPdfBuilder( RooWorkspace *&w, RooArgSet *&nuispara, TString currCategory )
{
  int cate = CatNameToIndex( currCategory );
  
  if( cate == 0 )
  {
    cout << "Building a 4th order Bernstein polynomials background model for category " << cate << endl;
    w->factory((TString)"RooBernstein::bkgPdf(atlas_invMass,{pconst[1],p0[0.1,-10,10],p1[0.1,-10,10],p2[0.1,-10,10],p3[0.1,-10,10]})");
    nuispara->add(*w->var("p0"));
    nuispara->add(*w->var("p1"));
    nuispara->add(*w->var("p2"));
    nuispara->add(*w->var("p3"));
  }
  else if( ( cate >= 1 && cate <= 9 ) || ( cate >= 12 && cate <= 20 ) )// for first 9 cos(theta*) categories, 2nd order exponential poly works:
  {
    cout << "Building exponentiated 2nd order polynomial background model " << endl;
    w->factory("EXPR::bkgPdf('exp(@1*(@0-100)/100.0+@2*(@0-100)*(@0-100)/10000.0 )',atlas_invMass, p0[-0.02,-500000.,-0.00005], p1[-0.25,-1000.5,1000.5])" );
    nuispara->add(*w->var("p0"));
    nuispara->add(*w->var("p1"));
  }
  else// simple exponential works fine
  {
    cout << "Building exponential background model " << endl;
    w->factory("EXPR::bkgPdf('exp(@1*(@0-100)/100.0)',atlas_invMass, p0[-0.02,-500000.,-0.00005])" );
    nuispara->add(*w->var("p0"));
  }
  
  w->factory("atlas_nbkg[500,0,1000000]");
  nuispara->add(*w->var("atlas_nbkg"));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// spurious_signal:

double DMWorkspace::spurious_signal( TString currCategory )
{
  // first entry [0] is inclusive 7 and 8 TeV!
  double spurious[23] = { 1.6, 0.17, 0.12, 0.03, 0.03, 0.05, 0.15, 0.08, 0.20, 0.04, 0.03, 0.03, 0.13, 0.10, 0.02, 0.02, 0.04, 0.1, 0.06, 0.14, 0.02, 0.02, 0.02 };
  int cate = CatNameToIndex( currCategory ); 
  double result = ( currCategory.Contains("8TeV") ) ? spurious[cate]*luminosity_8TeV : spurious[cate]*luminosity_7TeV;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// NPmaker:

void DMWorkspace::NPmaker( const char* varname, double setup[5], RooArgSet *&nuispara, RooArgSet *&constraints, RooArgSet *&globobs, RooArgSet *&expected )
{
  double sigma    = setup[0];
  double sigmalow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  double nonATLAS = setup[4];
  
  RooWorkspace* w = new RooWorkspace(varname);
  if( sigmalow > 0 ) 
  {
    cout << " Set up nuisance parameter for an asymmetric uncertainty " << endl;
    
    RooRealVar* var = new RooRealVar(varname,varname,0,-5,5);
    if( nonATLAS != 0 )
    { 
      TString atlasNPname = (TString)"atlas_"+varname;
      var->SetName(atlasNPname); 
      var->SetTitle(atlasNPname);
    }
    RooRealVar* beta_var = new RooRealVar((TString)"beta_"+varname,(TString)"beta_"+varname,beta);
    RooProduct* var_times_beta = new RooProduct(varname+(TString)"_times_beta",varname+(TString)"_times_beta",RooArgSet(*var,*beta_var));
    vector<double> sigma_var_high, sigma_var_low;
    sigma_var_high.push_back( 1+sigma );
    sigma_var_low.push_back( 1-sigmalow );
    RooArgList nuiList(*var_times_beta);
    RooStats::HistFactory::FlexibleInterpVar atlas_expected_var("atlas_expected_"+(TString)varname,"atlas_expected_"+(TString)varname,nuiList,nominal,sigma_var_low,sigma_var_high);
    w->import(atlas_expected_var);
    if( nonATLAS == 0 )
    {
      cout << " Nuisance parameter is shared between ATLAS and CMS " << endl;
      w->factory((TString)"RooGaussian::atlas_nui_"+(TString)varname+(TString)"(R_"+(TString)varname+(TString)"[0,-5,5],"+(TString)varname+(TString)",1)");
    }
    else
    {
      w->factory((TString)"RooGaussian::atlas_nui_"+(TString)varname+(TString)"(R_atlas_"+(TString)varname+(TString)"[0,-5,5],atlas_"+(TString)varname+(TString)",1)");
    }
  }
  else if( sigmalow == -999 )
  {
    cout << " Set up nuisance parameter with a Gaussian constraint term, parameter name : " << varname << endl;
    TString sigma_value=Form("%f", sigma);
    TString beta_value=Form("%f", beta);
    TString nominal_value=Form("%f", nominal);
    w->factory((TString)"sum::atlas_expected_"+(TString)varname+(TString)"(nominal_"+(TString)varname+"["+nominal_value+(TString)"] , prod::uncer_"+(TString)varname+(TString)"( prod::"+varname+(TString)"_times_beta(atlas_"+(TString)varname+(TString)"[ 0 , -5 , 5 ] ,beta_"+varname+(TString)"["+beta_value+(TString)"]), sigma_"+(TString)varname+(TString)"["+sigma_value+(TString)" ]))");
    w->factory("RooGaussian::atlas_nui_"+(TString)varname+(TString)"(R_atlas_"+(TString)varname+(TString)"[0,-5,5],atlas_"+(TString)varname+(TString)",1)");
    
  }
  else if( sigmalow<0 && sigmalow != -999)
  {
    TString beta_value=Form("%f", beta);
    TString log_kappa_value=Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    TString nominal_value=Form("%f", nominal );
    TString avalue=Form("%f", fabs(sigma/sigmalow) );
    
    cout << " Set up nuisance parameter with a Bifuricated Gaussian constraint term, parameter name : " << varname << endl;
    cout << " The asymmetric factor is " << avalue<< endl;
    w->factory((TString)"atlas_log_kappa_value_"+(TString)varname+"["+(TString)log_kappa_value+(TString)"]") ;
    if( nonATLAS == 0 )
    {
      w->factory("RooExponential::atlas_expTerm_"+(TString)varname+"(prod::"+varname+(TString)"_times_beta("+(TString)varname+(TString)"[ 0 , -5 , 5 ], beta_"+varname+(TString)"["+beta_value+(TString)"]),atlas_log_kappa_value_"+(TString)varname+")");}
    else
      w->factory("RooExponential::atlas_expTerm_"+(TString)varname+"(prod::"+varname+(TString)"_times_beta(atlas_"+(TString)varname+(TString)"[ 0 , -5 , 5 ], beta_"+varname+(TString)"["+beta_value+(TString)"]),atlas_log_kappa_value_"+(TString)varname+")");
    
    w->factory((TString)"prod::atlas_expected_"+(TString)varname+"(atlas_expTerm_"+(TString)varname+",nominal_"+(TString)varname+"["+(TString)nominal_value+(TString)"])");
    if( nonATLAS == 0 )
    {
      w->factory((TString)"RooBifurGauss::atlas_nui_"+varname+(TString)"(R_"+varname+(TString)"[0,-5,5],"+varname+(TString)",1,"+avalue+(TString)")");
    }
    else
      w->factory((TString)"RooBifurGauss::atlas_nui_"+varname+(TString)"(R_atlas_"+varname+(TString)"[0,-5,5],"+varname+(TString)",1,"+avalue+(TString)")");
  }
  
  else
  {
    cout << " Set up a nuisance parameter with a logNormal constraint, varname: "<< varname << endl;
    TString beta_value=Form("%f", beta);
    TString log_kappa_value=Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    TString nominal_value=Form("%f", nominal );
    w->factory((TString)"atlas_log_kappa_value_"+(TString)varname+"["+(TString)log_kappa_value+(TString)"]") ;
    if( nonATLAS != 0 )
      w->factory("RooExponential::atlas_expTerm_"+(TString)varname+"(prod::"+varname+(TString)"_times_beta(atlas_"+(TString)varname+(TString)"[ 0 , -5 , 5 ], beta_"+varname+(TString)"["+beta_value+(TString)"]),atlas_log_kappa_value_"+(TString)varname+")");
    else if( nonATLAS == 0 )
      w->factory("RooExponential::atlas_expTerm_"+(TString)varname+"(prod::"+varname+(TString)"_times_beta("+varname+(TString)"[ 0 , -5 , 5 ], beta_"+varname+(TString)"["+beta_value+(TString)"]),atlas_log_kappa_value_"+(TString)varname+")");
    w->factory((TString)"prod::atlas_expected_"+(TString)varname+"(atlas_expTerm_"+(TString)varname+",nominal_"+(TString)varname+"["+(TString)nominal_value+(TString)"])");
    if( nonATLAS != 0 )
      w->factory("RooGaussian::atlas_nui_"+(TString)varname+"(R_atlas_"+(TString)varname+"[0,-5,5],atlas_"+(TString)varname+",1)");
    else
    {
      cout << " Set up constraint term for " << varname << endl;
      w->factory("RooGaussian::atlas_nui_"+(TString)varname+"(R_"+(TString)varname+"[0,-5,5],"+(TString)varname+",1)");
    }
  }
  
  if( nonATLAS == 0 ) nuispara->add(*w->var(varname));
  else nuispara->add(*w->var("atlas_"+(TString)varname));
  cout << " Now, adding constraint term " << "atlas_nui_"<<varname<< endl;
  constraints->add(*w->pdf("atlas_nui_"+(TString)varname));
  if( nonATLAS ==0 ) globobs->add(*w->var("R_"+(TString)varname));
  else globobs->add(*w->var("R_atlas_"+(TString)varname));
  expected->add(*w->function("atlas_expected_"+(TString)varname));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// shapeNPmaker:

void DMWorkspace::shapeNPmaker( const char* varnameNP, const char* proc, double setup[5], RooArgSet *&nuispara, RooArgSet *&constraints, RooArgSet *&globobs, RooArgSet *&expected )
{
  // The "shape" NP maker will give variables used in parameterization process dependent name, but keep the same name for the nuisance parameter, and the global observables.
  
  double sigma    = setup[0];
  double sigmalow = setup[1];
  double beta     = setup[2];
  double nominal  = setup[3];
  double nonATLAS = setup[4];
  TString varname = (TString)varnameNP + (TString)proc;
  
  RooWorkspace* w = new RooWorkspace(varname);
  //----------------------------------------//
  // Asymmetric uncertainty:
  if( sigmalow > 0 )
  {
    cout << " Set up nuisance parameter for an asymmetric uncertainty " << endl;
    RooRealVar* var = new RooRealVar(varnameNP,varnameNP,0,-5,5);
    if( nonATLAS != 0 )
    { 
      TString atlasNPname = (TString)"atlas_"+varnameNP;
      var->SetName(atlasNPname); 
      var->SetTitle(atlasNPname); 
    }
    
    RooRealVar* beta_var = new RooRealVar((TString)"beta_"+varname,(TString)"beta_"+varname,beta);
    RooProduct* var_times_beta = new RooProduct(varname+(TString)"_times_beta",varname+(TString)"_times_beta",RooArgSet(*var,*beta_var));
    vector<double> sigma_var_high, sigma_var_low;
    sigma_var_high.push_back( 1+sigma );
    sigma_var_low.push_back( 1-sigmalow );
    RooArgList nuiList(*var_times_beta);
    RooStats::HistFactory::FlexibleInterpVar atlas_expected_var("atlas_expected_"+(TString)varname,"atlas_expected_"+(TString)varname,nuiList,nominal,sigma_var_low,sigma_var_high);
    w->import(atlas_expected_var);
    
    if( nonATLAS == 0 )
    {
      cout << " Nuisance parameter is shared between ATLAS and CMS " << endl;
      w->factory((TString)"RooGaussian::atlas_nui_"+(TString)varnameNP+(TString)"(R_"+(TString)varnameNP+(TString)"[0,-5,5],"+(TString)varnameNP+(TString)",1)");
    }
    else
    {
      w->factory((TString)"RooGaussian::atlas_nui_"+(TString)varnameNP+(TString)"(R_atlas_"+(TString)varnameNP+(TString)"[0,-5,5],atlas_"+(TString)varnameNP+(TString)",1)");
    }
  }
  //----------------------------------------//
  // Gaussian uncertainty:
  else if( sigmalow == -999 )
  {
    cout << " Set up nuisance parameter with a Gaussian constraint term " << endl;
    TString sigma_value=Form("%f", sigma);
    TString beta_value=Form("%f", beta);
    TString nominal_value=Form("%f", nominal );  
    w->factory((TString)"sum::atlas_expected_"+(TString)varname+(TString)"(nominal_"+(TString)varname+"["+nominal_value+(TString)"]  , prod::uncer_"+(TString)varname+(TString)"( prod::"+varname+(TString)"_times_beta(atlas_"+(TString)varnameNP+(TString)"[ 0 ,-5 , 5 ] ,beta_"+varname+(TString)"["+beta_value+(TString)"]), sigma_"+(TString)varname+(TString)"["+sigma_value+(TString)" ]))");
    w->factory("RooGaussian::atlas_nui_"+(TString)varnameNP+(TString)"(R_atlas_"+(TString)varnameNP+(TString)"[0,-5,5],atlas_"+(TString)varnameNP+(TString)",1)");
  }
  //----------------------------------------//
  // Other case?
  else
  {
    TString beta_value=Form("%f", beta);
    TString log_kappa_value=Form("%f", sqrt( log( 1+pow(sigma,2)) ) );
    TString nominal_value=Form("%f", nominal );
    w->factory((TString)"atlas_log_kappa_value_"+(TString)varname+"["+(TString)log_kappa_value+(TString)"]") ;
    w->factory("RooExponential::atlas_expTerm_"+(TString)varname+"(prod::"+varname+(TString)"_times_beta(atlas_"+(TString)varnameNP+(TString)"[ 0 , -5 , 5 ], beta_"+varname+(TString)"["+beta_value+(TString)"]),atlas_log_kappa_value_"+(TString)varname+")");
    w->factory((TString)"prod::atlas_expected_"+(TString)varname+"(atlas_expTerm_"+(TString)varname+",nominal_"+(TString)varname+"["+(TString)nominal_value+(TString)"])");
    
    if( nonATLAS != 0 )
      w->factory("RooGaussian::atlas_nui_"+(TString)varnameNP+"(R_atlas_"+(TString)varnameNP+"[0,-5,5],atlas_"+(TString)varnameNP+",1)");
    else
    {
      cout << " Set up constraint term for " << varnameNP << endl;
      w->factory("RooGaussian::atlas_nui_"+(TString)varnameNP+"(R_"+(TString)varnameNP+"[0,-5,5],"+(TString)varnameNP+",1)");
    }
  }
  
  // declare the NP and constraint term only when it's not declared in the workspace to avoid duplication.   
  if( nonATLAS == 0 ) nuispara->add(*w->var(varnameNP));
  else nuispara->add(*w->var("atlas_"+(TString)varnameNP));
  constraints->add(*w->pdf("atlas_nui_"+(TString)varnameNP));
  if( nonATLAS == 0 ) globobs->add(*w->var("R_"+(TString)varnameNP));
  else globobs->add(*w->var("R_atlas_"+(TString)varnameNP));
  expected->add(*w->function("atlas_expected_"+(TString)varname));
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CreateAsimovData:

void DMWorkspace::CreateAsimovData( TString currCategory, RooWorkspace* wchannel, RooDataSet *obsdata, RooRealVar wt, double xmin, double DMMyyRangeHi, int epsilon, TString option )
{
  TString spin = ( epsilon == 1 ) ? "0p" : "2p";
  cout << "CreateAsimovData( " << spin << " )" << endl;
  
  int npoints_Asimov = 275;
  
  // This is the dataset to be returned:
  RooDataSet *AsimovData = new RooDataSet( Form("asimovdatabinned%s",spin.Data()), Form("asimovdatabinned%s",spin.Data()), RooArgSet(*wchannel->var("atlas_invMass_"+currCategory),wt), WeightVar(wt) );
    
  // Load the PDF from the workspace:
  //wchannel->Print("v");
  RooAbsPdf *current_pdf = (RooAbsPdf*)(wchannel->pdf("modelSB_"+currCategory));
  double initial_epsilon = (wchannel->var("epsilon"))->getVal();
  (wchannel->var("epsilon"))->setVal(epsilon);
  double initial_mu;
  if( option.Contains("decorrmu") )
  {
    initial_mu = (wchannel->var(Form("mu_%s",currCategory.Data())))->getVal();
    (wchannel->var(Form("mu_%s",currCategory.Data())))->setVal(1.0);
  }
  else
  {
    if( currCategory.Contains("7TeV") )
    {
      initial_mu = (wchannel->var("mu_7TeV"))->getVal();
      (wchannel->var("mu_7TeV"))->setVal(1.0);
    }
    else if( currCategory.Contains("8TeV") )
    {
      initial_mu = (wchannel->var("mu_8TeV"))->getVal();
      (wchannel->var("mu_8TeV"))->setVal(1.0);
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
    (wchannel->var("atlas_invMass_"+currCategory))->setRange("range_Integral", mass_value-(0.5*width), mass_value+(0.5*width));
    RooAbsReal *integral = (RooAbsReal*)current_pdf->createIntegral(RooArgSet(*wchannel->var("atlas_invMass_"+currCategory)), NormSet(*wchannel->var("atlas_invMass_"+currCategory)), Range("range_Integral"));
    double weight_value = total_BkgEvents * integral->getVal();
    count_Asimov += weight_value;
    (wchannel->var("atlas_invMass_"+currCategory))->setVal(mass_value);
    wt.setVal(weight_value);
    AsimovData->add( RooArgSet( *wchannel->var("atlas_invMass_"+currCategory), wt ), weight_value );
  }
  if( fabs((count_Asimov-obsdata->sumEntries())/count_Asimov) > 0.04 ){ cout << "Bad Asimov Data: D=" << obsdata->sumEntries() << " A=" << count_Asimov << endl; exit(0); }
  wchannel->import(*AsimovData);
  (wchannel->var("epsilon"))->setVal(initial_epsilon);

  if( option.Contains("decorrmu") ) (wchannel->var(Form("mu_%s",currCategory.Data())))->setVal(initial_mu);
  else
  {
    if( currCategory.Contains("7TeV") )(wchannel->var("mu_7TeV"))->setVal(initial_mu);
    else if( currCategory.Contains("8TeV") )(wchannel->var("mu_8TeV"))->setVal(initial_mu);
  }
}

/**
   Plot the fits produced by the specified model.
   @param plotOptions - options for what fits to plot etc.
   @returns void
*/
void DMWorkspace::plotFit(TString plotOptions)
{
  cout << "plotBackgroundOnlyFit( " << currCategory << " )" << endl;
  TCanvas *c = new TCanvas();
  RooPlot* frame =  (*wchannel->var("atlas_invMass_"+currCategory)).frame(55);
  wchannel->data("obsdata")->plotOn(frame);
  (*wchannel->pdf("model_"+currCategory)).plotOn(frame, LineColor(2));
  (*wchannel->pdf("model_"+currCategory)).plotOn(frame,Components( (*wchannel->pdf("bkgPdf_"+currCategory)) ) , LineColor(4));
  double chi2 = frame->chiSquare() ;
  frame->SetYTitle("Events / GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  frame->Draw();
  
  TLatex lresult3;
  lresult3.SetNDC();
  lresult3.SetTextColor(1);
  lresult3.DrawLatex(0.5,0.78, currCategory);
  
  system(Form("mkdir -vp %s/figures/",outputDir.Data()));
  PrintCanvas(c, Form("%s/figures/data_fit_%s",outputDir.Data(),currCategory.Data()));
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
    if( !name.Contains("atlas_nbkg") && !name.Contains("p0") && !name.Contains("p1") && !name.Contains("p2") )
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
