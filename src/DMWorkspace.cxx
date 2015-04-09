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
  
  //--------------------------------------//
  // Read tables of ESS and Res and store values:
  //ess_tool = new ESSReader( file_name_ESS_values, nCategories);
  //res_tool = new ResReader( file_name_Res_values, nCategories);
  //ss_tool  = new SigShapeReader( file_name_SS_values, nCategories);
  
  // Instantiate the signal parameterization class using the observable:
  currSigParam = new DMSigParam(jobName, cateScheme, "FromFile");
  
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
  currSigParam->addSigToCateWS(currWS, essList, resList, "SM", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ggH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "VBF", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "WH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ZH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "ttH", currCateIndex);
  currSigParam->addSigToCateWS(currWS, essList, resList, "bbH", currCateIndex);
    
  // Instantiate the background parameterization class using the observable:
  currBkgModel = new DMBkgModel(jobName, cateScheme, "FromFile",
				currWS->var("m_yy"));
  // Construct the background PDF:
  backgroundPdfBuilder(nuisParamsBkg, currCateName);
  
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
////////// signalPdfBuilder:

/**
   Need to figure out the best interface method for the SigParam class. 
*/
void DMWorkspace::signalPdfBuilder( RooWorkspace *&w, vector<double> value, vector<TString> parNamesESS, vector<TString> parNamesRes, TString procname )
{
  //----------------------------------------//
  // Create list of ess to multiply:
  TString listESS = "";
  for (int i_e = 0; i_e < (int)parNamesESS.size(); i_e++) {
    TString atlas_exp_name_ess = Form("atlas_expected_%s",
				      parNamesESS[i_e].Data());
    if (!(bool)w->obj(atlas_exp_name_ess)) {
      w->factory(Form("%s[1]",atlas_exp_name_ess.Data()));
    }
    
    if (i_e < ((int)parNamesESS.size()-1)) {
      listESS.Append(Form("%s,",atlas_exp_name_ess.Data()));//added comma
    }
    else {
      listESS.Append(Form("%s",atlas_exp_name_ess.Data()));//no comma
    }
  }
  
  //----------------------------------------//
  // Create list of res to multiply:
  // Important difference from ESS: it is atlas_expected_mRes+procname, where procname = _inc,...
  TString listRes = "";
  for (int i_r = 0; i_r < (int)parNamesRes.size(); i_r++) {
    TString atlas_exp_name_res = Form("atlas_expected_%s",parNamesRes[i_r].Data());
    // fix this here:
    if (!(bool)w->obj(atlas_exp_name_res)) {
      w->factory(Form("%s%s[1]",atlas_exp_name_res.Data(),procname.Data()));
    }
    if (i_r < ((int)parNamesRes.size()-1)) {
      listRes.Append(Form("%s%s,",atlas_exp_name_res.Data(),procname.Data()));
    }
    else {
      listRes.Append(Form("%s%s",atlas_exp_name_res.Data(),procname.Data()));
    }
  }
  
  TString mHiggs = Form("%f", value[2]);
  TString mResVal = Form("%f", value[3]);
  TString tailAlpha = Form("%f", value[4]);
  TString mTail = Form("%f", value[6]);
  TString sigTail = Form("%f", value[7]/value[3]);
  TString frac = Form("%f", value[8]);

  // Previous code before modifying resolution systematics:
  //w->factory((TString)"RooCBShape::peakPdf"+procname+(TString)"(atlas_invMass , prod::mHiggs"+procname+(TString)"(mHiggs0"+procname+(TString)"["+mHiggs+(TString)"],"+listESS+(TString)") , atlas_expected_mRes"+procname+(TString)", tailAlpha"+procname+(TString)"["+tailAlpha+(TString)"] , 10)");
  //w->factory((TString)"RooGaussian::tailPdf"+procname+(TString)"(atlas_invMass, prod::mTail"+procname+(TString)"(mTail0"+procname+(TString)"["+mTail+(TString)"],"+listESS+(TString)+"), prod::sigTail"+procname+(TString)"(atlas_expected_mRes"+procname+(TString)","+sigTail+"))");
  
  w->factory((TString)"RooCBShape::peakPdf"+procname+(TString)"(atlas_invMass, prod::mHiggs"+procname+(TString)"(mHiggs0"+procname+(TString)"["+mHiggs+(TString)"],"+listESS+(TString)"), prod::mRes"+procname+(TString)"(mRes0"+procname+(TString)"["+mResVal+(TString)"],"+listRes+(TString)"), tailAlpha"+procname+(TString)"["+tailAlpha+(TString)"] , 10)");
  w->factory((TString)"RooGaussian::tailPdf"+procname+(TString)"(atlas_invMass, prod::mTail"+procname+(TString)"(mTail0"+procname+(TString)"["+mTail+(TString)"],"+listESS+(TString)+"), prod::sigTail"+procname+(TString)"(mRes0"+procname+(TString)"["+mResVal+(TString)"],"+listRes+(TString)","+sigTail+"))");
  // the implementation of sigTail above scales the resolution of the CB component to that of the GA component.
  w->factory((TString)"SUM::signalPdf"+procname+(TString)"(frac"+procname+(TString)"["+frac+(TString)"]*peakPdf"+procname+(TString)",tailPdf"+procname+(TString)")");
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// backgroundPdfBuilder:

// Also need to figure out how to interface with BkgModel class. Probably need to pass it the workspace. 

//void DMWorkspace::backgroundPdfBuilder( RooWorkspace *&w, RooArgSet *&nuispara, TString currCateName )
void DMWorkspace::backgroundPdfBuilder(RooArgSet *nuisParams) {
  
  RooAbsPdf *cateBkgPdf = DMBkgModel->getCateBkgPdf(currCateIndex);
  currWS->add(cateBkgPdf);
  // now add the associated parameters;
  RooArgSet *cateBkgPdfArgs = DMBkgModel->getCateBkgPars(currCateIndex);
  
  // THE LINE BELOW MUST BE FIXED. Maybe get from currWS by name... 
  nuisParams->add();
  
  //w->factory((TString)"RooBernstein::bkgPdf(m_yy,{pconst[1],p0[0.1,-10,10],p1[0.1,-10,10],p2[0.1,-10,10],p3[0.1,-10,10]})");
  //nuisParams->add(*currWS->var("p0"));
  //nuisParams->add(*currWS->var("p1"));
  
  // Declare the background normalziation parameter:
  w->factory("nBkg[100,0,1000000]");
  nuisParams->add(*currWS->var("nBkg"));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// spurious_signal:

double DMWorkspace::spurious_signal();
{
  // first entry [0] is inclusive 7 and 8 TeV!
  double spurious[23] = { 1.6, 0.17, 0.12, 0.03, 0.03, 0.05, 0.15, 0.08, 0.20, 0.04, 0.03, 0.03, 0.13, 0.10, 0.02, 0.02, 0.04, 0.1, 0.06, 0.14, 0.02, 0.02, 0.02 };
  return spurious[currCateIndex] * analysisLuminosity;
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
