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
   text
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
    
  // Make new or load old workspace:
  if (options.Contains("FromFile")) loadWSFromFile();
  else createNewWS();
  return;
}

/**
   text
*/
void DMWorkspace::createNewWS() {
  
  //--------------------------------------//
  // Determine which categories to turn on:
  int nchannel_switches_current = nchannel_switches_7TeV+nchannel_switches_8TeV+1;
  cout << "There are " << nchannel_switches_current << " input flags." << endl;
  cout << "........................................" << endl;
  vector<TString> CN;
  if( channel_switch[0] == 1 ) CN.push_back("inclusive");
  for( int i_f = 1; i_f < nchannel_switches_current; i_f++ )
  {
    int energy = 7; int bin_number = i_f;
    if( i_f > nchannel_switches_7TeV ){ bin_number = i_f - nchannel_switches_7TeV; energy = 8; }
    if( channel_switch[i_f] == 1 ) CN.push_back(Form("costs%i_%iTeV",bin_number,energy));
  }
  // Get the number of channels:
  nchannels = CN.size();
  cout << nchannels << " total 7+8 TeV categories to be included:" << endl;
  for( int i_c = 0; i_c < nchannels; i_c++ ) cout << "    " << CN[i_c].Data() << endl;
  cout << "........................................" << endl;
  cout << "Luminosity_7TeV: " << luminosity_7TeV << " ifb" << endl;
  cout << "Luminosity_8TeV: " << luminosity_8TeV << " ifb" << endl;
  cout << "........................................" << endl;
  cout << "Making workspace for the " << Spin2Type << " alternative signal." << endl;
  cout << " " << endl;
  
  //--------------------------------------//
  // Read tables of ESS and Res and store values:
  ess_tool = new ESSReader( file_name_ESS_values_7TeV, file_name_ESS_values_8TeV, nchannel_switches_7TeV, nchannel_switches_8TeV );
  res_tool = new ResReader( file_name_Res_values_7TeV, file_name_Res_values_8TeV, nchannel_switches_7TeV, nchannel_switches_8TeV );
  ss_tool  = new SigShapeReader( file_name_SS_values_7TeV, file_name_SS_values_8TeV, nchannel_switches_7TeV, nchannel_switches_8TeV );
  
  //--------------------------------------//
  // Initialize classes relevant to workspace:
  // Everything for simultaneous fit:
  vector<string> catName;
  RooWorkspace* w[nchannels];
  RooCategory* channellist = new RooCategory("channellist","channellist");
  RooWorkspace* combination = new RooWorkspace("combination");
  combination->importClassCode();
  RooSimultaneous CombinedPdf("CombinedPdf","",*channellist);
  
  // Sets of parameters:
  RooArgSet* nuisanceParameters = new RooArgSet();
  RooArgSet* globalObservables = new RooArgSet();
  RooArgSet* Observables = new RooArgSet();
  RooArgSet* constraints = new RooArgSet();

  // maps for datasets:
  map<string,RooDataSet*> datasetMap;
  map<string,RooDataSet*> datasetMap_binned;
  
  map<string,RooDataSet*> datasetMap_asimov0p;
  map<string,RooDataSet*> datasetMap_asimov2p;
  
  //--------------------------------------//
  // Loop over channels:
  cout << "Beginning loop over channels." << endl;
  for( int i_c = 0; i_c < nchannels; i_c++ )
  {
    cout << "  Channel name " << CN[i_c] << endl;
    w[i_c] = GenerateSingleChannel( CN[i_c], CN[0], option );
    channellist->defineType( CN[i_c] );
    
    // Only add the model part to simultaneous pdf.
    CombinedPdf.addPdf(*w[i_c]->pdf("model_"+CN[i_c]),CN[i_c]);
    // should constraints be included? for some reason it was commented...
    //constraints->add( (RooArgSet*)(w[i_c]->pdf("modelSB"+CN[i_c])->getAllConstraints() ));
    nuisanceParameters->add(*w[i_c]->set("nuisanceParameters"));
    globalObservables->add(*w[i_c]->set("globalObservables"));
    Observables->add(*w[i_c]->set("Observables"));
    
    catName.push_back((string)CN[i_c]);
    datasetMap[catName[i_c]] = (RooDataSet*)w[i_c]->data("obsdata");
    datasetMap_binned[catName[i_c]] = (RooDataSet*)w[i_c]->data("obsdatabinned");
    datasetMap_asimov0p[catName[i_c]] = (RooDataSet*)w[i_c]->data("asimovdatabinned0p");
    datasetMap_asimov2p[catName[i_c]] = (RooDataSet*)w[i_c]->data("asimovdatabinned2p");
  }
  
  combination->import(CombinedPdf);
  if( option.Contains("decorrmu") )
  {
    for( int i_c = 0; i_c < nchannels; i_c++ )
      nuisanceParameters->add(*combination->var(Form("mu_%s",CN[i_c].Data())));
  }
  else
  {
    nuisanceParameters->add(*combination->var("mu_7TeV"));
    nuisanceParameters->add(*combination->var("mu_8TeV"));
  }
  
  combination->defineSet("nuisanceParameters",*nuisanceParameters);
  combination->defineSet("Observables",*Observables);
  combination->defineSet("globalObservables",*globalObservables);
  combination->defineSet("poi",RooArgSet(*combination->var("epsilon")));   
  
  RooRealVar wt("wt","wt",1);
  RooArgSet *args = new RooArgSet();
  args->add(*Observables);
  args->add(wt);
  RooDataSet* obsData = new RooDataSet("obsData","combined data ",*args, Index(*channellist), Import(datasetMap), WeightVar(wt));
  RooDataSet* obsDatabinned = new RooDataSet("obsDatabinned","combined data binned",*args, Index(*channellist), Import(datasetMap_binned), WeightVar(wt));
  RooDataSet* asimovData0p = new RooDataSet("asimovData0p","combined asimov 0p data",*args, Index(*channellist), Import(datasetMap_asimov0p), WeightVar(wt));
  RooDataSet* asimovData2p = new RooDataSet("asimovData2p","combined asimov 2p data",*args, Index(*channellist), Import(datasetMap_asimov2p), WeightVar(wt));
  combination->import(*obsData);
  combination->import(*obsDatabinned);
  combination->import(*asimovData0p);
  combination->import(*asimovData2p);
  
  // Define the ModelConfig:
  ModelConfig *mconfig = new ModelConfig("mconfig",combination);
  mconfig->SetPdf(*combination->pdf("CombinedPdf"));
  mconfig->SetObservables( *combination->set("Observables"));
  mconfig->SetParametersOfInterest( (*combination->set("poi")) );
  mconfig->SetNuisanceParameters( (*combination->set("nuisanceParameters")) );
  mconfig->SetGlobalObservables( (*combination->set("globalObservables")) );
  combination->import(*mconfig);
  
  //--------------------------------------//
  // Start profiling the data:
  cout << "Start profiling data" << endl;
  RooRealVar *poi = (RooRealVar*)mconfig->GetParametersOfInterest()->first();
  RooArgSet* poiAndNuis=new RooArgSet();
  poiAndNuis->add(*mconfig->GetNuisanceParameters());
  poiAndNuis->add(*poi);
  poiAndNuis->Print();
  combination->saveSnapshot("paramsOrigin",*poiAndNuis);
  RooAbsPdf *pdf = mconfig->GetPdf();
  
  // Profile to spin 0:
  cout << "Profile to spin 0" << endl;
  poi->setVal(1);
  poi->setConstant(true);
  RooFitResult* res0p;
  if( option.Contains("fitasimov0p") ) res0p = pdf->fitTo(*asimovData0p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else if( option.Contains("fitasimov2p") ) res0p = pdf->fitTo(*asimovData2p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else res0p = pdf->fitTo(*obsData, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  
  PlotNuisParams( *combination->set("nuisanceParameters"), "Spin0" );
  
  double nll0p = res0p->minNll();
  combination->saveSnapshot("paramsProfile0p",*poiAndNuis);
  combination->loadSnapshot("paramsOrigin");
  
  // Profile to spin 2:
  cout << "Profile to spin 2" << endl;
  poi->setVal(0);
  poi->setConstant(true);
  RooFitResult* res2p;
  if( option.Contains("fitasimov0p") ) res2p = pdf->fitTo(*asimovData0p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else if( option.Contains("fitasimov2p") ) res2p = pdf->fitTo(*asimovData2p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else res2p = pdf->fitTo(*obsData, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  
  PlotNuisParams( *combination->set("nuisanceParameters"), "Spin2" );
  
  double nll2p = res2p->minNll();
  combination->saveSnapshot("paramsProfile2p",*poiAndNuis);
  combination->loadSnapshot("paramsOrigin");
  
  // Epsilon is free in the fit:
  poi->setVal(0.5);
  poi->setConstant(false);
  // pdf->fitTo(*obsData, Hesse(false), Minos(false), PrintLevel(0));
  if( option.Contains("fitasimov0p") ) pdf->fitTo(*asimovData0p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else if( option.Contains("fitasimov2p") ) pdf->fitTo(*asimovData2p, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  else pdf->fitTo(*obsData, Hesse(false), Minos(false), PrintLevel(0), Save(true));
  
  PlotNuisParams( *combination->set("nuisanceParameters"), "Free" );
  
  combination->saveSnapshot("paramsProfilemix",*poiAndNuis);
  combination->loadSnapshot("paramsOrigin");
  cout << setprecision(10) << nll0p << " " << nll2p << " " << 2*(nll0p-nll2p) << endl;
  combination->writeToFile(Form("%s/combinedWS.root",output_directory.Data()));
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GenerateSingleChannel:

RooWorkspace* DMWorkspace::newChannelWS(TString channelname) {
  // check the energy corresponding to this category:
  int energy = ( channelname.Contains("7TeV") ) ? 7 : 8;
  // get the index corresponding to the category:
  int channel_index = CatNameToIndex( channelname ) - 1;
  
  // mass range:
  double xmin = 105, xmax = 160;
  
  //--------------------------------------//
  // The bools that control the systematic uncertainties:
  bool inclusive = channelname == "inclusive";
  bool channel_constraints_attached = ( channelname == leadingchannel );
  bool m_noess = option.Contains("noess");
  bool m_nores = option.Contains("nores");
  bool m_noss  = option.Contains("noss");
  bool m_nobgm = option.Contains("nobgm");
  bool m_nomig = option.Contains("nomig");
  bool m_nosys = option.Contains("nosys");
  bool m_decorr_mu = option.Contains("decorrmu");
  if( m_nosys ) cout << "GenerateSingleChannel(): ALL systematics = OFF" << endl;
  else
  {
    if( m_noess ) cout << "GenerateSingleChannel(): energy scale systematics = OFF" << endl;
    if( m_nores ) cout << "GenerateSingleChannel(): resolution systematics   = OFF" << endl;
    if( m_noss  ) cout << "GenerateSingleChannel(): shape systematics        = OFF"  << endl;
    if( m_nobgm ) cout << "GenerateSingleChannel(): background systematics   = OFF"  << endl;
    if( m_nomig ) cout << "GenerateSingleChannel(): migration systematics    = OFF"  << endl;
  }
  if( m_decorr_mu ) cout << "GenerateSingleChannel(): Mu in each bin will be decorrelated." << endl;
  bool m_nonorm = true;
  bool s_normalization = false; // this should never be on. just left in case of a radical change...
  if( m_nosys )
  {
    m_noess = true;
    m_nores = true;
    m_noss  = true;
    m_nobgm = true;
    m_nomig = true;
  }
  
  //--------------------------------------//
  // Create the individual channel workspace:
  RooWorkspace* w = new RooWorkspace(channelname);
  
  // Reading signal shape inputs for each individual process
  vector<double> value_shape = readinput( channelname, MCSampleName[0] ); // spin0
  vector<double> value_0p    = readinput( channelname, MCSampleName[0] ); // spin0
  vector<double> value_2p    = readinput( channelname, Spin2Type ); // current spin2
  
  // nuispara:
  RooArgSet *nuispara = new RooArgSet();
  RooArgSet *nuispara_null = new RooArgSet();
  RooArgSet *nuispara_bkg = new RooArgSet();
  RooArgSet *nuispara_uncorrelated = new RooArgSet();
  RooArgSet *nuispara_proc = new RooArgSet();
  // constraints:
  RooArgSet *constraints = new RooArgSet();
  RooArgSet *constraints_null = new RooArgSet();
  RooArgSet *constraints_bias = new RooArgSet();
  RooArgSet *constraints_proc = new RooArgSet();
  // globobs:
  RooArgSet *globobs = new RooArgSet();
  RooArgSet *globobs_null = new RooArgSet();
  RooArgSet *globobs_proc = new RooArgSet();
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

  // double setup[5] = { 0.25, 0, 1, 1, 1 };
  // array setup[5] is used to configure a nuisance parameter
  // [0]    [1]       [2]   [3]      [4]
  // sigma, sigmalow, beta, nominal, nonATLAS
  
  //--------------------------------------//
  // SYSTEMATICS: Normalization (NEVER USE!)
  if( !m_nonorm )
  {
    double setup_lumi[5] = { 0.036, 0, 1, 1, 1 };
    NPmaker( "lumi", setup_lumi, *&nuispara, *&constraints, *&globobs, *&expected );
    double setup_trigger[5] = { 0.005, 0, 1, 1, 1 };
    NPmaker( "trigger", setup_trigger, *&nuispara, *&constraints, *&globobs, *&expected );
    double setup_isEM[5] = { 0.0526, 0, 1, 1, 1 };
    NPmaker( "isEM", setup_isEM, *&nuispara, *&constraints, *&globobs, *&expected );
    double setup_iso[5] = { 0.004, 0, 1, 1, 1 };
    NPmaker( "iso", setup_iso, *&nuispara, *&constraints, *&globobs, *&expected );
    double setup_ESCALE[5] = { 0.003, 0, 1, 1, 1 };
    NPmaker( "ESCALE", setup_ESCALE, *&nuispara, *&constraints, *&globobs, *&expected );
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Migration
  if( !m_nomig )
  {
    int bin_index = CatNameToBinNumber( channelname );
    int number_SS_sources = ss_tool->GetNumberOfSources(energy);
    // loop over ss sources.
    for( int i_s = 0; i_s < number_SS_sources; i_s++ )
    {
      TString current_SS_source_name = ss_tool->GetNameOfSource( i_s, energy );
      TString ss_np_name = Form("shape_%s",current_SS_source_name.Data());
      double current_ss_value = ss_tool->GetValue( current_SS_source_name, bin_index, energy );
      int current_ss_sign = ss_tool->GetSign( current_SS_source_name, bin_index, energy );
      
      // Asymmetric migration uncertainties:
      double setup_ss_current[5] = { current_ss_value, 0, current_ss_sign, 1, 1 };
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
      if( current_SS_source_name.Contains(MCSampleName[0]) )
      	NPmaker( ss_np_name, setup_ss_current, nuispara, constraints, globobs, expected_spin0p );
     
      // spin2 shape systematics:
      else if( current_SS_source_name.Contains(Spin2Type) )
	NPmaker( ss_np_name, setup_ss_current, nuispara, constraints, globobs, expected_spin2p );
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Spurious signal
  if( !m_nobgm )
  {
    double ss_events = spurious_signal( channelname );
    double setup_bias[5] = { ss_events, -999, 1, 0, 1 }; // use Gaussian constraint term
    NPmaker( "bias", setup_bias, *&nuispara_uncorrelated, *&constraints_bias, *&globobs, *&expected_bias );
  }
  else w->factory("atlas_expected_bias[0]");
  
  //--------------------------------------//
  // SYSTEMATICS: Resolution:
  vector<TString> res_names; res_names.clear();
  if( !m_nores )
  {
    // an energy-dependent index, 0 inclusive, 1-11 categories:
    int bin_index = CatNameToBinNumber( channelname );
    double setup_AllRes[5] = { 0.0, 0, 1, 1, 1 };
    
    // loop overs Res sources:
    int number_Res_sources = res_tool->GetNumberOfSources(energy);
    for( int i_s = 0; i_s < number_Res_sources; i_s++ )
    {
      TString current_Res_source_name = res_tool->GetNameOfSource( i_s, energy );
      TString full_res_name = Form("EM_%s_%iTeV",current_Res_source_name.Data(),energy);
      res_names.push_back(full_res_name);
      double current_sys_value = res_tool->GetValue( current_Res_source_name, bin_index, energy );
      int current_sys_sign = res_tool->GetSign( current_Res_source_name, bin_index, energy );
      setup_AllRes[0] = current_sys_value;
      setup_AllRes[2] = current_sys_sign;
      
      // Create NP for _inc as well as the spin0 and current spin2 processes:
      shapeNPmaker( full_res_name, "_inc", setup_AllRes, *&nuispara, *&constraints, *&globobs, *&expected_shape );
      shapeNPmaker( full_res_name, Form("_%s",MCSampleName[0].Data()), setup_AllRes, *&nuispara, *&constraints, *&globobs, *&expected_shape );
      shapeNPmaker( full_res_name, Form("_%s",Spin2Type.Data()), setup_AllRes, *&nuispara, *&constraints, *&globobs, *&expected_shape );
    }
  }
  
  //--------------------------------------//
  // SYSTEMATICS: Energy-scale
  vector<TString> ess_names; ess_names.clear();
  if( !m_noess )
  {
    // an energy-dependent index, 0 inclusive, 1-11 categories:
    int bin_index = CatNameToBinNumber( channelname );
    double setup_AllESS[5] = { 0.0, 0, 1, 1, 1 };
    // loop overs ESS sources:
    int number_ESS_sources = ess_tool->GetNumberOfSources(energy);
    for( int i_s = 0; i_s < number_ESS_sources; i_s++ )
    {
      TString current_ESS_source_name = ess_tool->GetNameOfSource( i_s, energy );
      TString full_ess_name = Form("EM_%s_%iTeV",current_ESS_source_name.Data(),energy);
      ess_names.push_back(full_ess_name);
      double current_sys_value = ess_tool->GetValue( current_ESS_source_name, bin_index, energy );
      int current_sys_sign = ess_tool->GetSign( current_ESS_source_name, bin_index, energy );
      setup_AllESS[0] = current_sys_value;
      setup_AllESS[2] = current_sys_sign;
      NPmaker( full_ess_name, setup_AllESS, *&nuispara, *&constraints, *&globobs, *&expected_shape );
    }
  }
  
  //--------------------------------------//
  // Parameters of interest (POIs):
  // 7 TeV:
  RooRealVar *mu_7TeV = new RooRealVar("mu_7TeV","mu_7TeV",1,-100,100);
  RooRealVar *mu_BR_gg_7TeV = new RooRealVar("mu_BR_gg_7TeV","mu_BR_gg_7TeV",1,-100,100);
  RooRealVar *mu_ggF_7TeV = new RooRealVar("mu_ggF_7TeV","mu_ggF_7TeV",1,-100,100);
  RooRealVar *mu_VBF_7TeV = new RooRealVar("mu_VBF_7TeV","mu_VBF_7TeV",1,-100,100);
  RooRealVar *mu_WH_7TeV = new RooRealVar("mu_WH_7TeV","mu_WH_7TeV",1,-100,100);
  RooRealVar *mu_ZH_7TeV = new RooRealVar("mu_ZH_7TeV","mu_ZH_7TeV",1,-100,100);
  RooRealVar *mu_VH_7TeV = new RooRealVar("mu_VH_7TeV","mu_VH_7TeV",1,-100,100);
  RooRealVar *mu_tH_7TeV = new RooRealVar("mu_tH_7TeV","mu_tH_7TeV",1,-100,100);
  RooRealVar *mu_ttH_7TeV = new RooRealVar("mu_ttH_7TeV","mu_ttH_7TeV",1,-100,100);
  RooRealVar *mu_VBFVH_7TeV = new RooRealVar("mu_VBFVH_7TeV","mu_VBFVH_7TeV",1,-100,100);
  // 8 TeV:
  RooRealVar *mu_8TeV = new RooRealVar("mu_8TeV","mu_8TeV",1,-100,100);
  RooRealVar *mu_BR_gg_8TeV = new RooRealVar("mu_BR_gg_8TeV","mu_BR_gg_8TeV",1,-100,100);
  RooRealVar *mu_ggF_8TeV = new RooRealVar("mu_ggF_8TeV","mu_ggF_8TeV",1,-100,100);
  RooRealVar *mu_VBF_8TeV = new RooRealVar("mu_VBF_8TeV","mu_VBF_8TeV",1,-100,100);
  RooRealVar *mu_WH_8TeV = new RooRealVar("mu_WH_8TeV","mu_WH_8TeV",1,-100,100);
  RooRealVar *mu_ZH_8TeV = new RooRealVar("mu_ZH_8TeV","mu_ZH_8TeV",1,-100,100);
  RooRealVar *mu_VH_8TeV = new RooRealVar("mu_VH_8TeV","mu_VH_8TeV",1,-100,100);
  RooRealVar *mu_tH_8TeV = new RooRealVar("mu_tH_8TeV","mu_tH_8TeV",1,-100,100);
  RooRealVar *mu_ttH_8TeV = new RooRealVar("mu_ttH_8TeV","mu_ttH_8TeV",1,-100,100);
  RooRealVar *mu_VBFVH_8TeV = new RooRealVar("mu_VBFVH_8TeV","mu_VBFVH_8TeV",1,-100,100);
  if( channelname.Contains("7TeV") )
  {
    expected->add( RooArgSet( *mu_7TeV, *mu_BR_gg_7TeV ) );
    expected_proc_ggF->add( RooArgSet( *mu_ggF_7TeV, *mu_tH_7TeV) );
    expected_proc_VBF->add( RooArgSet( *mu_VBF_7TeV, *mu_VBFVH_7TeV ) );
    expected_proc_WH->add( RooArgSet( *mu_WH_7TeV, *mu_VH_7TeV, *mu_VBFVH_7TeV ) );
    expected_proc_ZH->add( RooArgSet( *mu_ZH_7TeV, *mu_VH_7TeV, *mu_VBFVH_7TeV ) );
    expected_proc_ttH->add( RooArgSet( *mu_ttH_7TeV, *mu_tH_7TeV ) );
  }
  else if( channelname.Contains("8TeV") )
  {
    expected->add( RooArgSet( *mu_8TeV, *mu_BR_gg_8TeV ) );
    expected_proc_ggF->add( RooArgSet( *mu_ggF_8TeV, *mu_tH_8TeV) );
    expected_proc_VBF->add( RooArgSet( *mu_VBF_8TeV, *mu_VBFVH_8TeV ) );
    expected_proc_WH->add( RooArgSet( *mu_WH_8TeV, *mu_VH_8TeV, *mu_VBFVH_8TeV ) );
    expected_proc_ZH->add( RooArgSet( *mu_ZH_8TeV, *mu_VH_8TeV, *mu_VBFVH_8TeV ) );
    expected_proc_ttH->add( RooArgSet( *mu_ttH_8TeV, *mu_tH_8TeV ) );
  }
  
  // =========================== End of declaration of POIs =================== //
  
  RooProduct expectation_proc_VHttH("expectation_proc_VHttH","expectation_proc_VHttH",*expected_proc_VHttH);
  RooProduct expectation_common("expectation_common","expectation_common",*expected);
  RooProduct expectation_spin0p("expectation_spin0p","expectation_spin0p",*expected_spin0p);
  RooProduct expectation_spin2p("expectation_spin2p","expectation_spin2p",*expected_spin2p);
  RooProduct expectation_proc_ggF("expectation_proc_ggF","expectation_proc_ggF",*expected_proc_ggF);
  RooProduct expectation_proc_VBF("expectation_proc_VBF","expectation_proc_VBF",*expected_proc_VBF);
  RooProduct expectation_proc_WH("expectation_proc_WH","expectation_proc_WH",*expected_proc_WH);
  RooProduct expectation_proc_ZH("expectation_proc_ZH","expectation_proc_ZH",*expected_proc_ZH);
  RooProduct expectation_proc_ttH("expectation_proc_ttH","expectation_proc_ttH",*expected_proc_ttH);
  
  // Spurious signal term will assume the shape of "inclusive" pdf.
  w->import(expectation_proc_VHttH);
  w->import(expectation_common);
  w->import(expectation_spin0p);
  w->import(expectation_spin2p);
  
  w->import(expectation_proc_ggF);
  w->import(expectation_proc_VBF);
  w->import(expectation_proc_WH);
  w->import(expectation_proc_ZH);
  w->import(expectation_proc_ttH);
  
  w->import(*expected_shape);
  w->import(*expected_bias);
  
  // Declare the observables:
  w->factory(Form("atlas_invMass[%f,%f]",xmin,xmax));
  w->defineSet("Observables","atlas_invMass");
  
  // Construct the signal and background PDFs:
  w->factory("epsilon[0.5,0,1]");
  w->factory("sum::epsilon_min_1(plusone[1.], prod::mineps(minusone[-1.],epsilon))");
  signalPdfBuilder( *&w, value_shape, ess_names, res_names, "_inc"    ); // inclusive signal pdf
  signalPdfBuilder( *&w, value_shape, ess_names, res_names, "_spin0p" ); // inclusive signal pdf
  signalPdfBuilder( *&w, value_shape, ess_names, res_names, "_spin2p" ); // ggF signal pdf
  backgroundPdfBuilder( *&w, *&nuispara_bkg, channelname);
  
  nuispara_uncorrelated->add(*nuispara_bkg);
  
  // build the signal normalization
  TString nSM_0p = Form("%f", value_0p[1]);
  TString nSM_2p = Form("%f", value_2p[1]);
  cout << "  nSM_0p for " << channelname << " = " << nSM_0p << endl;
  cout << "  nSM_2p for " << channelname << " = " << nSM_2p << endl;
  
  
  // Normalization for each process follows such pattern:
  // mu*isEM*lumi*migr => expectation_common
  w->factory((TString)"prod::atlas_nsig_spin0p(atlas_nSM_spin0p["+nSM_0p+(TString)"],expectation_common,expectation_spin0p,epsilon)");
  w->factory((TString)"prod::atlas_nsig_spin2p(atlas_nSM_spin2p["+nSM_2p+(TString)"],expectation_common,expectation_spin2p,epsilon_min_1)");
  w->factory("SUM::modelSB(atlas_nsig_spin0p*signalPdf_spin0p,atlas_nsig_spin2p*signalPdf_spin2p,atlas_expected_bias*signalPdf_inc,atlas_nbkg*bkgPdf)");
  w->Print();
  
    
  if( channelname == leadingchannel )
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
  // all uncorrelated nuisance parameters should be included in nuispara_uncorrelated.
  TString correlated;
  
  if( m_decorr_mu ) correlated = "epsilon";
  else
  {
    if( channelname.Contains("8TeV") ) correlated = "mu_8TeV,epsilon,mu_BR_gg_8TeV";
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
  cout << " For channel " << channelname << " the following variables will not be renamed : " << correlated << endl;
  
  // sub-channel labeling
  // import the workspace w to another workspace and add channelname as a suffix to all nodes and variables of w.
  // the correlated nuisance parameters and their respective global observables will not be renamed.
  RooWorkspace* wchannel = new RooWorkspace("wchannel"+channelname);
  wchannel->import( (*w->pdf("model")), RenameAllNodes(channelname), RenameAllVariablesExcept(channelname,correlated), Silence() );
  
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
  //   From nuispara_uncorrelated:
  cout << " Adding uncorrelated nuisance parameters to nuisanceParameters RooArgSet" << endl;
  TIterator *iter_nui_uncorrelated = nuispara_uncorrelated->createIterator();
  RooRealVar* parg_nui_uncorrelated = NULL;
  while( (parg_nui_uncorrelated = (RooRealVar*)iter_nui_uncorrelated->Next()) )
  {
    TString name_of_nuisance = parg_nui_uncorrelated->GetName()+(TString)"_"+channelname;
    nuisance_wchannel->add( *(RooRealVar*)wchannel->obj(name_of_nuisance) );
  }
  
  // The following are nps from background pdf, which don't have constraints: 
  //   From nuispara_bkg:
  RooArgSet* nuispara_bkg_wchannel = new RooArgSet();
  TIterator *iter_nui_bkg = nuispara_bkg->createIterator();
  RooRealVar* parg_nui_bkg = NULL;
  while( (parg_nui_bkg = (RooRealVar*)iter_nui_bkg->Next()) )
  {
    TString name_of_parameter = parg_nui_bkg->GetName()+(TString)"_"+channelname;
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
  while( (parg_global = (RooRealVar*)iter_global->Next()) )//&& (channelname==channel_constraints_attached) )
  {
    TString name_of_global = parg_global->GetName()+(TString)"_"+channelname;
    cout << " Channel Name " << channelname << " getting global observable " << parg_global->GetName() << endl;
    
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
  TIterator *iter_observable = w->set("Observables")->createIterator();
  RooRealVar *parg_observable;
  while( (parg_observable = (RooRealVar*)iter_observable->Next()) )
  {
    TString name_of_observable = parg_observable->GetName()+(TString)"_"+channelname;
    if( (bool)wchannel->obj(name_of_observable) == true )
      observable_wchannel->add( *(RooRealVar*)wchannel->obj(name_of_observable) );
    else
      observable_wchannel->add( *(RooRealVar*)wchannel->obj(parg_observable->GetName()) );
  }
  
  RooArgSet* muconstants_wchannel = new RooArgSet();
  if( channelname.Contains("7TeV") )
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
  else if( channelname.Contains("8TeV") )
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
  wchannel->defineSet("Observables", *observable_wchannel);
  wchannel->defineSet("nuisanceParameters", *nuisance_wchannel);
  wchannel->defineSet("globalObservables", *global_wchannel);
  
  //--------------------------------------//
  // Import data set and set up the background related nuisance parameter values
  dataInputDir = Form("%s/%s/mass_points_%iTeV/",master_output.Data(),jobname.Data(),energy);
  RooDataSet *obsdata = RooDataSet::read(dataInputDir+(TString)"mass_"+channelname+(TString)".txt",RooArgList(*wchannel->var("atlas_invMass_"+channelname)));
  obsdata->SetNameTitle("obsdata","obsdata");
  (*wchannel->var("atlas_nbkg_"+channelname) ).setVal(obsdata->numEntries() );
  (*wchannel->pdf("bkgPdf_"+channelname)).fitTo( *obsdata, Minos(RooArgSet(*nuispara_bkg_wchannel ) ) );
  (*wchannel->var("atlas_nbkg_"+channelname) ).setVal(obsdata->numEntries() );
  nuispara_bkg_wchannel->Print("v");
  wchannel->import(*obsdata);
  
  //--------------------------------------//
  // Create a binned data set:
  RooRealVar wt("wt","wt",1);
  RooArgSet* obs_plus_wt = new RooArgSet();
  obs_plus_wt->add(wt);
  obs_plus_wt->add(*wchannel->var("atlas_invMass_"+channelname));
  
  // histogram to store binned data:
  TH1F* h_data = new TH1F("h_data", "", 240, xmin, xmax );
  RooArgSet* obs = (RooArgSet*)obsdata->get();
  RooRealVar* xdata = (RooRealVar*)obs->find("atlas_invMass_"+channelname);
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
    wchannel->var("atlas_invMass_"+channelname)->setVal( mass_val );
    double weight = h_data->GetBinContent(ibin);
    wt.setVal(weight);
    obsdatabinned->add( RooArgSet(*wchannel->var("atlas_invMass_"+channelname), wt), weight );
    mass_val += 0.25;
  }
  wchannel->import(*obsdatabinned);
  
  //--------------------------------------//
  // Create a binned Asimov dataset:
  // Create Asimov spin 0+ data:
  CreateAsimovData( channelname, wchannel, obsdata, wt, xmin, xmax, 1, option );
  // Create Asimov spin 2+ data:
  CreateAsimovData( channelname, wchannel, obsdata, wt, xmin, xmax, 0, option );
  
  //--------------------------------------//
  // Plot the single-channel fit:
  plotBackgroundOnlyFit( wchannel, channelname );
  
  delete h_data;
  return wchannel;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// readinput:

vector<double> DMWorkspace::readinput( TString channelname, TString signal_name )
{
  paramInputDir = Form("%s/%s/parameterization_%iTeV/FinalSignal/",master_output.Data(),jobName.Data(),CatNameToEnergy(channelname));
  TString SignalFileName = Form("%s/fitpars_%s_%s.txt",paramInputDir.Data(),channelname.Data(),signal_name.Data());
  cout << "Reading file " << SignalFileName.Data() << endl;
  ifstream file_to_read(SignalFileName.Data(),ios::in);
  assert(file_to_read);
  double value[9];
  while( !file_to_read.eof() )
    file_to_read >> value[0] >> value[1] >> value[2] >> value[3] >> value[4] >> value[5] >> value[6] >> value[7] >> value[8];
  file_to_read.close();
  
  // The signal yield we use now correspond to 1 fb-1. Scaling to current luminosity.
  if( channelname.Contains("7TeV") ) value[1] *= luminosity_7TeV;
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

void DMWorkspace::backgroundPdfBuilder( RooWorkspace *&w, RooArgSet *&nuispara, TString channelname )
{
  int cate = CatNameToIndex( channelname );
  
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

double DMWorkspace::spurious_signal( TString channelname )
{
  // first entry [0] is inclusive 7 and 8 TeV!
  double spurious[23] = { 1.6, 0.17, 0.12, 0.03, 0.03, 0.05, 0.15, 0.08, 0.20, 0.04, 0.03, 0.03, 0.13, 0.10, 0.02, 0.02, 0.04, 0.1, 0.06, 0.14, 0.02, 0.02, 0.02 };
  int cate = CatNameToIndex( channelname ); 
  double result = ( channelname.Contains("8TeV") ) ? spurious[cate]*luminosity_8TeV : spurious[cate]*luminosity_7TeV;
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

void DMWorkspace::CreateAsimovData( TString channelname, RooWorkspace* wchannel, RooDataSet *obsdata, RooRealVar wt, double xmin, double xmax, int epsilon, TString option )
{
  TString spin = ( epsilon == 1 ) ? "0p" : "2p";
  cout << "CreateAsimovData( " << spin << " )" << endl;
  
  int npoints_Asimov = 275;
  
  // This is the dataset to be returned:
  RooDataSet *AsimovData = new RooDataSet( Form("asimovdatabinned%s",spin.Data()), Form("asimovdatabinned%s",spin.Data()), RooArgSet(*wchannel->var("atlas_invMass_"+channelname),wt), WeightVar(wt) );
    
  // Load the PDF from the workspace:
  //wchannel->Print("v");
  RooAbsPdf *current_pdf = (RooAbsPdf*)(wchannel->pdf("modelSB_"+channelname));
  double initial_epsilon = (wchannel->var("epsilon"))->getVal();
  (wchannel->var("epsilon"))->setVal(epsilon);
  double initial_mu;
  if( option.Contains("decorrmu") )
  {
    initial_mu = (wchannel->var(Form("mu_%s",channelname.Data())))->getVal();
    (wchannel->var(Form("mu_%s",channelname.Data())))->setVal(1.0);
  }
  else
  {
    if( channelname.Contains("7TeV") )
    {
      initial_mu = (wchannel->var("mu_7TeV"))->getVal();
      (wchannel->var("mu_7TeV"))->setVal(1.0);
    }
    else if( channelname.Contains("8TeV") )
    {
      initial_mu = (wchannel->var("mu_8TeV"))->getVal();
      (wchannel->var("mu_8TeV"))->setVal(1.0);
    }
  }
  
  // use fit result or integrals and sidebands to get the estimate of the background:
  double total_BkgEvents = obsdata->sumEntries();
  double width = ( xmax - xmin ) / ((double)npoints_Asimov);
  
  // loop over the number of asimov points:
  double count_Asimov = 0.0;
  for( int i_p = 0; i_p < npoints_Asimov; i_p++ )
  {
    double mass_value = xmin + ( 0.5 * width ) + ( width * (double)i_p );
    (wchannel->var("atlas_invMass_"+channelname))->setRange("range_Integral", mass_value-(0.5*width), mass_value+(0.5*width));
    RooAbsReal *integral = (RooAbsReal*)current_pdf->createIntegral(RooArgSet(*wchannel->var("atlas_invMass_"+channelname)), NormSet(*wchannel->var("atlas_invMass_"+channelname)), Range("range_Integral"));
    double weight_value = total_BkgEvents * integral->getVal();
    count_Asimov += weight_value;
    (wchannel->var("atlas_invMass_"+channelname))->setVal(mass_value);
    wt.setVal(weight_value);
    AsimovData->add( RooArgSet( *wchannel->var("atlas_invMass_"+channelname), wt ), weight_value );
  }
  if( fabs((count_Asimov-obsdata->sumEntries())/count_Asimov) > 0.04 ){ cout << "Bad Asimov Data: D=" << obsdata->sumEntries() << " A=" << count_Asimov << endl; exit(0); }
  wchannel->import(*AsimovData);
  (wchannel->var("epsilon"))->setVal(initial_epsilon);

  if( option.Contains("decorrmu") ) (wchannel->var(Form("mu_%s",channelname.Data())))->setVal(initial_mu);
  else
  {
    if( channelname.Contains("7TeV") )(wchannel->var("mu_7TeV"))->setVal(initial_mu);
    else if( channelname.Contains("8TeV") )(wchannel->var("mu_8TeV"))->setVal(initial_mu);
  }
}

/**
   Plot the fits produced by the specified model.
   @param plotOptions - options for what fits to plot etc.
   @returns void
*/
void DMWorkspace::plotFit(TString plotOptions)
{
  cout << "plotBackgroundOnlyFit( " << channelname << " )" << endl;
  TCanvas *c = new TCanvas();
  RooPlot* frame =  (*wchannel->var("atlas_invMass_"+channelname)).frame(55);
  wchannel->data("obsdata")->plotOn(frame);
  (*wchannel->pdf("model_"+channelname)).plotOn(frame, LineColor(2));
  (*wchannel->pdf("model_"+channelname)).plotOn(frame,Components( (*wchannel->pdf("bkgPdf_"+channelname)) ) , LineColor(4));
  double chi2 = frame->chiSquare() ;
  frame->SetYTitle("Events / GeV");
  frame->SetXTitle("M_{#gamma#gamma} [GeV]");
  frame->Draw();
  
  TLatex lresult3;
  lresult3.SetNDC();
  lresult3.SetTextColor(1);
  lresult3.DrawLatex(0.5,0.78, channelname);
  
  system(Form("mkdir -vp %s/figures/",outputDir.Data()));
  PrintCanvas(c, Form("%s/figures/data_fit_%s",outputDir.Data(),channelname.Data()));
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
