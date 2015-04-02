//////////////////////////////////////////
//                                      //
//  spinv1_makespace.cc                 //
//                                      //
//  Authors: Andrew Hard                //
//           Haichen Wang               //
//           Hongtao Yang               //
//                                      //
//  Date: 09/04/2014                    //
//                                      //
//  This program generates the 10x1D    //
//  spin analysis workspace from inputs //
//  for the 7 and 8 TeV spin analyses.  //
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

#include "spinv1_makespace.hh"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GenerateSingleChannel:

RooWorkspace* GenerateSingleChannel( TString channelname, TString leadingchannel = "costs1_7TeV", TString option = "" );

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// main:

int main( int argc, char **argv )
{
  // Check all inputs provided:
  if( argc < 4 ){ cout<<"Usage: ./<exe> <jobname> <options> <Spin2Type>"<<endl; return 0; }
  SetAtlasStyle();
  jobname = argv[1];
  TString option = argv[2];
  option.ToLower();
  Spin2Type = argv[3];
  
  output_directory = option.Contains("decorrmu") ? Form("%s/%s/ws_decorrmu_%s",master_output.Data(),jobname.Data(),Spin2Type.Data()) : Form("%s/%s/ws_%s",master_output.Data(),jobname.Data(),Spin2Type.Data());
  
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

RooWorkspace* GenerateSingleChannel( TString channelname, TString leadingchannel, TString option )
{
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
