//////////////////////////////////////////
//                                      //
//  spinv1_makespace.hh                 //
//                                      //
//  Authors: Andrew Hard                //
//           Haichen Wang               //
//           Hongtao Yang               //
//                                      //
//  Date: 09/04/2014                    //
//                                      //
//  Used to produce official workspace  //
//  for Higgs to gamma gamma 2013 spin  //
//  analysis publication.               //
//                                      //
//  v1 implements the pT category       //
//                                      //
//////////////////////////////////////////

#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "CommonFunc.h"
#include "statistics.hh"

// Systematic Uncertainty readers:
#include "ESSReader.hh"
#include "ResReader.hh"
#include "SigShapeReader.hh"

// GLOBAL SETTINGS DEFINED HERE:
#include "spinv1_Master.hh"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;

TString jobname;
TString Spin2Type;

int nchannels = 22;

// locations of the input files are determined by the master directory (defined in spin_Master.hh) and in the code.
TString dataInputDir = "./";
TString paramInputDir = "./";
TString xsecInputDir = "./";
TString output_directory;
RooArgSet* empty = new RooArgSet();

ESSReader *ess_tool;
ResReader* res_tool;
SigShapeReader* ss_tool;
  
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CatNameToEnergy:

int CatNameToEnergy( TString channelname )
{
  int result = 8;
  if( channelname.Contains("7TeV") ) result = 7;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CatNameToIndex:

int CatNameToIndex( TString channelname )
{
  int category_index = -1;
  if( channelname == "inclusive" ) category_index = 0;
  else
  {
    for( int i_c = 1; i_c <= nchannels; i_c++ )
    {
      if( channelname.Contains(Form("costs%i_7TeV",i_c)) ) category_index = i_c;
      else if( channelname.Contains(Form("costs%i_8TeV",i_c)) ) category_index = i_c+11;
    }
  }
  return category_index;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CatNameToBinNumber: (modulo sqrt(s))

int CatNameToBinNumber( TString channelname )
{
  int category_bin = -1;
  if( channelname == "inclusive" ) category_bin = 0;
  else
  {
    for( int i_b = 1; i_b <= nchannels; i_b++ )
    {
      if( channelname.Contains(Form("costs%i_7TeV",i_b)) || channelname.Contains(Form("costs%i_8TeV",i_b)) ) category_bin = i_b;
    }
  }
  return category_bin;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// backgroundPdfBuilder:

void backgroundPdfBuilder( RooWorkspace *&w, RooArgSet *&nuispara, TString channelname )
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

double spurious_signal( TString channelname )
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

void NPmaker( const char* varname, double setup[5], RooArgSet *&nuispara, RooArgSet *&constraints, RooArgSet *&globobs, RooArgSet *&expected )
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
////////// readinput:

vector<double> readinput( TString channelname, TString signal_name )
{
  paramInputDir = Form("%s/%s/parameterization_%iTeV/FinalSignal/",master_output.Data(),jobname.Data(),CatNameToEnergy(channelname));
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
////////// shapeNPmaker:

void shapeNPmaker( const char* varnameNP, const char* proc, double setup[5], RooArgSet *&nuispara, RooArgSet *&constraints, RooArgSet *&globobs, RooArgSet *&expected )
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
////////// signalPdfBuilder:

void signalPdfBuilder( RooWorkspace *&w, vector<double> value, vector<TString> ess_parnames, vector<TString> res_parnames, TString procname )
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
////////// plotBackgroundOnlyFit:

void plotBackgroundOnlyFit( RooWorkspace* wchannel, TString channelname )
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
  
  system(Form("mkdir -vp %s/figures/",output_directory.Data()));
  PrintCanvas(c, Form("%s/figures/data_fit_%s",output_directory.Data(),channelname.Data()));
  delete c;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CreateAsimovData:

void CreateAsimovData( TString channelname, RooWorkspace* wchannel, RooDataSet *obsdata, RooRealVar wt, double xmin, double xmax, int epsilon, TString option )
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PlotNuisParams:

void PlotNuisParams( RooArgSet nuis, TString signal_type )
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
  can->Print( Form("%s/figures/nuisparams_%s.eps",output_directory.Data(),signal_type.Data()) );
  can->Print( Form("%s/figures/nuisparams_%s.png",output_directory.Data(),signal_type.Data()) );
  can->Print( Form("%s/figures/nuisparams_%s.C",output_directory.Data(),signal_type.Data()) );
  delete can;
}
