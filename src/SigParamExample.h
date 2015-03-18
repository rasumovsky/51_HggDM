/////////////////////////////////////////
//                                     //
//  spinv1_parameterization.hh         //
//                                     //
//  Author: Andrew Hard                //
//  Date: 21/01/2014                   //
//                                     //
//  Used by:                           //
//    spinv1_parameterization.cc       //
//                                     //
/////////////////////////////////////////

#include "../inc/CommonHead.h"
#include "../inc/CommonFunc.h"
#include "../inc/RooFitHead.h"
#include "../inc/statistics.hh"

#include "../inc/ntup_SpinPub.hh"
#include "../src/spinv1_Master.hh"
#include "../src/SigShapeReader.hh"

using namespace std;
using namespace RooFit;
using namespace CommonFunc;

int energy;

ntup_SpinPub* p;

TString output_directory;

///////////////////////////////////////////////////////////////////////////////
// For the parameterizaion:

double lumi = 1000.0;
const int nmasspoint = 11;
const int nprocess = 5;
int signalmass[nmasspoint] = {100,105,110,115,120,125,130,135,140,145,150};

TString jobname;
TString sample_list;

// get the number of categories:
int ncat;

vector<int> read_runs;
vector<double> read_SumWt;

TChain* chain;
int numev;

TTree* tree[50][nmasspoint];
double mass[50] = {0.0};
double weight[50] = {0.0};

double all_w[50][nmasspoint] = {{0}};
double mode_w[50][nprocess][nmasspoint] = {{{0}}};
double counter[50][nprocess] = {{0.0}};

///////////////////////////////////////////////////////////////////////////////
// For the signal yields:

// validation histograms:
TH1F* h_costs_interf[NumMCSamples];
TH1F* h_costs_nointerf[NumMCSamples];
TH1F* h_costs_VarHi[NumMCSamples];
TH1F* h_costs_VarLo[NumMCSamples];

TCanvas* can;

// yield information for categories and spin types:
double total_yield[NumMCSamples] = {0.0};
double binned_yield[NumMCSamples][50] = {{0.0}};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SpinSampleType:

TString SpinSampleType()
{
  TString result = "";
  // get generator:
  if( p->Run == 189380 ) result = "MG5_Spin2_kg1_kq1";
  else if( p->Run == 189381 ) result = "MG5_Spin2_kg1_kq0";
  else if( p->Run == 189382 ) result = "MG5_Spin2_kg05_kq1";
  else if( p->Run == 181762 ) result = "AMC_Spin2_kg1_kq1";
  else result = "POW_Spin0";
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// get_event:

void get_event( int i ) 
{
  if ( p->LoadTree(i) < 0 ) 
  { 
    cout << "\nProblem in LoadTree." << "\nEntry: " << i << endl;
    exit(0);
  }
  p->fChain->GetEntry(i);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PrintProgressBar:

void PrintProgressBar( int index, int total )
{
  if( index%10000 == 0 )
  {
    TString print_bar = " [";
    for( int bar = 0; bar < 20; bar++ )
    {
      double current_fraction = double(bar) / 20.0;
      if( double(index)/double(total) > current_fraction ) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.*(double(index)/double(total));
    TString text = Form("%s %2.2f ",print_bar.Data(),percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// OutputProcess:

void OutputProcess( int cat_number, int mass, double ggf_w, double vbf_w, double tth_w, double wh_w, double zh_w, double mu_CB, double sigma_CB, double alpha_CB, double n_CB, double mu_GA, double sigma_GA, double frac )
{
  ofstream outfile;
  double truthmass = double(mass)/10.;
  outfile.open(Form("%s/ggf/ggf_cate_fit_%d_%d.txt",output_directory.Data(),mass,cat_number),ios::out);
  outfile<<truthmass<<"\t"<<ggf_w<<"\t"<<mu_CB<<"\t"<<sigma_CB<<"\t"<<alpha_CB<<"\t"<<n_CB<<"\t"<<mu_GA<<"\t"<<sigma_GA<<"\t"<<frac<<endl;
  outfile.close();
  
  outfile.open(Form("%s/vbf/vbf_cate_fit_%d_%d.txt",output_directory.Data(),mass,cat_number),ios::out);
  outfile<<truthmass<<"\t"<<vbf_w<<"\t"<<mu_CB<<"\t"<<sigma_CB<<"\t"<<alpha_CB<<"\t"<<n_CB<<"\t"<<mu_GA<<"\t"<<sigma_GA<<"\t"<<frac<<endl;
  outfile.close();
  
  outfile.open(Form("%s/tth/tth_cate_fit_%d_%d.txt",output_directory.Data(),mass,cat_number),ios::out);
  outfile<<truthmass<<"\t"<<tth_w<<"\t"<<mu_CB<<"\t"<<sigma_CB<<"\t"<<alpha_CB<<"\t"<<n_CB<<"\t"<<mu_GA<<"\t"<<sigma_GA<<"\t"<<frac<<endl;
  outfile.close();
  
  outfile.open(Form("%s/wh/wh_cate_fit_%d_%d.txt",output_directory.Data(),mass,cat_number),ios::out);
  outfile<<truthmass<<"\t"<<wh_w<<"\t"<<mu_CB<<"\t"<<sigma_CB<<"\t"<<alpha_CB<<"\t"<<n_CB<<"\t"<<mu_GA<<"\t"<<sigma_GA<<"\t"<<frac<<endl;
  outfile.close();
  
  outfile.open(Form("%s/zh/zh_cate_fit_%d_%d.txt",output_directory.Data(),mass,cat_number),ios::out);
  outfile<<truthmass<<"\t"<<zh_w<<"\t"<<mu_CB<<"\t"<<sigma_CB<<"\t"<<alpha_CB<<"\t"<<n_CB<<"\t"<<mu_GA<<"\t"<<sigma_GA<<"\t"<<frac<<endl;
  outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PrintMassCate:

void PrintMassCate( Int_t mass, Int_t cat, Color_t color = kBlack )
{
  TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextColor(color);
  l.DrawLatex( 0.2, 0.85, Form("m_{#gamma#gamma}=%dGeV Category%d", mass, cat) );
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetWeight:

float GetWeight( double lumi, bool useInterference, bool useSetMass )
{
  double weight_result = 1.0;
  if( !(p->TYPE_DATA) )
  {
    double XS;
    if( useSetMass )
    {
      // updated 15/10/2013 from YR@8TeV for mH=126.5GeV:
      //double XSvals_8TeV[5] = { 18.82, 1.558, 0.1247, 0.6767, 0.4000 };
      //double XSvals_7TeV[5] = { 14.77, 1.206, 0.08326, 0.5555, 0.3227 };
      
      // updated 08/04/2014 from YR@8TeV for mH=125.6GeV:
      double XSvals_8TeV[5] = { 19.09, 1.572, 0.1274, 0.6931, 0.4091 };
      double XSvals_7TeV[5] = { 14.99, 1.214, 0.08508, 0.5688, 0.3299 };
            
      if( energy == 7 ) XS = XSvals_7TeV[p->MCtype];
      else if( energy == 8 ) XS = XSvals_8TeV[p->MCtype];
      else{ cout << "Undefined cross-section. Exiting." << endl; exit(0); }
    }
    else XS = p->xsection;
    double BR = p->branching_ratio;
    
    weight_result = lumi * p->MCweight * p->PUweight * p->ZVTXweight * p->pT_weight * XS * BR / p->NEVENT_weighted;
    if( useInterference ) weight_result *= p->interference_weight_CTS;
  }
  
  //cout << "Run=" << p->Run << " Event=" << p->Event << " MCweight=" << p->MCweight << " PUweight=" <<  p->PUweight << " ZVTXweight=" << p->ZVTXweight << " pT_weight=" << p->pT_weight << " NEVENT_weighted=" << p->NEVENT_weighted;
  return weight_result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetCTSPTtBin:

int GetCTSPTtBin( double cts, double pt )
{
  int result = -1;
  int bmax = 0;
  for( int i_b = 0; i_b < 10; i_b++ )
  {
    if( i_b > bmax ) bmax = i_b;
    double bin_lo = ((double)i_b) / 10;
    double bin_hi = ((double)i_b+1.0) / 10;
    if( cts > bin_lo && cts < bin_hi ){ result = i_b; break; }
  }
  if( pt > 125000 && pt < 300000 ) result = bmax + 1;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// DoSignalParameterization:

void DoSignalParameterization()
{
  // Loop over events to get parameterization:
  cout << "We have "<< numev << " events to be processed." << endl;
  for( int index = 0; index < numev; index++ )
  {
    p->fChain->GetEntry(index);
    PrintProgressBar( index, numev );
        
    // Only use Standard Model for parameterization:
    TString sample = SpinSampleType();
    if( !sample.Contains("POW") ) continue;
    
    // Check pT (Relative cuts) and all other cuts:
    if( !(p->flag_pre && p->flag_PID && p->flag_iso && p->flag_RelPt_35_25) ) continue;
    
    int icat = GetCTSPTtBin( p->costhetastar_CS, p->PT ) + 1;
    
    mass[icat] = p->mass_PV_EM/1000.0;
    // always use luminosity of 1fb-1 for signal parameterization:
    weight[icat] = GetWeight( lumi, true, false );
    
    mass[0] = mass[icat];
    weight[0] = weight[icat];
    
    int MCtype = p->MCtype;
    for( int m = 0; m < nmasspoint; m++ )
    {
      int mass_used = signalmass[m];
      if( p->Higgs_truth_mass == mass_used )
      {
	tree[icat][m]->Fill();
	all_w[icat][m] += weight[icat];
	mode_w[icat][MCtype][m] += weight[icat];
	if( mass_used == 125 ) counter[icat][MCtype] += weight[icat];
	
	tree[0][m]->Fill();
	all_w[0][m] += weight[0];
	mode_w[0][MCtype][m] += weight[0];
	if( mass_used == 125 ) counter[0][MCtype] += weight[0];
      }
    }
  } 
  
  // loop over the categories:
  for( int icat = 0; icat <= ncat; icat++ )
  {
    // Observables:
    RooRealVar v_mass("mass", "mass", 80, 180);
    RooRealVar v_weight("weight", "weight", -1e20, 1e20);
    
    // RooDataSets
    RooDataSet *data[nmasspoint];
    for( int n = 0; n < nmasspoint; n++ )
    {
      int mass_used2 = signalmass[n];
      char nameRD[20];
      sprintf(nameRD,"data%i",mass_used2);
      data[n] = new RooDataSet(nameRD, "", RooArgSet( v_mass, v_weight ), WeightVar( v_weight ), Import(*tree[icat][n]) ) ;
    }
    
    //--------------------------------------//
    // Roo variables:
    RooRealVar alpha("alpha","alpha", 2., 0.5, 5.);
    RooRealVar n_CB("n_CB", "n_CB", 10.);
    RooRealVar k_GA("k_GA", "k_GA",2., 1., 5.);
    RooRealVar frac("frac", "frac", 0.98, 0.8, 1.);
  
    RooRealVar a_mu("a_mu","a_mu",-0.3, -0.5, 0.5);
    RooRealVar b_mu("b_mu","b_mu",-0.001, -0.01, 0.01);
    RooRealVar a_sigma("a_sigma","a_sigma",1.5, 0.5, 5);
    RooRealVar b_sigma("b_sigma","b_sigma",0.01, 0.001, 0.5);
    
    //--------------------------------------//
    // For each mass point:
    RooFormulaVar *mu[20], *sigma[20], *sGA[20];
    RooCBShape *CB[20];
    RooGaussian *GA[20];
    RooAddPdf *Signal[20];
    
    double x_points[nmasspoint];
    double camila[nmasspoint];
    for( int i = 0; i < nmasspoint; i++ )
    {
      x_points[i] = signalmass[i];
      camila[i] = signalmass[i] - 125;
    }
    
    for( int i = 0; i < nmasspoint; i++ )
    {
      mu[i] = new RooFormulaVar(Form("mu_%d",signalmass[i]), Form("%f+@0+@1*%f",x_points[i],camila[i]), RooArgList(a_mu, b_mu));
      sigma[i] = new RooFormulaVar(Form("sigma_%d",signalmass[i]), Form("@0+@1*%f",camila[i]), RooArgList(a_sigma, b_sigma));
      sGA[i] = new RooFormulaVar(Form("sGA_%d",signalmass[i]),"@0*@1", RooArgList(k_GA, *sigma[i]));
      CB[i] = new RooCBShape(Form("CB_%d",signalmass[i]), Form("CB_%d",signalmass[i]), v_mass, *mu[i], *sigma[i], alpha, n_CB);
      GA[i] = new RooGaussian(Form("GA_%d",signalmass[i]), Form("GA_%d",signalmass[i]), v_mass, *mu[i], *sGA[i]);
      Signal[i] = new RooAddPdf(Form("Signal_%d",signalmass[i]), Form("Signal_%d",signalmass[i]), *CB[i], *GA[i], frac);
    }
    
    RooCategory sample("sample", "sample");
    map<string,RooDataSet*> dataMap;
    
    for( int s = 0; s < nmasspoint; s++ )
    {
      int massval = signalmass[s];
      string massname = Form("m_%i",massval);
      sample.defineType(Form("m_%i",massval));
      dataMap[massname] = data[s];
    }
    
    RooDataSet combData("combData", "combData", RooArgSet(v_mass, v_weight), WeightVar(v_weight), Index(sample), Import(dataMap));
  
    RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);
    
    // add pdfs for each mass point to RooSimultaneous:
    for( int s = 0; s < nmasspoint; s++ ) simPdf.addPdf( *Signal[s], Form("m_%d",signalmass[s]) );
    
    combData.Print();
    simPdf.Print();
    
    statistics::setDefaultPrintLevel(0);
    RooNLLVar *nll=(RooNLLVar*)simPdf.createNLL(combData);
    statistics::minimize(nll);
    
    TF1* f_tot = new TF1("f_tot","pol3",100,150);
    TGraph* tot_graph = new TGraph(nmasspoint, x_points, all_w[icat]);
    tot_graph->Fit(f_tot);
  
    TF1* f_ggf = new TF1("f_ggf","pol3",100,150);
    TGraph* ggf_graph = new TGraph(nmasspoint, x_points, mode_w[icat][0]);
    ggf_graph->Fit(f_ggf);
  
    TF1* f_vbf = new TF1("f_vbf","pol3",100,150);
    TGraph* vbf_graph = new TGraph(nmasspoint, x_points, mode_w[icat][1]);
    vbf_graph->Fit(f_vbf);
  
    TF1* f_tth = new TF1("f_tth","pol3",100,150);
    TGraph* tth_graph = new TGraph(nmasspoint, x_points, mode_w[icat][2]);
    tth_graph->Fit(f_tth);
  
    TF1* f_wh = new TF1("f_wh","pol3",100,150);
    TGraph* wh_graph = new TGraph(nmasspoint, x_points, mode_w[icat][3]);
    wh_graph->Fit(f_wh);
    
    TF1* f_zh = new TF1("f_zh","pol3",100,150);
    TGraph* zh_graph = new TGraph(nmasspoint, x_points, mode_w[icat][4]);
    zh_graph->Fit(f_zh);
  
    cout << "The number of events in each tree are: \n";
    for( int o = 0; o < nmasspoint; o++ )
    {
      int massp = signalmass[o];
      string pout = Form("M=%iGeV", massp);
      cout << pout << tree[icat][o]->GetEntries() << endl;
    }
    
    int nbin = 50;
    TCanvas* c1 = new TCanvas("c1","c1");
    c1->cd();
    
    RooPlot* frame[nmasspoint];
    for( int f = 0; f < nmasspoint; f++ )
    {
      int mass_bin = signalmass[f];
      char cutname2[10];
      sprintf(cutname2,"Mass_%i",mass_bin);
      frame[f] = v_mass.frame(Bins(nbin), Title(cutname2),Range(mass_bin-10,mass_bin+10) );
      char cutname[30];
      sprintf(cutname,"sample==sample::m_%i",mass_bin);
      combData.plotOn(frame[f],Cut(cutname));
      char cutname3[10];
      sprintf(cutname3,"m_%i",mass_bin);
      simPdf.plotOn(frame[f],Slice(sample,cutname3),ProjWData(sample,combData)) ;
      //simPdf.paramOn(frame[f]);
      frame[f]->GetXaxis()->SetTitle("m_{#gamma#gamma}[GeV]");
      frame[f]->Draw();
      PrintMassCate( mass_bin, icat );
      TString outputCanvas=Form("%s/Plots/m_%i_%d",output_directory.Data(),mass_bin,icat);
      PrintCanvas(c1, outputCanvas);
    }
    
    ofstream outfile;
    for( int i = 0; i < nmasspoint; i++ )
    {
      outfile.open(Form("%s/all/cate_fit_%d_%d.txt",output_directory.Data(),signalmass[i]*10,icat),ios::out);
      outfile<<Form("%d\t",signalmass[i])<<all_w[icat][i]<<"\t"<<mu[i]->getVal()<<"\t"<<sigma[i]->getVal()<<"\t"<<alpha.getVal()<<"\t"<<n_CB.getVal()<<"\t"<<mu[i]->getVal()<<"\t"<<sGA[i]->getVal()<<"\t"<<frac.getVal()<<endl;
      outfile.close();
      OutputProcess( icat, signalmass[i]*10, mode_w[icat][0][i], mode_w[icat][1][i], mode_w[icat][2][i], mode_w[icat][3][i], mode_w[icat][4][i], mu[i]->getVal(), sigma[i]->getVal(), alpha.getVal(), n_CB.getVal(), mu[i]->getVal(), sGA[i]->getVal(), frac.getVal() );
    }
    
    outfile.open(Form("%s/param_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<a_mu.getVal()<<"\t"<<b_mu.getVal()<<"\t"<<a_sigma.getVal()<<"\t"<<b_sigma.getVal()<<"\t"<<alpha.getVal()<<"\t"<<n_CB.getVal()<<"\t"<<k_GA.getVal()<<"\t"<<frac.getVal()<<endl;
    outfile.close();
    
    outfile.open(Form("%s/all/yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_tot->GetParameter(0)<<"\t"<<f_tot->GetParameter(1)<<"\t"<<f_tot->GetParameter(2)<<"\t"<<f_tot->GetParameter(3)<<endl;
    outfile.close();
  
    outfile.open(Form("%s/ggf/ggf_yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_ggf->GetParameter(0)<<"\t"<<f_ggf->GetParameter(1)<<"\t"<<f_ggf->GetParameter(2)<<"\t"<<f_ggf->GetParameter(3)<<endl;
    outfile.close();
  
    outfile.open(Form("%s/vbf/vbf_yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_vbf->GetParameter(0)<<"\t"<<f_vbf->GetParameter(1)<<"\t"<<f_vbf->GetParameter(2)<<"\t"<<f_vbf->GetParameter(3)<<endl;
    outfile.close();
  
    outfile.open(Form("%s/tth/tth_yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_tth->GetParameter(0)<<"\t"<<f_tth->GetParameter(1)<<"\t"<<f_tth->GetParameter(2)<<"\t"<<f_tth->GetParameter(3)<<endl;
    outfile.close();
  
    outfile.open(Form("%s/wh/wh_yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_wh->GetParameter(0)<<"\t"<<f_wh->GetParameter(1)<<"\t"<<f_wh->GetParameter(2)<<"\t"<<f_wh->GetParameter(3)<<endl;
    outfile.close();
  
    outfile.open(Form("%s/zh/zh_yield_%d.txt",output_directory.Data(),icat),ios::out);
    outfile<<icat<<"\t"<<f_zh->GetParameter(0)<<"\t"<<f_zh->GetParameter(1)<<"\t"<<f_zh->GetParameter(2)<<"\t"<<f_zh->GetParameter(3)<<endl;
    outfile.close();
  }
  
  cout << " " << endl;
  cout << "Parameterization Complete. Printing yield information." << endl;
  cout << " " << endl;
  double sum[5] = {0};
  for( int c = 0; c <= ncat; c++ )
  {
    cout << "  Category: " << c << endl;
    for( int m = 0; m < nprocess; m++ )
    {
      cout << "    Process " << m << " " << counter[c][m] << endl;
      if( c != 0 ) sum[m] += counter[c][m];
    }
  }
  cout << " " << endl;
  cout << "Sum per process: " << endl;
  for( int m = 0; m < nprocess; m++ ) cout << "  Process " << m << " " << sum[m] << endl;
  
  cout << " " << endl;
  cout << "Finished DoSignalParameterization() process." << endl;
  cout << " " << endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetYield:

double GetYield( double mass, double p0, double p1, double p2, double p3 )
{
  double result = p0 + (p1*mass) + (p2*mass*mass) + (p3*mass*mass*mass);
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// main:

int InterpolateSignal( int cat_number )
{
  cout << "InterpolateSignal( " << cat_number << " ) called" << endl;
  
  //--------------------------------------//
  // Open files to get inclusive parameters:
  // parameterized from -25 to 25, not 100-150;
  ifstream infile;
  infile.open(Form("%s/param_%d.txt",output_directory.Data(),cat_number),ios::in);
  double null, a_mu, b_mu, a_sigma, b_sigma, alpha, n, k_ga, frac;
  infile >> null >> a_mu >> b_mu >> a_sigma >> b_sigma >> alpha >> n >> k_ga >> frac;
  infile.close();
  
  ifstream yield_file;
  yield_file.open(Form("%s/all/yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double p0, p1, p2, p3;
  yield_file >> null >> p0 >> p1 >> p2 >> p3;
  
  //--------------------------------------//
  // Inputs per process:
  ifstream ggf_yield_file;
  ggf_yield_file.open(Form("%s/ggf/ggf_yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double ggf_p0, ggf_p1, ggf_p2, ggf_p3;
  ggf_yield_file >> null >> ggf_p0 >> ggf_p1 >> ggf_p2 >> ggf_p3;

  ifstream vbf_yield_file;
  vbf_yield_file.open(Form("%s/vbf/vbf_yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double vbf_p0, vbf_p1, vbf_p2, vbf_p3;
  vbf_yield_file >> null >> vbf_p0 >> vbf_p1 >> vbf_p2 >> vbf_p3;
  
  ifstream tth_yield_file;
  tth_yield_file.open(Form("%s/tth/tth_yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double tth_p0, tth_p1, tth_p2, tth_p3;
  tth_yield_file >> null >> tth_p0 >> tth_p1 >> tth_p2 >> tth_p3;
  
  ifstream wh_yield_file;
  wh_yield_file.open(Form("%s/wh/wh_yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double wh_p0, wh_p1, wh_p2, wh_p3;
  wh_yield_file >> null >> wh_p0 >> wh_p1 >> wh_p2 >> wh_p3;
  
  ifstream zh_yield_file;
  zh_yield_file.open(Form("%s/zh/zh_yield_%d.txt",output_directory.Data(),cat_number),ios::in);
  double zh_p0, zh_p1, zh_p2, zh_p3;
  zh_yield_file >> null >> zh_p0 >> zh_p1 >> zh_p2 >> zh_p3;
  
  ofstream outfile;
  ofstream ggf_outfile;
  ofstream vbf_outfile;
  ofstream tth_outfile;
  ofstream wh_outfile;
  ofstream zh_outfile;
  
  //--------------------------------------//
  // Interpolate for masses in steps of 0.5:
  for( float mass = 100; mass < 150.5; mass += 0.5 )
  {
    int mass_print = (int)(10*mass);
    double mass_fit = mass - 125.0;
    double mu_print = mass + a_mu + (b_mu * mass_fit);
    double sigma_print = a_sigma + (b_sigma * mass_fit);
    if( mass != (int)mass )
    {
      outfile.open(Form("%s/all/cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      outfile<<mass<<"\t"<<GetYield(mass,p0,p1,p2,p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      outfile.close();
      
      ggf_outfile.open(Form("%s/ggf/ggf_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      ggf_outfile<<mass<<"\t"<<GetYield(mass,ggf_p0,ggf_p1,ggf_p2,ggf_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      ggf_outfile.close();
      
      vbf_outfile.open(Form("%s/vbf/vbf_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      vbf_outfile<<mass<<"\t"<<GetYield(mass,vbf_p0,vbf_p1,vbf_p2,vbf_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      vbf_outfile.close();
      
      tth_outfile.open(Form("%s/tth/tth_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      tth_outfile<<mass<<"\t"<<GetYield(mass,tth_p0,tth_p1,tth_p2,tth_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      tth_outfile.close();
      
      wh_outfile.open(Form("%s/wh/wh_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      wh_outfile<<mass<<"\t"<<GetYield(mass,wh_p0,wh_p1,wh_p2,wh_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      wh_outfile.close();
      
      zh_outfile.open(Form("%s/zh/zh_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      zh_outfile<<mass<<"\t"<<GetYield(mass,zh_p0,zh_p1,zh_p2,zh_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      zh_outfile.close();
    }
    else
    {
      outfile.open(Form("%s/all/cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      outfile<<mass<<"\t"<<GetYield(mass,p0,p1,p2,p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      outfile.close();
      
      ggf_outfile.open(Form("%s/ggf/ggf_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      ggf_outfile<<mass<<"\t"<<GetYield(mass,ggf_p0,ggf_p1,ggf_p2,ggf_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      ggf_outfile.close();
      
      vbf_outfile.open(Form("%s/vbf/vbf_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      vbf_outfile<<mass<<"\t"<<GetYield(mass,vbf_p0,vbf_p1,vbf_p2,vbf_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      vbf_outfile.close();
      
      tth_outfile.open(Form("%s/tth/tth_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      tth_outfile<<mass<<"\t"<<GetYield(mass,tth_p0,tth_p1,tth_p2,tth_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      tth_outfile.close();
      
      wh_outfile.open(Form("%s/wh/wh_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      wh_outfile<<mass<<"\t"<<GetYield(mass,wh_p0,wh_p1,wh_p2,wh_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      wh_outfile.close();
      
      zh_outfile.open(Form("%s/zh/zh_cate_fit_%i_%i.txt",output_directory.Data(),mass_print,cat_number),ios::out);
      zh_outfile<<mass<<"\t"<<GetYield(mass,zh_p0,zh_p1,zh_p2,zh_p3)<<"\t"<<mu_print<<"\t"<<sigma_print<<"\t"<<alpha<<"\t"<<n<<"\t"<<mu_print<<"\t"<<k_ga*(sigma_print)<<"\t"<<frac<<endl;
      zh_outfile.close();
    }
  }
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// plot_CosThetaStar:

void plot_CosThetaStar()
{
  can->Clear();
  int colors[6] = {2,3,4,6,9,8};
  TLegend leg(0.2, 0.72, 0.6, 0.92);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.05);
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    if( i_s == 0 )
    {
      h_costs_interf[i_s]->GetXaxis()->SetTitle("Cos#theta* CS");
      h_costs_interf[i_s]->SetLineColor(colors[i_s]);
      h_costs_interf[i_s]->SetMarkerColor(colors[i_s]);
      h_costs_interf[i_s]->Draw("LEP");
      leg.AddEntry(h_costs_interf[i_s],MCSampleName[i_s],"l");
    }
    else
    {
      h_costs_nointerf[i_s]->GetXaxis()->SetTitle("Cos#theta* CS");
      h_costs_nointerf[i_s]->SetLineColor(colors[i_s]);
      h_costs_nointerf[i_s]->SetMarkerColor(colors[i_s]);
      h_costs_nointerf[i_s]->Draw("LEPSAME");
      leg.AddEntry(h_costs_nointerf[i_s],MCSampleName[i_s],"l");
    }
  }
  leg.Draw("SAME");
  can->Print(Form("%s/Plots/costhetastar_comparison.eps",output_directory.Data()));
  can->Clear();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetBinYields:

void GetBinYields( int number_of_bins )
{
  // scale all amounts to nominal PowHeg JP=0+:
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    double scale_factor = total_yield[0] / total_yield[i_s];
    for( int i_b = 0; i_b < number_of_bins; i_b++ )
    {
      binned_yield[i_s][i_b] = scale_factor * binned_yield[i_s][i_b];
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PrintFinalSignalInputs:

void PrintFinalSignalInputs()
{
  // Loop over inclusive files at 125.5 GeV:
  for( int i_b = 0; i_b < ncat; i_b++ )
  {
    ifstream signal_fileIN;
    signal_fileIN.open(Form("%s/all/cate_fit_1255_%i.txt",output_directory.Data(),i_b+1));
    double f_m;
    double f_v[8] = {0.0};
    while( !signal_fileIN.eof() ) signal_fileIN >> f_m >> f_v[0] >> f_v[1] >> f_v[2] >> f_v[3] >> f_v[4] >> f_v[5] >> f_v[6] >> f_v[7];
    
    ofstream signal_fileOUT[NumMCSamples];
    for( int i_s = 0; i_s < NumMCSamples; i_s++ )
    {
      double yield_current = binned_yield[i_s][i_b];
      signal_fileOUT[i_s].open(Form("%s/FinalSignal/fitpars_costs%i_%iTeV_%s.txt",output_directory.Data(),i_b+1,energy,MCSampleName[i_s].Data()));
      signal_fileOUT[i_s] << f_m << " " << yield_current << " " << f_v[1] << " " << f_v[2] << " " << f_v[3] << " " << f_v[4] << " " << f_v[5] << " " << f_v[6] << " " << f_v[7] << endl;
      cout << MCSampleName[i_s] << " " << f_m << " " << yield_current << " " << f_v[1] << " " << f_v[2] << " " << f_v[3] << " " << f_v[4] << " " << f_v[5] << " " << f_v[6] << " " << f_v[7] << endl;
      signal_fileOUT[i_s].close();
    }
    signal_fileIN.close();
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// NotZero:

bool NotZero( double value )
{
  bool result = ( fabs(value) > 0.000001 ) ? true : false;
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetSampleIndex:
int GetSampleIndex( TString sample_name )
{
  int result = -1;
  for( int i_n = 0; i_n < NumMCSamples; i_n++ )
  {
    if( sample_name.Contains(MCSampleName[i_n]) )
    {
      result = i_n;
      break;
    }
  }
  if( result == -1 )
  {
    cout << "GetSampleIndex( " << sample_name << " ) : No match for sample_name." << endl;
    exit(0);
  }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetYieldsAtSpecificMass:

void GetYieldsAtSpecificMass()
{
  int i_cutflow_counter[NumMCSamples][10] = {{0.0}};
  
  SigShapeReader* ss_tool = new SigShapeReader( file_name_SS_values_7TeV, file_name_SS_values_8TeV, nchannel_switches_7TeV, nchannel_switches_8TeV );
  
  // Make sure values are null:
  for( int i_t = 0; i_t < NumMCSamples; i_t++ )
  {
    total_yield[i_t] = 0.0;
    for( int i_b = 0; i_b < 50; i_b++ ) binned_yield[i_b][i_t] = 0.0;
  }
    
  // format histograms:
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    h_costs_interf[i_s] = new TH1F(Form("h_costs_interf_%s",MCSampleName[i_s].Data()), Form("h_costs_interf_%s",MCSampleName[i_s].Data()), ncat, 0.0, (double)ncat/10 );
    h_costs_interf[i_s]->Sumw2();
    h_costs_nointerf[i_s] = new TH1F(Form("h_costs_nointerf_%s",MCSampleName[i_s].Data()), Form("h_costs_nointerf_%s",MCSampleName[i_s].Data()), ncat, 0.0, (double)ncat/10 );
    h_costs_nointerf[i_s]->Sumw2();
    h_costs_VarHi[i_s] = new TH1F(Form("h_costs_VarHi_%s",MCSampleName[i_s].Data()), Form("h_costs_VarHi_%s",MCSampleName[i_s].Data()), ncat, 0.0, (double)ncat/10 );
    h_costs_VarHi[i_s]->Sumw2();
    h_costs_VarLo[i_s] = new TH1F(Form("h_costs_VarLo_%s",MCSampleName[i_s].Data()), Form("h_costs_VarLo_%s",MCSampleName[i_s].Data()), ncat, 0.0, (double)ncat/10 );
    h_costs_VarLo[i_s]->Sumw2();
  }
  
  //--------------------------------------//
  // Loop over events to get reweighted costhetastar distributions:
  cout << " We have " << numev << " events to be processed " << endl;
  for( int i = 0; i < numev; i++ )
  {
    get_event(i);
    PrintProgressBar( i, numev );
        
    TString sample_type = SpinSampleType();
    int sample_index = GetSampleIndex( sample_type );
        
    // Only use 125 GeV and 125.5 GeV samples:
    if( !(p->Higgs_truth_mass>124 && p->Higgs_truth_mass<127) ) continue;
    
    // Apply cuts and get a cut-flow table for each signal:
    i_cutflow_counter[sample_index][0]++;
    // Pre-selection:
    if( !p->flag_pre ) continue;
    i_cutflow_counter[sample_index][1]++;
    // Photon ID:
    if( !p->flag_PID ) continue;
    i_cutflow_counter[sample_index][2]++;
    // Isolation:
    if( !p->flag_iso ) continue;
    i_cutflow_counter[sample_index][3]++;
    
    // need to scale mass:
    // Master mass values found in spinv1_Master.hh:
    double mass_scale_factor = 1.0;
    if( sample_type.Contains("POW") ) mass_scale_factor = MasterHiggsMass / MasterPOWMass;
    else if( sample_type.Contains("MG5") ) mass_scale_factor = MasterHiggsMass / MasterMG5Mass;
    else if( sample_type.Contains("AMC") ) mass_scale_factor = MasterHiggsMass / MasterAMCMass;
    
    TLorentzVector photon_1st, photon_2nd, diphoton;
    double e1_corr = mass_scale_factor * p->ph_Ecorr_1st;
    double e2_corr = mass_scale_factor * p->ph_Ecorr_2nd;
    double pt1_corr = e1_corr / cosh(p->ph_eta_1st);
    double pt2_corr = e2_corr / cosh(p->ph_eta_2nd);
    photon_1st.SetPtEtaPhiE( pt1_corr, p->ph_eta_1st, p->ph_phi_1st, e1_corr );
    photon_2nd.SetPtEtaPhiE( pt2_corr, p->ph_eta_2nd, p->ph_phi_2nd, e2_corr );
    diphoton = photon_1st + photon_2nd;
    double diphoton_mass = diphoton.M()/1000;
    double pt_mgg_ratio1 = (photon_1st.Pt()/1000.) / diphoton_mass;
    double pt_mgg_ratio2 = (photon_2nd.Pt()/1000.) / diphoton_mass;
    
    // Cut on pT:
    if( pt_mgg_ratio1 < 0.35 || pt_mgg_ratio2 < 0.25 ) continue;
    i_cutflow_counter[sample_index][4]++;
    // Cut on diphoton mass:
    if( diphoton.M()/1000.0 < 105.0 || diphoton.M()/1000.0 > 160.0 ) continue;
    i_cutflow_counter[sample_index][5]++;
    
    // Recalculate cos(theta*):
    double Q = diphoton.M();
    double Qt = diphoton.Pt();
    double costhetastarCS = TMath::Abs( 1 / (Q*sqrt(Q*Q+Qt*Qt))*((photon_1st.E()+photon_1st.Z())*(photon_2nd.E()-photon_2nd.Z()) - (photon_1st.E()-photon_1st.Z())*(photon_2nd.E()+photon_2nd.Z()) ) ); 
    
    // get the costhetastar bin:
    int bin_CTS = GetCTSPTtBin( costhetastarCS, diphoton.Pt() );
    
    double evt_weight;
    if( sample_type.Contains("POW") ) evt_weight = GetWeight( lumi, true, true );
    else // Adding the interference weight by hand!:
    {
      evt_weight = GetWeight( lumi, false, true );
      double interference_value = ss_tool->GetValue( Form("%s_interf",sample_type.Data()), bin_CTS, energy );
      int interference_sign = ss_tool->GetSign( Form("%s_interf",sample_type.Data()), bin_CTS, energy );
      evt_weight *= ( 1.0 + (interference_sign*interference_value) );
    }
    double evt_weight_nointerf = GetWeight( lumi, false, true );
    
    // Special use of the weights:
    if( sample_type.Contains("POW") )
    {
      h_costs_interf[sample_index]->AddBinContent( bin_CTS+1, evt_weight );
      h_costs_nointerf[sample_index]->Fill( bin_CTS+1, evt_weight_nointerf );
      total_yield[sample_index] += evt_weight;
      binned_yield[sample_index][bin_CTS] += evt_weight;
    }
    else // don't use interference nominally for spin 2:
    {
      h_costs_interf[sample_index]->AddBinContent( bin_CTS+1, evt_weight );
      h_costs_nointerf[sample_index]->Fill( bin_CTS+1, evt_weight_nointerf );
      total_yield[sample_index] += evt_weight_nointerf;
      binned_yield[sample_index][bin_CTS] += evt_weight_nointerf;
    }
  }
  
  //--------------------------------------------//
  // End of event loop!
  
  // Get the Yield in each category:
  GetBinYields( ncat );
  // Normalize all distributions:
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    h_costs_interf[i_s]->Scale(1.0/h_costs_interf[i_s]->Integral());
    h_costs_nointerf[i_s]->Scale(1.0/h_costs_nointerf[i_s]->Integral());
  }
  
  // Get the interference systematics:
  int nbinsc = h_costs_interf[0]->GetNbinsX();
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    for( int i_b = 0; i_b <= nbinsc; i_b++ )
    {
      double interference_diff =  h_costs_interf[i_s]->GetBinContent(i_b) - h_costs_nointerf[i_s]->GetBinContent(i_b);
      h_costs_VarHi[i_s]->SetBinContent( i_b, h_costs_interf[i_s]->GetBinContent(i_b) + interference_diff );
      h_costs_VarLo[i_s]->SetBinContent( i_b, h_costs_interf[i_s]->GetBinContent(i_b) - interference_diff );
      
      h_costs_interf[i_s]->SetBinError( i_b, interference_diff );
      h_costs_nointerf[i_s]->SetBinError( i_b, interference_diff );
    }
    // scale the systematics histograms to the same yield as the nominal:
    h_costs_VarHi[i_s]->Scale( h_costs_interf[i_s]->Integral() / h_costs_VarHi[i_s]->Integral() );
    h_costs_VarLo[i_s]->Scale( h_costs_interf[i_s]->Integral() / h_costs_VarLo[i_s]->Integral() );
  }
  
  //--------------------------------------------//
  // Plot distributions:
  can = new TCanvas("can","can");
  plot_CosThetaStar();
    
  //--------------------------------------------//
  // Create Root file to save histograms:
  TFile *sf = new TFile(Form("%s/FinalSignal/spin_parameterization_histograms.root",output_directory.Data()),"RECREATE");
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    h_costs_interf[i_s]->Write();
    h_costs_nointerf[i_s]->Write();
    h_costs_VarHi[i_s]->Write();
    h_costs_VarLo[i_s]->Write();
  }
  sf->Write();
  sf->Close();
  
  //--------------------------------------------//
  // Then open signal param files, combine with new yields:
  PrintFinalSignalInputs();
  
  //--------------------------------------------//
  // Summarize data after event loop
  cout << ". . . . . . . . . . . . . . . . . " << endl;
  cout << "  " << endl;
  cout << "Yields for each sample directly from MC: " << endl;
  for( int i_s = 0; i_s < NumMCSamples; i_s++ ) cout << "  " << MCSampleName[i_s] << " \t= " << total_yield[i_s] << endl;
  cout << "  " << endl;
  
  cout << ". . . . . . . . . . . . . . . . . " << endl;
  cout << "Bin-by-bin yields for each sample: " << endl;
  double checksum[NumMCSamples] = {0.0};
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    cout << "  " << MCSampleName[i_s] << " \t= ";
    for( int i_b = 0; i_b < ncat; i_b++ )
    {
      cout << binned_yield[i_s][i_b] << " ";
      checksum[i_s] += binned_yield[i_s][i_b];
    }
    cout << "-> SUM: " << checksum[0] << endl;
  }
  cout << " " << endl;
  
  cout << ". . . . . . . . . . . . . . . . . " << endl;
  cout << "Printing cutflows:" << endl;
  for( int i_s = 0; i_s < NumMCSamples; i_s++ )
  {
    cout << "  " << MCSampleName[i_s] << " \t: ";
    for( int i_c = 0; i_c < 6; i_c++ ) cout << i_cutflow_counter[i_s][i_c] << " \t";
  }
  cout << " " << endl;
  cout << " " << endl;
}
