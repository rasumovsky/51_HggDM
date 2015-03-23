/////////////////////////////////////////
//                                     //
//  SimultaneousBkg.cc                 //
//                                     //
//  Created: Andrew Hard 07/04/2013    //
//                                     //
//  Perform various fits to background //
//  and to the S+B in sidebands or the //
//  signal region in order to find the //
//  systematics on the background.     //
//                                     //
/////////////////////////////////////////

#include "SimultaneousBkg.hh"
#include "TStyle.h"

int main (int argc, char **argv)
{ 
  // Root macros:
  SetAtlasStyle();
  
  //--------------------------------------//
  // Check arguments:
  if( argc < 2 )
  {
    printf("\nUsage: %s <jobname> <input list>\n\n",argv[0]);
    exit(0);
  }
  TString input_directory = argv[1];
  cout << "Reading data points from: " << endl;
  cout << "   --> " << input_directory << endl;
  cout << " " << endl;
  
  // Create output directory:
  output_directory = Form("output/SimultaneousBkg/");
  system(Form("mkdir -vp %s",output_directory.Data()));
  system(Form("rm %s/*.eps",output_directory.Data()));
  
  // Initialize the canvas:
  TCanvas *can = new TCanvas("can","can",2400,900);
  
  // Create a counter for the categories & regions
  double counter[NBINS][4] = {{0.0}};
  
  //--------------------------------------//
  // Create the masses and data sets for the fits:
  RooRealVar *mass[NBINS][4];
  RooDataSet *data[NBINS][4];
  
  double bounds[4][2] = { { 105.0,        signal_band1 },
			  { signal_band2, 160.0        },
			  { 105.0,        160.0        },
			  { signal_band1, signal_band2 } };
  
  //--------------------------------------//
  // loop over categories:
  for( int i_b = 0; i_b < NBINS; i_b++ )
  {
    cout << "Currently working on bin " << i_b << endl;
        
    // loop over regions:
    for( int i_r = 0; i_r < 4; i_r++ )
    {
      mass[i_b][i_r] = new RooRealVar( Form("mass_%i_%i",i_b,i_r), Form("mass_%i_%i",i_b,i_r), bounds[i_r][0], bounds[i_r][1] );
      data[i_b][i_r] = new RooDataSet( Form("data_%i_%i",i_b,i_r), Form("data_%i_%i",i_b,i_r), *mass[i_b][i_r] ); 
    }
    
    // Read input masses from text files:
    ifstream fdata( Form("%s/mass_cate_%i.txt",input_directory.Data(),i_b+1) );
    double readout_mass = 0.0;
    while( fdata >> readout_mass )
    {
      if( readout_mass < 105.0 ) continue;
      else if( readout_mass > 160.0 ) continue;
      
      // for a test:
      //if( readout_mass > 122 && readout_mass < 130 ) continue;
      
      // Region: SB1
      else if( readout_mass < signal_band1 )
      {
	mass[i_b][0]->setVal(readout_mass);
	data[i_b][0]->add(*mass[i_b][0]);
	counter[i_b][0]++;
      }
      // Region: SB2
      else if( readout_mass > signal_band2 )
      {
	mass[i_b][1]->setVal(readout_mass);
	data[i_b][1]->add(*mass[i_b][1]);
	counter[i_b][1]++;
      }
      // Region: SR
      else if( readout_mass > signal_band1 && readout_mass < signal_band2 )
      {
	mass[i_b][3]->setVal(readout_mass);
	data[i_b][3]->add(*mass[i_b][3]);
	counter[i_b][3]++;
      }
      // Region: ALL
      mass[i_b][2]->setVal(readout_mass);
      data[i_b][2]->add(*mass[i_b][2]);
      counter[i_b][2]++;
    }
    fdata.close();
  }
  
  //--------------------------------------------//
  // Prepare for loop over fit types:
  // clear vectors:
  vector_fitnames.clear();
  vector_background.clear();
  vector_coeff_Background.clear();
  vector_coeff_Signal.clear();
  vector_nbkg.clear();
  
  // Load signal parameterizations from file:
  LoadSignalParameterization("/xdata05/andrew/SpinFiles/From_Twiki_Moriond/CamilaSignal/textFiles");
  
  // Create fit names:
  TString fit_types[3] = { "Bkg_RegALL",
			   "Bkg_RegSB",
			   "SigBkg_RegALL" };
  
  TString fit_functions[6] = { "Exppol_ORDER1",
			       "Exppol_ORDER2",
			       "Exppol_ORDER3",
			       "Bern_ORDER2",
			       "Bern_ORDER3",
			       "Bern_ORDER4" };
  
  for( int i_t = 0; i_t < 3; i_t++ )
    for( int i_f = 0; i_f < 6; i_f++ )
      vector_fitnames.push_back(Form("%s_%s",fit_types[i_t].Data(),fit_functions[i_f].Data()));
  
  //--------------------------------------------//
  // Loop over fit types:
  for( int i_v = 0; i_v < (int)vector_fitnames.size(); i_v++ )
  {
    // Then create the variables for the fits:
    RooRealVar *CB_mu[NBINS];
    RooRealVar *CB_sigma[NBINS];
    RooRealVar *CB_alpha[NBINS];
    RooRealVar *CB_n[NBINS];
    RooRealVar *GA_sigma[NBINS];
    RooRealVar *fraction[NBINS];
    
    RooRealVar p0( "p0", "p0", 1 );
    RooRealVar p1( "p1", "p1", 0.1, 0.0, 10.0 );
    RooRealVar p2( "p2", "p2", 0.1, 0.0, 10.0 );
    RooRealVar p3( "p3", "p3", 0.1, 0.0, 10.0 );
    RooRealVar p4( "p4", "p4", 0.1, 0.0, 10.0 );
    RooRealVar p5( "p5", "p5", 0.1, 0.0, 10.0 );
    RooRealVar p6( "p6", "p6", 0.1, 0.0, 10.0 );
    RooRealVar p7( "p7", "p7", 0.1, 0.0, 10.0 );
    RooRealVar p8( "p8", "p8", 0.1, 0.0, 10.0 );
    
    RooConstVar min("min","min",105.0);
    RooConstVar max("max","max",160.0);
    
    RooRealVar c1("c1", "c1", 0.0, -1.0, 1.0 );
    RooRealVar c2("c2", "c2", 0.0, -1.0, 1.0 );
    RooRealVar c3("c3", "c3", 0.0, -1.0, 1.0 );
    RooRealVar c4("c4", "c4", 0.0, -1.0, 1.0 );
    RooRealVar c5("c5", "c5", 0.0, -1.0, 1.0 );
    RooRealVar c6("c6", "c6", 0.0, -1.0, 1.0 );
    RooRealVar c7("c7", "c7", 0.0, -1.0, 1.0 );
    RooRealVar c8("c8", "c8", 0.0, -1.0, 1.0 );
    
    // Declare array of functions:
    RooCBShape *CB[NBINS][3] = {{NULL}};
    RooGaussian *GA[NBINS][3] = {{NULL}};
    RooAddPdf *Signal[NBINS][3] = {{NULL}};
        
    RooBernsteinM *Bern1[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern2[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern3[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern4[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern5[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern6[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern7[NBINS][3] = {{NULL}};
    RooBernsteinM *Bern8[NBINS][3] = {{NULL}};
        
    RooGenericPdf *Exppol1[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol2[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol3[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol4[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol5[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol6[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol7[NBINS][3] = {{NULL}};
    RooGenericPdf *Exppol8[NBINS][3] = {{NULL}};
    
    RooAbsPdf *Background[NBINS][3] = {{NULL}};
    RooExtendPdf *Background_Extended[NBINS][3] = {{NULL}};
    RooAbsPdf *Model[NBINS][3] = {{NULL}};
    RooAddPdf *AddModel[NBINS][3] = {{NULL}};
    
    //RooRealVar *fractionSB[NBINS];
    RooRealVar *coeff_Signal[NBINS];
    RooRealVar *coeff_Background[NBINS];
    RooRealVar *nbkg[NBINS];
    
    //----------------------------------------//
    // Prepare simultaneous fit:
    map<string,RooDataSet*> datasetMap;
    RooCategory *data_category = new RooCategory("data_category","data_category");
    RooSimultaneous simultaneousPDF("simultaneousPDF","simultaneousPDF",*data_category);
    simultaneousPDF.Print();
    RooArgSet *args = new RooArgSet();
    
    //----------------------------------------//
    // loop over the bins to define the functions:
    for( int i_b = 0; i_b < NBINS; i_b++ )
    { 
      int i_s = MassToIndex( 126.5 );
      CB_mu[i_b] = new RooRealVar(Form("CB_mu_%i",i_b),Form("CB_mu_%i",i_b),param_values[i_b][i_s][2]);
      CB_sigma[i_b] = new RooRealVar(Form("CB_sigma_%i",i_b),Form("CB_sigma_%i",i_b),param_values[i_b][i_s][3]);
      CB_alpha[i_b] = new RooRealVar(Form("CB_alpha_%i",i_b),Form("CB_alpha_%i",i_b),param_values[i_b][i_s][4]);
      CB_n[i_b] = new RooRealVar(Form("CB_n_%i",i_b),Form("CB_n_%i",i_b),param_values[i_b][i_s][5]);
      GA_sigma[i_b] = new RooRealVar(Form("GA_sigma_%i",i_b),Form("GA_sigma_%i",i_b),param_values[i_b][i_s][7]);
      fraction[i_b] = new RooRealVar(Form("fraction_%i",i_b),Form("fraction_%i",i_b),param_values[i_b][i_s][8]);
      
      //fractionSB[i_b] = new RooRealVar(Form("fractionSB_%i",i_b), Form("fractionSB_%i",i_b), 0.01, 0.0, 1.0);
      coeff_Signal[i_b] = new RooRealVar(Form("coeff_Signal_%i",i_b), Form("coeff_Signal_%i",i_b), 100, -10000.0, 10000);
      coeff_Background[i_b]  = new RooRealVar(Form("coeff_Background_%i",i_b), Form("coeff_Background_%i",i_b), 1000, 0.0, 100000);
      nbkg[i_b] = new RooRealVar(Form("nbkg_%i",i_b), Form("nbkg_%i",i_b), 10000, 0.0, 1000000);
      
      // loop over the regions:
      for( int i_r = 0; i_r < 3; i_r++ )
      { 
	// NEED TO CHANGE NAMES BASED ON CATEGORY
	CB[i_b][i_r] = new RooCBShape(Form("CB_%i_%i",i_b,i_r),Form("CB_%i_%i",i_b,i_r),*mass[i_b][i_r],*CB_mu[i_b],*CB_sigma[i_b],*CB_alpha[i_b],*CB_n[i_b]);
	GA[i_b][i_r] = new RooGaussian(Form("GA_%i_%i",i_b,i_r),Form("GA_%i_%i",i_b,i_r),*mass[i_b][i_r],*CB_mu[i_b],*GA_sigma[i_b]);
	Signal[i_b][i_r] = new RooAddPdf(Form("Signal_%i_%i",i_b,i_r),Form("Signal_%i_%i",i_b,i_r),RooArgList(*CB[i_b][i_r],*GA[i_b][i_r]),*fraction[i_b]);
	
	Bern1[i_b][i_r] = new RooBernsteinM(Form("Bern1_%i_%i",i_b,i_r), Form("Bern1_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1), &min, &max);
	Bern2[i_b][i_r] = new RooBernsteinM(Form("Bern2_%i_%i",i_b,i_r), Form("Bern2_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2), &min, &max);
	Bern3[i_b][i_r] = new RooBernsteinM(Form("Bern3_%i_%i",i_b,i_r), Form("Bern3_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3), &min, &max);
	Bern4[i_b][i_r] = new RooBernsteinM(Form("Bern4_%i_%i",i_b,i_r), Form("Bern4_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3,p4), &min, &max);
	Bern5[i_b][i_r] = new RooBernsteinM(Form("Bern5_%i_%i",i_b,i_r), Form("Bern5_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3,p4,p5), &min, &max);
	Bern6[i_b][i_r] = new RooBernsteinM(Form("Bern6_%i_%i",i_b,i_r), Form("Bern6_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3,p4,p5,p6), &min, &max);
	Bern7[i_b][i_r] = new RooBernsteinM(Form("Bern7_%i_%i",i_b,i_r), Form("Bern7_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3,p4,p5,p6,p7), &min, &max);
	Bern8[i_b][i_r] = new RooBernsteinM(Form("Bern8_%i_%i",i_b,i_r), Form("Bern8_%i_%i",i_b,i_r), *mass[i_b][i_r], RooArgList(p0,p1,p2,p3,p4,p5,p6,p7,p8), &min, &max);
	
	Exppol1[i_b][i_r] = new RooGenericPdf(Form("Exppol1_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) )",RooArgList(*mass[i_b][i_r],c1));
	Exppol2[i_b][i_r] = new RooGenericPdf(Form("Exppol2_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2));
	Exppol3[i_b][i_r] = new RooGenericPdf(Form("Exppol3_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3));
	Exppol4[i_b][i_r] = new RooGenericPdf(Form("Exppol4_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) + @4*(@0-100)*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3,c4));
	Exppol5[i_b][i_r] = new RooGenericPdf(Form("Exppol5_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) + @4*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @5*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3,c4,c5));
	Exppol6[i_b][i_r] = new RooGenericPdf(Form("Exppol6_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) + @4*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @5*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @6*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3,c4,c5,c6));
	Exppol7[i_b][i_r] = new RooGenericPdf(Form("Exppol7_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) + @4*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @5*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @6*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @7*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3,c4,c5,c6,c7));
	Exppol8[i_b][i_r] = new RooGenericPdf(Form("Exppol8_%i_%i",i_b,i_r),"TMath::Exp( @1*(@0-100) + @2*(@0-100)*(@0-100) + @3*(@0-100)*(@0-100)*(@0-100) + @4*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @5*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @6*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @7*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) + @8*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100)*(@0-100) )",RooArgList(*mass[i_b][i_r],c1,c2,c3,c4,c5,c6,c7,c8));
	
	// Background function:
	if( vector_fitnames[i_v].Contains("Bern_ORDER1") ){ Background[i_b][i_r] = Bern1[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER2") ){ Background[i_b][i_r] = Bern2[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER3") ){ Background[i_b][i_r] = Bern3[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER4") ){ Background[i_b][i_r] = Bern4[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER5") ){ Background[i_b][i_r] = Bern5[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER6") ){ Background[i_b][i_r] = Bern6[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER7") ){ Background[i_b][i_r] = Bern7[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Bern_ORDER8") ){ Background[i_b][i_r] = Bern8[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER1") ){ Background[i_b][i_r] = Exppol1[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER2") ){ Background[i_b][i_r] = Exppol2[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER3") ){ Background[i_b][i_r] = Exppol3[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER4") ){ Background[i_b][i_r] = Exppol4[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER5") ){ Background[i_b][i_r] = Exppol5[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER6") ){ Background[i_b][i_r] = Exppol6[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER7") ){ Background[i_b][i_r] = Exppol7[i_b][i_r]; }
	else if( vector_fitnames[i_v].Contains("Exppol_ORDER8") ){ Background[i_b][i_r] = Exppol8[i_b][i_r]; }
	
	//mass[i_b][i_r]->setRange(Form("rangeName_%i_%i",i_b,i_r),bounds[i_r][0],bounds[i_r][1]);
	
	// Make the background extendable:
	Background_Extended[i_b][i_r] = new RooExtendPdf(Form("Background_Extended_%i_%i",i_b,i_r),Form("Background_Extended_%i_%i",i_b,i_r),*Background[i_b][i_r],*nbkg[i_b]);
	
	// For signal + background fit:
	if( vector_fitnames[i_v].Contains("SigBkg") )
	{
	  //AddModel[i_b][i_r] = new RooAddPdf(Form("Model_%i_%i",i_b,i_r),Form("Model_%i_%i",i_b,i_r),RooArgList(*Signal[i_b][i_r],*Background[i_b][i_r]),*fractionSB[i_b]);
	  //AddModel[i_b][i_r] = new RooAddPdf(Form("Model_%i_%i",i_b,i_r),Form("Model_%i_%i",i_b,i_r),RooArgList(*Signal[i_b][i_r],*Background[i_b][i_r]),RooArgList(*coeff_Signal[i_b],*coeff_Background[i_b]));
	  AddModel[i_b][i_r] = new RooAddPdf(Form("Model_%i_%i",i_b,i_r),Form("Model_%i_%i",i_b,i_r),RooArgList(*Signal[i_b][i_r],*Background_Extended[i_b][i_r]),RooArgList(*coeff_Signal[i_b],*coeff_Background[i_b]));
	  Model[i_b][i_r] = AddModel[i_b][i_r];
	}
	else
	{
	  //Model[i_b][i_r] = Background[i_b][i_r];
	  Model[i_b][i_r] = Background_Extended[i_b][i_r];
	}
	
	data_category->defineType(Form("BIN%i_Reg%i",i_b,i_r));
	data_category->setRange(Form("rangeName_%i_%i",i_b,i_r),Form("BIN%i_Reg%i",i_b,i_r));
	
	// in case of sideband fit:
	if( vector_fitnames[i_v].Contains("RegSB") && ( i_r == 0 || i_r == 1 ) )
	{
	  simultaneousPDF.addPdf(*Model[i_b][i_r],Form("BIN%i_Reg%i",i_b,i_r));
	  datasetMap[Form("BIN%i_Reg%i",i_b,i_r)] = data[i_b][i_r];
	  args->add(RooArgSet(*mass[i_b][i_r]));
	}
	// fit to full range:
	else if( vector_fitnames[i_v].Contains("RegALL") && i_r == 2 )
	{
	  simultaneousPDF.addPdf(*Model[i_b][i_r],Form("BIN%i_Reg%i",i_b,i_r));
	  datasetMap[Form("BIN%i_Reg%i",i_b,i_r)] = data[i_b][i_r];
	  args->add(RooArgSet(*mass[i_b][i_r]));
	}
      }
    }
    
    RooDataSet *obsData = new RooDataSet("obsData", "combined data", *args, Index(*data_category), Import(datasetMap) );
    // Do the fit:
    simultaneousPDF.fitTo( *obsData, Hesse(0), Save(true), PrintLevel(FIT_PrintLevel), Verbose(FIT_Verbose) );
    
    //--------------------------------------------//
    // Count the background events per category:
    vector<double> current_bkg_vector; current_bkg_vector.clear();
    vector<double> current_vector_coeff_Background; current_vector_coeff_Background.clear();
    vector<double> current_vector_coeff_Signal; current_vector_coeff_Signal.clear();
    vector<double> current_vector_nbkg; current_vector_nbkg.clear();
    
    for( int i_b = 0; i_b < NBINS; i_b++ )
    {
      double bkg_prediction = 0;
      
      mass[i_b][2]->setRange("range_ALL",105.0,160.0);
      mass[i_b][2]->setRange("range_SR",signal_band1,signal_band2);
      //RooAbsReal *integral_Bkg_ALL = (RooAbsReal*)Background[i_b][2]->createIntegral(RooArgSet(*mass[i_b][2]), NormSet(*mass[i_b][2]), Range("range_ALL"));
      //RooAbsReal *integral_Bkg_SR = (RooAbsReal*)Background[i_b][2]->createIntegral(RooArgSet(*mass[i_b][2]), NormSet(*mass[i_b][2]), Range("range_SR"));
      RooAbsReal *integral_Bkg_ALL = (RooAbsReal*)Background_Extended[i_b][2]->createIntegral(RooArgSet(*mass[i_b][2]), NormSet(*mass[i_b][2]), Range("range_ALL"));
      RooAbsReal *integral_Bkg_SR = (RooAbsReal*)Background_Extended[i_b][2]->createIntegral(RooArgSet(*mass[i_b][2]), NormSet(*mass[i_b][2]), Range("range_SR"));
      
      // Fit in sidebands (requires categories):
      if( vector_fitnames[i_v].Contains("RegSB") )
      {
	// assumes a background-only fit for now:
	cout << "    Fit region: RegSB" << endl;
	double fraction_SR = integral_Bkg_SR->getVal() / integral_Bkg_ALL->getVal();
	bkg_prediction = ( fraction_SR / ( 1.0 - fraction_SR ) ) * ( data[i_b][0]->sumEntries() + data[i_b][1]->sumEntries() );
      }
      // Fit over full spectrum:
      else if( vector_fitnames[i_v].Contains("RegALL") )
      {
	cout << "    Fit region: RegALL" << endl;
	// signal+background fit:
	if( vector_fitnames[i_v].Contains("SigBkg") )
	{
	  cout << "    Fit model: SigBkg" << endl;
	  cout << "fractionSB.getVal(): " << coeff_Signal[i_b]->getVal() / coeff_Background[i_b]->getVal() << endl; //fractionSB[i_b]->getVal() << endl;
	  //bkg_prediction = ( integral_Bkg_SR->getVal() / integral_Bkg_ALL->getVal() ) * ( (1.0-fractionSB[i_b]->getVal()) * data[i_b][2]->sumEntries() );
	  bkg_prediction = ( integral_Bkg_SR->getVal() / integral_Bkg_ALL->getVal() ) * ( (1.0-(coeff_Signal[i_b]->getVal() / coeff_Background[i_b]->getVal())) * data[i_b][2]->sumEntries() );
	}
	// background-only fit:
	else if( vector_fitnames[i_v].Contains("Bkg") )
	{
	  cout << "    Fit model: Bkg" << endl;
	  bkg_prediction = ( integral_Bkg_SR->getVal() / integral_Bkg_ALL->getVal() ) * data[i_b][2]->sumEntries();
	}
      }
      current_bkg_vector.push_back(bkg_prediction);
      current_vector_coeff_Background.push_back(coeff_Background[i_b]->getVal());
      current_vector_coeff_Signal.push_back(coeff_Signal[i_b]->getVal());
      current_vector_nbkg.push_back(nbkg[i_b]->getVal());
    }
    vector_background.push_back(current_bkg_vector);
    vector_coeff_Background.push_back(current_vector_coeff_Background);
    vector_coeff_Signal.push_back(current_vector_coeff_Signal);
    vector_nbkg.push_back(current_vector_nbkg);

    //--------------------------------------------//
    // Plot the results:
    cout << "Starting to make the plots" << endl;
    can->Clear();
    can->Divide(5,2);
    
    // loop over bins to plot results:
    cout << "  DoFits::Plotting the result" << endl;
    for( int i_b = 0; i_b < NBINS; i_b++ )
    {
      mass[i_b][2]->setRange("range_SB1",105,signal_band1);
      mass[i_b][2]->setRange("range_SB2",signal_band2,160);
      mass[i_b][2]->setRange("range_ALL",105,160);
      
      can->cd(i_b+1);
      RooPlot *frame = mass[i_b][2]->frame();
      frame->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
      frame->GetXaxis()->Set(55,105,160);
      frame->GetYaxis()->SetTitle("Events / GeV");
      data[i_b][2]->plotOn(frame, LineColor(kBlack), MarkerColor(kBlack));
      if( vector_fitnames[i_v].Contains("RegALL") )
	Model[i_b][2]->plotOn(frame,Range("range_ALL"),NormRange("range_ALL"),LineColor(kRed));
      else if( vector_fitnames[i_v].Contains("RegSB") )
	Model[i_b][2]->plotOn(frame,Range("range_SB1,range_SB2"),NormRange("range_SB1,range_SB2"),LineColor(kRed));
      frame->Draw();
      
      TLine *line_SR1 = new TLine(); line_SR1->SetLineStyle(1); line_SR1->SetLineWidth(2); line_SR1->SetLineColor(kBlack);
      TLine *line_SR2 = new TLine(); line_SR2->SetLineStyle(1); line_SR2->SetLineWidth(2); line_SR2->SetLineColor(kBlack);
      line_SR1->DrawLine( signal_band1, frame->GetMinimum(), signal_band1, frame->GetMaximum() );
      line_SR2->DrawLine( signal_band2, frame->GetMinimum(), signal_band2, frame->GetMaximum() );
      can->Print(Form("%s/plot_%s.eps",output_directory.Data(),vector_fitnames[i_v].Data()));
    }
  }
  
  //--------------------------------------------//
  // Create histogram of the difference in background estimates (%):
  can->cd();
  can->Clear();
  can->Divide(5,2);
  TH1F *h_background_difference[NBINS+1];
  TH1F *h_bkgdiff_over_expsig[NBINS+1];
  
  h_background_difference[NBINS] = new TH1F(Form("h_background_difference_%i",NBINS),Form("h_background_difference_%i",NBINS),20,-10,10);
  h_background_difference[NBINS]->GetXaxis()->SetTitle("% difference from mean");
  h_background_difference[NBINS]->GetYaxis()->SetTitle("Entries");
  
  h_bkgdiff_over_expsig[NBINS] = new TH1F(Form("h_bkgdiff_over_expsig_%i",NBINS),Form("h_bkgdiff_over_expsig_%i",NBINS),20,-200,200);
  h_bkgdiff_over_expsig[NBINS]->GetXaxis()->SetTitle("Difference from mean over expected signal %");
  h_bkgdiff_over_expsig[NBINS]->GetYaxis()->SetTitle("Entries");
  
  // loop over bins:
  for( int i_b = 0; i_b < NBINS; i_b++ )
  {
    h_background_difference[i_b] = new TH1F(Form("h_background_difference_%i",i_b),Form("h_background_difference_%i",i_b),20,-10,10);
    h_background_difference[i_b]->GetXaxis()->SetTitle("% difference from mean");
    h_background_difference[i_b]->GetYaxis()->SetTitle("Entries");
    
    h_bkgdiff_over_expsig[i_b] = new TH1F(Form("h_bkgdiff_over_expsig_%i",i_b),Form("h_bkgdiff_over_expsig_%i",i_b),20,-200,200);
    h_bkgdiff_over_expsig[i_b]->GetXaxis()->SetTitle("Difference from mean over expected signal %");
    h_bkgdiff_over_expsig[i_b]->GetYaxis()->SetTitle("Entries");
    
    // loop over fits to get the mean value:
    double sum = 0.0;
    for( int i_v = 0; i_v < vector_fitnames.size(); i_v++ ) sum += vector_background[i_v][i_b];
    double average = sum / ( (double)vector_fitnames.size() );
    
    // then loop over fits, getting difference and adding it to the proper histogram:
    for( int i_v = 0; i_v < vector_fitnames.size(); i_v++ )
    {
      double difference = vector_background[i_v][i_b] - average;
      h_background_difference[i_b]->Fill( 100*difference/average );
      h_background_difference[NBINS]->Fill( 100*difference/average );
      
      int i_s = MassToIndex( 126.5 );
      double expected_signal = param_values[i_b][i_s][1];
      h_bkgdiff_over_expsig[i_b]->Fill( 100*difference/expected_signal );
      h_bkgdiff_over_expsig[NBINS]->Fill( 100*difference/expected_signal );
    }
    
    can->cd(i_b+1);
    h_background_difference[i_b]->Draw();
    
  }
  can->Print(Form("%s/hist_bkgdiff_over_mean_individual.eps",output_directory.Data()));
  
  can->Clear();
  can->Divide(5,2);
  for( int i_b = 0; i_b < NBINS; i_b++ )
  {
    can->cd(i_b+1);
    h_bkgdiff_over_expsig[i_b]->Draw();
  }
  can->Print(Form("%s/hist_bkgdiff_over_expsig_individual.eps",output_directory.Data()));
  
  TCanvas *can2 = new TCanvas("can2","can2",800,800);
  can2->cd();
  can2->Clear();
  h_background_difference[NBINS]->Draw();
  can2->Print(Form("%s/hist_bkgdiff_over_mean_Total.eps",output_directory.Data()));
  
  can2->Clear();
  h_bkgdiff_over_expsig[NBINS]->Draw();
  can2->Print(Form("%s/hist_bkgdiff_over_expsig_Total.eps",output_directory.Data()));
  
  
  /*
  //--------------------------------------------//
  // Loop over categories:
  for( int i_b = 0; i_b < NBINS; i_b++ )
  {
    can->Clear();
    
    mass[i_b][2]->setRange("range_SB1",105,signal_band1);
    mass[i_b][2]->setRange("range_SB2",signal_band2,160);
    mass[i_b][2]->setRange("range_ALL",105,160);
    
    RooPlot *frame = mass[i_b][2]->frame();
    frame->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    frame->GetXaxis()->Set(55,105,160);
    frame->GetYaxis()->SetTitle("Events / GeV");
    data[i_b][2]->plotOn(frame, LineColor(kBlack), MarkerColor(kBlack));
    // then loop over fit types:
    for( int i_v = 0; i_v < vector_fitnames.size(); i_v++ )
    {
      if( vector_fitnames[i_v].Contains("RegALL") )
	SAVE_PDFs[i_v][i_b]->plotOn(frame,Range("range_ALL"),NormRange("range_ALL"),LineColor(i_v+2));
      else if( vector_fitnames[i_v].Contains("RegSB") )
	SAVE_PDFs[i_v][i_b]->plotOn(frame,Range("range_SB1,range_SB2"),NormRange("range_SB1,range_SB2"),LineColor(i_v+2));
    }
    frame->Draw();
    
    TLine *line_SR1 = new TLine(); line_SR1->SetLineStyle(1); line_SR1->SetLineWidth(2); line_SR1->SetLineColor(kBlack);
    TLine *line_SR2 = new TLine(); line_SR2->SetLineStyle(1); line_SR2->SetLineWidth(2); line_SR2->SetLineColor(kBlack);
    line_SR1->DrawLine( signal_band1, frame->GetMinimum(), signal_band1, frame->GetMaximum() );
    line_SR2->DrawLine( signal_band2, frame->GetMinimum(), signal_band2, frame->GetMaximum() );
    can->Print(Form("%s/plot_BIN%i.eps",output_directory.Data(),i_b));
  }
  */
  
  if( vector_background.size() != vector_fitnames.size() ) cout << "Problem: Differing numbers of background estimates and fits!" << endl;
  
  //--------------------------------------------//
  // Summarize results after event loop:
  // increase precision of output:
  std::cout.precision(10);
  
  cout << "  " << endl; 
  cout << ". . . . . . . . . . . . . . . . . " << endl;
  cout << "  Printing background estimates by bin.  " << endl;
  cout << "  " << endl;
  
  for( int i_c = 0; i_c < NBINS; i_c++ )
  {
    cout << "  BIN " << i_c << endl;
    cout << "    NAME \t<< background estimate \t<< coeff_Background \t<< coeff_Signal \t<< nbkg << " << endl;
    for( int i_v = 0; i_v < (int)vector_fitnames.size(); i_v++ )
    {
      cout << "    " << vector_fitnames[i_v] << " \t" << vector_background[i_v][i_c] << " \t" << vector_coeff_Background[i_v][i_c] << " \t" << vector_coeff_Signal[i_v][i_c] << " \t" << vector_nbkg[i_v][i_c] << endl;
    }
  }
  
  cout << "  " << endl; 
  cout << ". . . . . . . . . . . . . . . . . " << endl;
  cout << "  Printing region information by category." << endl;
  cout << "  " << endl;
  for( int i_b = 0; i_b < NBINS; i_b++ )
    cout << "Category " << i_b+1 << " \tALL: " << counter[i_b][2] << " \tSB: " << counter[i_b][0] + counter[i_b][1] << " \tSR: " << counter[i_b][3] << endl;
  cout << "  " << endl;
  return 0;
}
