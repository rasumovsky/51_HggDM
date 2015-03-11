/////////////////////////////////////////
//                                     //
//  spin_masspoints.cc                 //
//                                     //
//  Author: Andrew Hard                //
//  Date: 08/10/2013                   //
//                                     //
//  Description:                       //
//  This code runs over the ntuples    //
//  for data to create the mass points //
//  to be used in the workspace. To be //
//  used for the 2013 publication of   //
//  the spin measurement of the Higgs- //
//  like particle decaying to two      //
//  photons with the ATLAS detector.   //
//                                     //
//  ./bin/spin_masspoints <1> <2>      //
//    <1> output name and jobnames.    //
//    "#Bins"                          //
//                                     //
/////////////////////////////////////////

#include "spinv1_masspoints.hh"

int main( int argc, char **argv )
{
  // Root macros:
  SetAtlasStyle();
  
  // Check arguments:
  if( argc < 4 ){ printf("\nUsage: %s <jobname> <file_list> <energy>\n\n",argv[0]); exit(0); }
  jobname = argv[1];
  TString sample_list = argv[2];
  energy = atoi(argv[3]);
  
  // Print input file:
  cout << "Reading data points from: " << endl;
  system(Form("cat %s",sample_list.Data()));
  cout << " " << endl;
  
  //--------------------------------------//
  // Get configuration from options:
  int NbinsCTS = 11;
  
  //--------------------------------------//
  // Make TChain and find number of included events:
  TChain* chain = MakeChain("Hgg_tree", sample_list, "badfile");
  p = new ntup_SpinPub(chain);
  int numev = p->fChain->GetEntries();
  
  //--------------------------------------//
  // Create output text files for spin_masspoints:
  output_directory = Form("%s/%s/mass_points_%iTeV/",master_output.Data(),jobname.Data(),energy);
  system(Form("mkdir -vp %s",output_directory.Data()));
  system(Form("mkdir -vp %s/plots",output_directory.Data()));
  ofstream mass_files[20];
  for( int i_c = 0; i_c < NbinsCTS; i_c++ )
  {
    mass_files[i_c].open(Form("%s/mass_costs%i_%iTeV.txt",output_directory.Data(),i_c+1,energy));
    mass_files[i_c].setf(ios::showpoint);
    mass_files[i_c].setf(ios::fixed);
    mass_files[i_c].precision(14);
  }
  
  ofstream inclusive_file;
  inclusive_file.open("temp_list.txt");
  
  // Relative pt histograms:
  h_relpt1 = new TH1F("h_relpt1","h_relpt1",30,0.35,1);
  h_relpt2 = new TH1F("h_relpt2","h_relpt2",30,0.25,1);
  h2_relpt = new TH2D("h2_relpt","h2_relpt",30,0.35,1,30,0.25,1);
  h_sinhdeta = new TH1F("h_sinhdeta","h_sinhdeta",30,-8,8);
  h_relpt1->Sumw2();
  h_relpt2->Sumw2();
  h_sinhdeta->Sumw2();
  
  InitializeKinematicHists();
  
  // Counting (as a cross-check printout:
  int category_counter[20] = {0.0};
  
  int cutflow_counter[10] = {0.0};
  
  //--------------------------------------//
  // Loop over events:
  cout << "We have "<< numev << " events to be processed." << endl;
  for( int index = 0; index < numev; index++ )
  {
    p->fChain->GetEntry(index);
    PrintProgressBar( index, numev );
    
    cutflow_counter[0]++;
    if( p->flag_pre )
    {
      cutflow_counter[1]++;
      if( p->flag_PID )
      {
	cutflow_counter[2]++;
	if( p->flag_iso )
	{
	  cutflow_counter[3]++;
	  if( p->flag_mgg )
	  {
	    cutflow_counter[4]++;
	  }
	}
      }
    }
    
    // Apply quality flags:
    if( !p->flag_pre || !p->flag_PID || !p->flag_iso || !p->flag_mgg ) continue;
    
    // do pT cut here:
    if( !p->flag_RelPt_35_25 ) continue;
    cutflow_counter[5]++;

    h_relpt1->Fill( p->ph_pt_corr_1st/p->mass_PV_EM );
    h_relpt2->Fill( p->ph_pt_corr_2nd/p->mass_PV_EM );
    h2_relpt->Fill( p->ph_pt_corr_1st/p->mass_PV_EM, p->ph_pt_corr_2nd/p->mass_PV_EM );
    h_sinhdeta->Fill( TMath::SinH( p->ph_eta_1st-p->ph_eta_2nd) );
    
    // get the costhetastar bin:
    //int bin_CTS = GetCTSBin( NbinsCTS, p->costhetastar_CS );
    int bin_CTS = GetCTSPTBin( NbinsCTS, p->costhetastar_CS, p->PT );
   
    // Fill the mass files and counter:
    double mass = p->mass_PV_EM/1000.0;
    
    // pT cut:
    if( p->PT > 300000 ) continue;
    cutflow_counter[6]++;
   
    // only include masses above 105 GeV:
    if( mass < 105 ) continue;
    cutflow_counter[7]++;
    
    mass_files[bin_CTS] << mass << endl;
    inclusive_file << p->Run << " " << p->Event << endl;
    category_counter[bin_CTS]++;

    TLorentzVector p1, p2;
    p1.SetPtEtaPhiE( p->ph_pt_corr_1st/1000, p->ph_eta_1st, p->ph_phi_1st, p->ph_Ecorr_1st/1000 );
    p2.SetPtEtaPhiE( p->ph_pt_corr_2nd/1000, p->ph_eta_2nd, p->ph_phi_2nd, p->ph_Ecorr_2nd/1000 );
    FillKinematicHists( p1, p2, fabs(p->costhetastar_CS), p->pt_t/1000, 1.0 );
  }
  
  //--------------------------------------//
  // Close mass point files:
  for( int i_c = 0; i_c < NbinsCTS; i_c++ ) mass_files[i_c].close();
  inclusive_file.close();
  
  //--------------------------------------//
  // Then compare two files:
  //CompareTwoFiles( "Hongtao_list.txt", "temp_list.txt" );
  
  // plot histograms:
  can = new TCanvas("can","can");
  can->cd();
  DrawRelPtPlots();
  SaveKinematicHists();
  
  //--------------------------------------//
  // Print summary:
  ofstream file_summary;
  file_summary.open(Form("%s/event_summary.txt",output_directory.Data()));
  
  cout << " " << endl;
  cout << "Mass points have been generated successsfully." << endl;
  cout << " " << endl;
  cout << "Results for " << jobname << endl;
  
  file_summary << " " << endl;
  file_summary << "Mass points have been generated successsfully." << endl;
  file_summary << " " << endl;
  file_summary << "Results for " << jobname << endl;
  int category_sum = 0;
  for( int i_c = 0; i_c < NbinsCTS; i_c++ )
  {
    cout << "  category " << i_c+1 << " : " << category_counter[i_c] << endl;
    file_summary << "  category " << i_c+1 << " : " << category_counter[i_c] << endl;
    category_sum += category_counter[i_c];
  }
  cout << "  Sum : " << category_sum << endl;
  cout << " " << endl;
  
  file_summary << "  Sum : " << category_sum << endl;
  file_summary << " " << endl;
  file_summary.close();

  // Print the cutflow:
  for( int i_c = 0; i_c < 10; i_c++ )
    cout << "  cut " << i_c << "  " << cutflow_counter[i_c] << endl;
  
  return 0;
}
