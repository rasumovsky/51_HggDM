////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This program is useful as an interface to the H->diphoton + DM analysis   //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

int main( int argc, char **argv )
{
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <jobname> <option>\n\n",argv[0]);
    exit(0);
  }
  
  TString MasterJobName = argv[1];
  TString option = argv[2];
  
  //--------------------------------------//
  // Set the analysis components to execute:
  bool MakeInputs     = option.Contains("MakeInputs");
  bool MakeSBPlots    = option.Contains("MakeSBPlots");
  bool MakeWorkspaces = option.Contains("MakeWorkspaces");
  bool ResubmitWSJobs = option.Contains("ResubmitWSJobs");
  bool RunInParallel  = true;
  
  //--------------------------------------//
  // Compile all relevant executables:
  TString exe_inputs     = "NPPV1_hf_inputs";     // "NPPV1_inputs";
  TString exe_sbplots    = "NPPV1_sbplots";       // same
  TString exe_yields     = "NPPV1_yields";        // same
  TString exe_makespace  = "NPPV1_histfactory";   // "NPPV1_makespace";
  TString exe_checkjobs  = "NPPV1_checkjobs";     // same
  
  if( MakeInputs     ){ MakeExe( exe_inputs ); }
  if( MakeSBPlots    ){ MakeExe( exe_sbplots ); MakeExe( exe_yields ); }
  if( MakeWorkspaces ){ MakeExe( exe_makespace ); }
  if( ResubmitWSJobs ){ MakeExe( exe_makespace ); MakeExe( exe_checkjobs ); }
   
  //--------------------------------------//
  // Some options for the more complicated steps:
  TString input_options = "event_based";//"object_based";
  TString ws_option = "muneg_noplot";//"fitasimov_muneg_noplot";
  
  //--------------------------------------//
  // Step 1: Make input histograms:
  if( MakeInputs )
  {
    cout << "NPPV1_Master: Step 1 - Make input histograms." << endl;
    for( int i_l = 0; i_l < nlambda_values; i_l++ )
    {
      int current_lambda = lambda_values[i_l];
      system(Form("./bin/%s %s %s %i", exe_inputs.Data(), MasterJobName.Data(), input_options.Data(), current_lambda));
    }
  }
  
  //--------------------------------------//
  // Step 2: S&B Plots & signal yield
  if( MakeSBPlots )
  {
    cout << "NPPV1_Master: Step 2 - Make plots of S&B." << endl;
    int chosen_lambda = lambda_values[16];
    int chosen_lifetime = ps_lifetimes[5];
    int chosen_lambda2 = lambda_values[18];
    int chosen_lifetime2 = ps_lifetimes[3];
    // signal and background plots;
    // this can be asimov1, asimov0, or data:
    system(Form("./bin/%s %s %i %i %i %s", exe_sbplots.Data(), MasterJobName.Data(),ncategories,chosen_lambda,chosen_lifetime,sbplot_option.Data()));
    system(Form("./bin/%s %s %i %i %i %s", exe_sbplots.Data(), MasterJobName.Data(),ncategories,chosen_lambda2,chosen_lifetime2,sbplot_option.Data()));
    // yield plot for signal
    system(Form("./bin/%s %s %i", exe_yields.Data(), MasterJobName.Data(), ncategories));
  }
  
  return 0;
}
