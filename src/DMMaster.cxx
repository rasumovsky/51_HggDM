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

/**
   This is the main DMMaster method:
*/
int main( int argc, char **argv ) {
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <jobname> <option>\n\n",argv[0]);
    exit(0);
  }
  
  TString MasterJobName = argv[1];
  TString MasterOption = argv[2];
  
  //--------------------------------------//
  // Set the analysis components to execute:
  bool MakeMassPoints = MasterOption.Contains("MakeMassPoints");
  bool MakeSigParam = MasterOption.Contains("MakeSigParam");
  
  //--------------------------------------//
  // Compile all relevant executables:
  TString exe_MassPoints = "DMMassPoints";
  TString exe_SigParam = "DMSigParam";
  if (MakeMassPoints) { MakeExe( exe_MassPoints ); }
  if (MakeSigParam) { MakeExe( exe_SigParam ); }
  
  //--------------------------------------//
  // Step 1: Make mass points:
  if (MakeMassPoints) {
    cout << "DMMaster: Step 1 - Make mass points." << endl;
  }
  
  //--------------------------------------//
  // Step 2: Make the signal parameterization:
  if (MakeSigParam) {
    cout << "DMMaster: Step 2 - Make signal parameterization." << endl;
  }
  
  return 0;
}

/**
   Compiles the executable required by the job options.
*/
void MakeExe(TString exename) {
  
  // recompile all executables before running...
  system(Form("rm bin/%s",exename.Data()));
  system(Form("make bin/%s",exename.Data()));
}

/**
   Submits the workspace jobs via bsub.
*/
void SubmitWSViaBsub(TString executable_name, TString executable_jobname, TString executable_option, int executable_lambda, int executable_lifetime) {
  
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",ws_jobscript.Data()));
  
  //change permissions for local input files:!!!!!!!!!!!!!
  //system(Form("chmod +x %s/%s/workspace_files/combinedWS.root",master_output.Data(),executable_jobname.Data()));
  TString ws_input_directory = Form("%s/%s",master_output.Data(),executable_jobname.Data());
  
  system(Form("cp -f %s %s/ws_jobfile.sh",ws_jobscript.Data(),exe.Data()));
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file=Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile=Form("%s/out/%s_%i_%i.out",dir.Data(),executable_jobname.Data(),executable_lambda,executable_lifetime);
  TString name_errfile=Form("%s/err/%s_%i_%i.err",dir.Data(),executable_jobname.Data(),executable_lambda,executable_lifetime);
  
  // Here you define the arguments for the job script:
  TString name_jscript=Form("%s %s %s %s %s %i %i",ws_jobscript.Data(),executable_jobname.Data(),input_file.Data(),executable_option.Data(),executable_name.Data(),executable_lambda,executable_lifetime);
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}

/**
   Submits toy MC jobs via bsub.
*/
void SubmitToysViaBsub(TString executable_name, TString executable_jobname, TString executable_option, int executable_seed, int executable_toys_per_job, int chosen_lambda, int chosen_lifetime) {
  
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",toy_jobscript.Data()));
  system(Form("chmod +x %s",toy_jobscriptCorr.Data()));
  system(Form("chmod +x %s/%s/workspaces/rootfiles/workspace_NPP_%iTeV_%ips.root",master_output.Data(),executable_jobname.Data(),chosen_lambda,chosen_lifetime));
  
  system(Form("cp -f %s %s/toy_jobfile.sh",toy_jobscript.Data(),exe.Data()));
  
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file = Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile = Form("%s/out/%s_%i.out",dir.Data(),executable_jobname.Data(),executable_seed);
  TString name_errfile = Form("%s/err/%s_%i.err",dir.Data(),executable_jobname.Data(),executable_seed);
  
  // Here you define the arguments for the job script:
  TString name_jscript;
  if( executable_option.Contains("correlation") )
    name_jscript = Form("%s %s %s %s %s %i %i %i %i %i", toy_jobscriptCorr.Data(), executable_jobname.Data(), input_file.Data(), executable_option.Data(), executable_name.Data(), executable_seed, executable_toys_per_job, chosen_lambda, chosen_lifetime, toy_percent_timing_corr1 );
  else
    name_jscript = Form("%s %s %s %s %s %i %i %i %i", toy_jobscript.Data(), executable_jobname.Data(), input_file.Data(), executable_option.Data(), executable_name.Data(), executable_seed, executable_toys_per_job, chosen_lambda, chosen_lifetime );
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}
