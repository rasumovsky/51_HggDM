////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonHead.h"

using namespace std;

/**
   ALL Important global values or file locations for the analysis should be
   stored in this header file. 
 */

bool doBlind = false;

double analysis_luminosity = 20.3;

TString master_input = "/afs/cern.ch/work/a/ahard/files_NPP/GlobalInputs";
TString master_output = "/afs/cern.ch/work/a/ahard/files_NPP/FullAnalysis";

// Ntuple locations:
TString ntuple_input_background_gamma = "/afs/cern.ch/work/a/ahard/files_NPP/GlobalInputs/list_background_gamma.txt";

// Signal cross-sections file:
TString cross_sections_file = "/afs/cern.ch/work/a/ahard/files_NPP/GlobalInputs/cross_sections_8TeV.txt";

// Various bsub job scripts:
TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/41_NPP_WS/scripts/ws_jobfile.sh";
TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/41_NPP_WS/scripts/toy_jobfile.sh";

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// MakeExe:

void MakeExe( TString exename )
{
  // recompile all executables before running...
  system(Form("rm bin/%s",exename.Data()));
  system(Form("make bin/%s",exename.Data()));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitViaBsub:

void SubmitWSViaBsub( TString executable_name, TString executable_jobname, TString executable_option, int executable_lambda, int executable_lifetime )
{
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitToysViaBsub:

void SubmitToysViaBsub( TString executable_name, TString executable_jobname, TString executable_option, int executable_seed, int executable_toys_per_job, int chosen_lambda, int chosen_lifetime )
{
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
  
  if( executable_option.Contains("correlation") )
  {
    //system(Form("cp -f %s %s/toy_jobfile.sh",toy_jobscriptCorr.Data(),exe.Data()));
    system(Form("cp -f %s %s/toy_jobfileCorrelation.sh",toy_jobscriptCorr.Data(),exe.Data()));
  }
  else
  {
    system(Form("cp -f %s %s/toy_jobfile.sh",toy_jobscript.Data(),exe.Data()));
    //system(Form("cp -f %s %s/toy_jobfileCorrelation.sh",toy_jobscript.Data(),exe.Data()));
  }
  
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitCLViaBsub:

void SubmitCLViaBsub( TString executable_name, TString executable_jobname, int chosen_lambda, int chosen_lifetime )
{
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s_CL",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",CL_jobscript.Data()));
  system(Form("chmod +x %s/%s/workspaces/rootfiles/workspace_NPP_%iTeV_%ips.root",master_output.Data(),executable_jobname.Data(),chosen_lambda,chosen_lifetime));
  
  system(Form("cp -f %s %s/CL_jobfile.sh",CL_jobscript.Data(),exe.Data()));
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file = Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile = Form("%s/out/%s_%iTeV_%ips.out", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  TString name_errfile = Form("%s/err/%s_%iTeV_%ips.err", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  
  // Here you define the arguments for the job script:
  TString name_jscript = Form("%s %s %s %s %i %i", CL_jobscript.Data(), executable_jobname.Data(), input_file.Data(), executable_name.Data(), chosen_lambda, chosen_lifetime );
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitMuLimitViaBsub:

void SubmitMuLimitViaBsub( TString executable_name, TString executable_jobname, int chosen_lambda, int chosen_lifetime, TString limit_joboption )
{
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s_MuLimit",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",MuLimit_jobscript.Data()));
  system(Form("chmod +x %s/%s/workspaces/rootfiles/workspace_NPP_%iTeV_%ips.root",master_output.Data(),executable_jobname.Data(),chosen_lambda,chosen_lifetime));
  
  system(Form("cp -f %s %s/MuLimit_jobfile.sh",MuLimit_jobscript.Data(),exe.Data()));
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file = Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile = Form("%s/out/%s_%iTeV_%ips.out", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  TString name_errfile = Form("%s/err/%s_%iTeV_%ips.err", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  
  // Here you define the arguments for the job script:
  TString name_jscript = Form("%s %s %s %s %i %i %s", MuLimit_jobscript.Data(), executable_jobname.Data(), input_file.Data(), executable_name.Data(), chosen_lambda, chosen_lifetime, limit_joboption.Data() );
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitP0ViaBsub:

void SubmitP0ViaBsub( TString executable_name, TString executable_jobname, int chosen_lambda, int chosen_lifetime )
{
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s_p0",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",p0_jobscript.Data()));
  system(Form("chmod +x %s/%s/workspaces/rootfiles/workspace_NPP_%iTeV_%ips.root",master_output.Data(),executable_jobname.Data(),chosen_lambda,chosen_lifetime));
  
  system(Form("cp -f %s %s/p0_jobfile.sh",p0_jobscript.Data(),exe.Data()));
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file = Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile = Form("%s/out/%s_%iTeV_%ips.out", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  TString name_errfile = Form("%s/err/%s_%iTeV_%ips.err", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  
  // Here you define the arguments for the job script:
  TString name_jscript = Form("%s %s %s %s %i %i", p0_jobscript.Data(), executable_jobname.Data(), input_file.Data(), executable_name.Data(), chosen_lambda, chosen_lifetime );
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SubmitGlobP0ViaBsub:

void SubmitGlobP0ViaBsub( TString executable_name, TString executable_jobname, int chosen_lambda, int chosen_lifetime, int maximum_toys )
{
  // Make directories for job info:
  TString dir = Form("/afs/cern.ch/work/a/ahard/jobfiles/%s_globp0",executable_jobname.Data());
  TString out = Form("%s/out",dir.Data());
  TString err = Form("%s/err",dir.Data());
  TString exe = Form("%s/exe",dir.Data());
  system(Form("mkdir -vp %s",out.Data()));
  system(Form("mkdir -vp %s",err.Data()));
  system(Form("mkdir -vp %s",exe.Data()));
  
  // create .tar file with everything:
  system(Form("tar zcf Cocoon.tar bin/%s",executable_name.Data()));
  system(Form("chmod +x %s",globp0_jobscript.Data()));
  system(Form("chmod +x %s/%s/workspaces/rootfiles/workspace_NPP_%iTeV_%ips.root",master_output.Data(),executable_jobname.Data(),chosen_lambda,chosen_lifetime));
  system(Form("chmod +x %s/%s/toyData_260TeV_1291ps/*.root",master_output.Data(),executable_jobname.Data()));
  
  system(Form("cp -f %s %s/globp0_jobfile.sh",globp0_jobscript.Data(),exe.Data()));
  system(Form("mv Cocoon.tar %s",exe.Data()));
  TString input_file = Form("%s/Cocoon.tar",exe.Data());
  TString name_outfile = Form("%s/out/%s_%iTeV_%ips.out", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  TString name_errfile = Form("%s/err/%s_%iTeV_%ips.err", dir.Data(), executable_jobname.Data(), chosen_lambda, chosen_lifetime);
  
  // Here you define the arguments for the job script:
  TString name_jscript = Form("%s %s %s %s %i %i %i", globp0_jobscript.Data(), executable_jobname.Data(), input_file.Data(), executable_name.Data(), chosen_lambda, chosen_lifetime, maximum_toys );
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s",name_outfile.Data(),name_errfile.Data(),name_jscript.Data()));
}
