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
//  MasterOption - Note: Each can be followed by the suffix "New"             //
//    - MassPoints                                                            //
//    - SigParam                                                              //
//    - BkgModel                                                              //
//    - Workspace                                                             //
//    - ToyMC                                                                 //
//    - CalcCLs                                                               //
//                                                                            //
//  Need to rethink the DMSigParam handling of the RooDataSet. Maybe we       //
//  should just hand it a RooDataSet?                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMaster.h"

/**
   This is the main DMMaster method:
*/
int main( int argc, char **argv ) {
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <MasterJobName> <MasterOption>\n\n",argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterJobName = argv[1];
  TString masterOption = argv[2];
    
  //--------------------------------------//
  // Analysis tools to be initialized depending on options:
  DMMassPoints *mp;// May need an array, unlike sig param. 
  DMSigParam *sp;
  RooRealVar *m_yy = new RooRealVar("m_yy","m_yy",DMMyyRangeLo, DMMyyRangeHi);
  RooCategory *cate = new RooCategory(Form("categories_%s", cateScheme.Data()),
				      Form("categories_%s", cateScheme.Data()));
  // Loop over categories to define categories:
  for (int i_c = 0; i_c < selector->getNCategories(newCateScheme); i_c++) {
    newCategories->defineType(Form("%s_%d",newCateScheme.Data(),i_c));
    //newCategories->setRange(Form("rangeName_",i_b,i_r),Form("%s_%d",cateScheme.Data(),i_c));
  }
    
  //--------------------------------------//
  // Step 1: Make or load mass points:
  if (masterOption.Contains("MassPoints")) {
    cout << "DMMaster: Step 1 - Make mass points." << endl;
    
    // NOTE: WANT TO LOOP OVER FILES OR SAMPLES...
    for (int i_s = 0; i_s < (int)sampleNames.size(); i_s++) {
     
      if (masterOption.Contains("MassPointsNew")) {
	mp = new DMMassPoints(masterJobName, file, cateScheme, 
			      "New", m_yy, cate);
      }
      else {
	mp = new DMMassPoints(masterJobName, file, cateScheme, 
			      "FromFile", m_yy, cate);
      }
    }
  }
  
  //--------------------------------------//
  // Step 2: Make or load the signal parameterization:
  if (masterOption.Contains("SigParam")) {
    cout << "DMMaster: Step 2 - Make signal parameterization." << endl;
  
    if (masterOption.Contains("SigParamNew")) {
      sp = new DMSigParam(masterJobName, cateScheme, "New", m_yy, cate);
    }
    else {
      sp = new DMSigParam(masterJobName, cateScheme, "FromFile", m_yy, cate);
    }
  }
  
  //--------------------------------------//
  // Step 3: Make the workspace:
  if (masterOption.Contains()) {
    cout << "DMMaster: Step 3 - Make the workspace for fits." << endl;
    
  }
  
  return 0;
}

/*
//   Compiles the executable required by the job options.

void MakeExe(TString exename) {
  
  // recompile all executables before running...
  system(Form("rm bin/%s",exename.Data()));
  system(Form("make bin/%s",exename.Data()));
}


//   Submits the workspace jobs via bsub.

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


//   Submits toy MC jobs via bsub.

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
