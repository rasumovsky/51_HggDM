////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMheckjobs.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/04/2015                                                          //
//                                                                            //
//  This class checks to see whether jobs of a particular type have finished. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMCheckJobs.h"

DMCheckJobs::DMCheckJobs(TString newJobName) {
  
  jobName = newJobName;
    
  // Save names and values of failed jobs:
  int nfailed_jobs = 0;
  vector<TString> incomplete_files; incomplete_files.clear();
  vector<int> incomplete_lambda; incomplete_lambda.clear();
  vector<int> incomplete_lifetime; incomplete_lifetime.clear();
  
  TString output_name;
  if( type == "WS" ) output_name = Form("%s/%s/workspaces/failed_ws_points.txt",master_output.Data(),jobname.Data());
  if( type == "CL" ) output_name = Form("%s/%s/CL_limits/failed_cl_points.txt",master_output.Data(),jobname.Data());
  if( type == "MU" ) output_name = Form("%s/%s/limits_mu/failed_mu_points.txt",master_output.Data(),jobname.Data());
  if( type == "p0" ) output_name = Form("%s/%s/p0_values/failed_p0_points.txt",master_output.Data(),jobname.Data());
  ofstream failed_points(output_name);
  
  //--------------------------------------//
  // Then loop over submissions to see whether output files exist:
  for( int i_l = 0; i_l < nlambda_values; i_l++ )
  {
    //if( type == "MU" && ( lambda_values[i_l] != specified_lambda ) ) continue;
    for( int i_t = 0; i_t < nlifetime_values; i_t++ )
    {
      int current_lambda = lambda_values[i_l];
      int current_lifetime = ps_lifetimes[i_t];
      TString file_name;
      if( type == "WS" ) file_name = Form( "workspace_NPP_%iTeV_%ips.root", current_lambda, current_lifetime );
      if( type == "CL" ) file_name = Form( "CL_values_%iTeV_%ips.txt", current_lambda, current_lifetime );
      if( type == "MU" ) file_name = Form( "text_CLs_%iTeV_%ips.txt", current_lambda, current_lifetime );
      if( type == "p0" ) file_name = Form( "p0_values_%iTeV_%ips.txt", current_lambda, current_lifetime );
      TString full_name;
      if( type == "WS" ) full_name = Form("%s/%s/workspaces/rootfiles/%s",master_output.Data(),jobname.Data(),file_name.Data());
      if( type == "CL" ) full_name = Form("%s/%s/CL_limits/single_files/%s",master_output.Data(),jobname.Data(),file_name.Data());
      if( type == "MU" ) full_name = Form("%s/%s/limits_mu/single_files/%s",master_output.Data(),jobname.Data(),file_name.Data());
      if( type == "p0" ) full_name = Form("%s/%s/p0_values/single_files/%s",master_output.Data(),jobname.Data(),file_name.Data());
      ifstream test_file(full_name);
      if( !test_file )
      {
	incomplete_files.push_back(file_name);
	incomplete_lambda.push_back(current_lambda);
	incomplete_lifetime.push_back(current_lifetime);
	failed_points << current_lambda << " " << current_lifetime << endl;
	nfailed_jobs++;
      }
    }
  }
  
  //--------------------------------------//
  // Print information on failed jobs:
  cout << "Failed to make the following files ( " << nfailed_jobs << " total )" << endl;
  for( int i_f = 0; i_f < nfailed_jobs; i_f++ )
  {
    cout << "  " << incomplete_lambda[i_f] << " TeV, \t" << incomplete_lifetime[i_f] << " ps, \t" << incomplete_files[i_f] << endl;
  }
  failed_points.close();
  cout << "End of NPPV1_checkjobs.cc" << endl;
  return 0;
}
