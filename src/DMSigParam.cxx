/////////////////////////////////////////
//                                     //
//  spinv1_parameterization.cc         //
//                                     //
//  Author: Andrew Hard                //
//  Date: 21/01/2014                   //
//                                     //
//  Description:                       //
//  This code implements the global    //
//  signal parameterization for the    //
//  the 2013 Higgs -> gamma gamma spin //
//  publication.                       //
//                                     //
//  Updated with the 11 categories (10 //
//  10 Cos(theta*) categories and a PT //
//  category.                          //
//                                     //
//  Inputs:                            //
//  1. output name                     //
//  2. list of root files              //
//                                     //
/////////////////////////////////////////

#include "spinv1_parameterization.hh"

int main (int argc, char **argv)
{
  // Root macros:
  SetAtlasStyle();
  
  // Check arguments:
  if( argc < 4 ){ printf("\nUsage: %s <jobname> <input list> <energy>\n\n",argv[0]); exit(0); }
  jobname = argv[1];
  sample_list = argv[2];
  energy = atoi(argv[3]);
  
  // SetNumberCategories:
  ncat = 11;
  
  // Print input file:
  cout << "Reading signal MC from: " << endl;
  system(Form("cat %s",sample_list.Data()));
  cout << " " << endl;
  
  // Create output directories:
  output_directory = Form("%s/%s/parameterization_%iTeV",master_output.Data(),jobname.Data(),energy);
  system(Form("mkdir -vp %s",output_directory.Data()));
  system(Form("mkdir -vp %s/Plots",output_directory.Data())); 
  system(Form("mkdir -vp %s/all",output_directory.Data()));
  system(Form("mkdir -vp %s/ggf",output_directory.Data()));
  system(Form("mkdir -vp %s/vbf",output_directory.Data()));
  system(Form("mkdir -vp %s/wh",output_directory.Data()));
  system(Form("mkdir -vp %s/zh",output_directory.Data()));
  system(Form("mkdir -vp %s/tth",output_directory.Data()));
  system(Form("mkdir -vp %s/FinalSignal",output_directory.Data()));
  
  //--------------------------------------//
  // Make TChain and find number of included events:
  chain = MakeChain("Hgg_tree", sample_list, "badfile");
  p = new ntup_SpinPub(chain);
  numev = p->fChain->GetEntries();
  
  //--------------------------------------//
  // Create TTrees for each mass point:
  for(int icat = 0; icat <= ncat; icat++ )
  {
    for( int t = 0; t < nmasspoint; t++ )
    {
      TString tname = Form("t_%d_%d", icat, signalmass[t]);
      tree[icat][t] = new TTree(tname,tname);
      tree[icat][t]->Branch("mass", &mass[icat], "mass/D");
      tree[icat][t]->Branch("weight", &weight[icat], "weight/D");
    }
  }
  
  //--------------------------------------//
  // Parameterization of SM signal in CTS bins:
  DoSignalParameterization();
  
  //--------------------------------------//
  // Interpolation of SM signal:
  cout << "Beginning interpolation loop." << endl;
  for( int i_c = 0; i_c <= ncat; i_c++ )
    InterpolateSignal( i_c );
  
  //--------------------------------------//
  // Then get the signal yields at the specific mass:
  GetYieldsAtSpecificMass(); 
  
  return 0;
}
