////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.h                                                          //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This header file stores all of the global information for the H->gg + DM  //
//  search with 13 TeV data in 2015. It also has all of the includes that are //
//  necessary for the analysis.                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Set True for final analysis on data:
bool doBlind = false;

// Luminosity in fb-1:
double analysis_luminosity = 20.3;

// Location of global input files:
TString masterInput = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs";

// Location of output directory:
TString masterOutput = "/afs/cern.ch/work/a/ahard/files_HDM/FullAnalysis";

// Ntuple locations:
//TString ntuple_input_background_gamma = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/list_background_gamma.txt";

// Signal cross-sections file:
//TString cross_sections_file = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/cross_sections_8TeV.txt";

// Various job scripts:
//TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
//TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";

double DMMassRangeLo = 105.0;
double DMMassRangeHi = 160.0;
