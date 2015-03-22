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

////////////////////////////////////////
// GLOBAL PARAMETERS
////////////////////////////////////////

// Set True for final analysis on data:
bool doBlind = false;

// Luminosity in fb-1:
double analysisLuminosity = 20.3;

double DMMyyRangeLo = 105.0;
double DMMyyRangeHi = 160.0;

int const nProdModes = 6;
TString sigProdModes[nProdModes];// = {"ggH","VBF","WH","ZH","ttH","bbH"};

//std::map<TString,TString> nameToSample;// = { {"ggH", "sampleName_ggH"}, 
// {"VBF", "sampleName_VBF"},
//					  {"WH", "sampleName_WH"},
//					  {"ZH", "sampleName_ZH"},
//					  {"ttH", "sampleName_ttH"},
//					  {"bbH","sampleName_bbH"} };

////////////////////////////////////////
// INPUT AND OUTPUT DIRECTORIES
////////////////////////////////////////

// Location of global input files:
TString masterInput = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs";

// Location of output directory:
TString masterOutput = "/afs/cern.ch/work/a/ahard/files_HDM/FullAnalysis";

////////////////////////////////////////
// FILE LOCATIONS
////////////////////////////////////////

// Ntuple locations:
std::map<TString,TString> nameToSample;
/*
 = { {"ggH","sampleName_ggH"},
					  {"VBF","sampleName_VBF"},
					  {"WH","~lkashif/public/for_andrew/H2yyMETAnalysis_WH/data-outputLabel/sample.root"},
					  {"ZH","~lkashif/public/for_andrew/H2yyMETAnalysis_ZH/data-outputLabel/sample.root"},
					  {"ttH","~lkashif/public/for_andrew/H2yyMETAnalysis_ttH/data-outputLabel/sample.root"},
					  {"shxx_gg_ms100_mx100","~lkashif/public/for_andrew/H2yyMETAnalysis_shxx_gg_ms100_mx100/data-outputLabel/sample.root"},
					  {"shxx_gg_ms100_mx500","~lkashif/public/for_andrew/H2yyMETAnalysis_shxx_gg_ms100_mx500/data-outputLabel/sample.root"}
};
*/


////////////////////////////////////////
// SCRIPT LOCATIONS
////////////////////////////////////////

// Various job scripts:
//TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
//TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";
