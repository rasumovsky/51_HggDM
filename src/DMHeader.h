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

/// c++ headers
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cassert>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
#include <map>
#include <cstdlib>
#include <cmath>

/// ROOT headers
#include <TROOT.h>
#include <TObject.h>
#include <TSelector.h>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TCollection.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>
#include <TH2.h>
#include <THStack.h>
#include <TBox.h>
#include <TLegend.h>
#include <TPad.h>
#include <TF1.h>
#include <TProfile.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TLine.h>
#include <TFrame.h>
#include <TAxis.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TVersionCheck.h>
#include <TGenPhaseSpace.h>
#include <TPaveStats.h>
#include <TPolyLine3D.h>
#include <TH2Poly.h>
#include <TGaxis.h>

/// ROOT math headers
#include <Math/QuantFuncMathCore.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <Math/Factory.h>

/// TMVA headers
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>

// Set True for final analysis on data:
bool doBlind = false;

// Luminosity in fb-1:
double analysis_luminosity = 20.3;

// Location of global input files:
TString master_input = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs";

// Location of output directory:
TString master_output = "/afs/cern.ch/work/a/ahard/files_HDM/FullAnalysis";

// Ntuple locations:
//TString ntuple_input_background_gamma = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/list_background_gamma.txt";

// Signal cross-sections file:
//TString cross_sections_file = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs/cross_sections_8TeV.txt";

// Various job scripts:
//TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
//TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";

