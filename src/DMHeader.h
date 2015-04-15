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
//  13/04/2015 NOTE                                                           //
//    Maybe we should make this a struct!                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////
//         GLOBAL PARAMETERS          //
////////////////////////////////////////

// Set True for final analysis on data:
bool doBlind = false;

// Luminosity in fb-1:
double analysisLuminosity = 20.3;

double DMMyyRangeLo = 105.0;
double DMMyyRangeHi = 160.0;

int const nProdModes = 6;
TString sigProdModes[nProdModes] = {"ggH","VBF","WH","ZH","ttH","bbH"};

////////////////////////////////////////
//    INPUT AND OUTPUT DIRECTORIES    //
////////////////////////////////////////

// Location of global input files:
TString masterInput = "/afs/cern.ch/work/a/ahard/files_HDM/GlobalInputs";
// Location of output directory:
TString masterOutput = "/afs/cern.ch/work/a/ahard/files_HDM/FullAnalysis";

////////////////////////////////////////
//           FILE LOCATIONS           //
////////////////////////////////////////

TString nameToFileList(TString name) {
  TString result = Form("%s/FileLists/",masterInput.Data());
  if (name.EqualTo("ggH")) {
    result += "list_H2yyMETAnalysis_ggH.txt";
  }
  else if (name.EqualTo("VBF")) {
    result += "list_H2yyMETAnalysis_VBF.txt";
  }
  else if (name.EqualTo("WH")) {
    result += "list_H2yyMETAnalysis_WH.txt";
  }
  else if (name.EqualTo("ZH")) {
    result += "list_H2yyMETAnalysis_ZH.txt";
  }
  else if (name.EqualTo("ttH")) {
    result += "list_H2yyMETAnalysis_ttH.txt";
  }
  else if (name.EqualTo("gg_gjet")) {
    result += "list_H2yyMETAnalysis_gg_gjet.txt";
  }
  else if (name.EqualTo("shxx_gg_ms100_mx100")) {
    result += "list_H2yyMETAnalysis_shxx_gg_ms100_mx100.txt";
  }
  else if (name.EqualTo("shxx_gg_ms100_mx500")) {
    result += "list_H2yyMETAnalysis_shxx_gg_ms100_mx500.txt";
  }
  else if (name.EqualTo("zphxx_gg_mzp100_mx100")) {
    result += "list_H2yyMETAnalysis_zphxx_gg_mzp100_mx100.txt";
  }
  else {
    std::cout << "nameToRootFile: Error! No corresponding file" << std::endl;
  }
  return result;
}

////////////////////////////////////////
//          BACKGROUND PDFS           //
////////////////////////////////////////

TString cateToBkgFunc(TString category) {
  TString result = "";
  result = "Exppol01";
  //Possibilities are "BernO1",... "BernO6", "ExppolO1",... "ExppolO6"
  return result;
}

////////////////////////////////////////
//          SCRIPT LOCATIONS          //
////////////////////////////////////////

//TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
//TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";
