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

TString cateToBkgFunc(TString category) {
  TString result = "";
  result = "Exppol01";
  //Possibilities are "BernO1",... "BernO6", "ExppolO1",... "ExppolO6"
  return result;
}

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

TString nameToRootFile(TString name) {
  TString result = "";
  if (name.EqualTo("ggH")) result = "sampleName_ggH";
  else if (name.EqualTo("VBF")) result = "sampleName_VBF";
  else if (name.EqualTo("WH")) result = "~lkashif/public/for_andrew/H2yyMETAnalysis_WH/data-outputLabel/sample.root";
  else if (name.EqualTo("ZH")) result = "~lkashif/public/for_andrew/H2yyMETAnalysis_ZH/data-outputLabel/sample.root";
  else if (name.EqualTo("ttH")) result = "~lkashif/public/for_andrew/H2yyMETAnalysis_ttH/data-outputLabel/sample.root";
  else if (name.EqualTo("bbH")) result = "sampleName_bbH";
  else std::cout << "nameToRootFile: Error! No corresponding file" << std::endl;
  return result;
}

TString nameToFileList(TString name) {
  TString result = Form("%s/FileLists/",masterInput.Data());
  if (name.EqualTo("ggH")) result += "fileList_ggH.txt";
  else if (name.EqualTo("VBF")) result += "fileList_VBF.txt";
  else if (name.EqualTo("WH")) result += "fileList_WH.txt";
  else if (name.EqualTo("ZH")) result += "fileList_ZH.txt";
  else if (name.EqualTo("ttH")) result += "fileList_ttH.txt";
  else if (name.EqualTo("bbH")) result += "fileList_bbH.txt";
  else std::cout << "nameToRootFile: Error! No corresponding file" << std::endl;
  return result;
}

////////////////////////////////////////
//          SCRIPT LOCATIONS          //
////////////////////////////////////////

//TString ws_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/ws_jobfile.sh";
//TString toy_jobscript = "/afs/cern.ch/user/a/ahard/work_directory/analysis/51_HDM/scripts/toy_jobfile.sh";
