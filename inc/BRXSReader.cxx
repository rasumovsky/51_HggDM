////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BRXSReader.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/03/2015                                                          //
//                                                                            //
//  This class is designed to load the cross-sections and branching-ratios    //
//  for the SM Higgs boson at a variety of masses, for sqrt(s) = 13 TeV.      //
//                                                                            //
//  13 TeV Cross-section values from:                                         //
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
//                                                                            //
//  Branching-ratio values from:                                              //
// https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3    //
//                                                                            //
//  After initializing the tool, cross-sections and branching ratios can be   //
//  accessed using the methods below.                                         //
//                                                                            //
//  double getXS(int mass, TString production, TString value);                //
//    production = "ggH", "VBF", "WH", "ZH", "ttH", "bbH"                     //
//    type = "XS" for nominal cross-section value                             //
//           "+QCD"   for the XS + QCD uncertainty in %                       //
//           "-QCD"   for the XS - QCD uncertainty in %                       //
//           "+PDF"   for the XS + PDF uncertainty in %                       //
//           "-PDF"   for the XS - PDF uncertainty in %                       //
//                                                                            //
//  double getBR(int mass, TString decay, TString value);                     //
//    decay = "gg", "gammagamma", "Zgamma", "WW", "ZZ",                       //
//            "bb", "tauttau", "mumu", "cc", "ss", "tt"                       //
//    type = "+ERR"   for the value of the BR + total error in %              //
//           "-ERR"   for the value of the BR - total error in %              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "BRXSReader.h"

/**
   Initializes the tool and loads XS, BR values from files. 
*/
BRXSReader::BRXSReader(TString inputDirectory) {
  
  directory = inputDirectory;
  valuesXS.clear();
  valuesBR.clear();
  
  // Open cross-section files and store values.
  loadXS("ggH");
  loadXS("VBF");
  loadXS("WH");
  loadXS("ZH");
  loadXS("ttH");
  loadXS("bbH");
  
  // Open BR files and store values.
  loadBR("2bosons");
  loadBR("2fermions");
}

/**
   Returns the branching-ratio, or related uncertainties, for a decay mode
   at a particular Higgs mass. 
*/
float BRXSReader::getBR(double mass, TString decay, TString value) {
  TString currKey = getMapKey(mass, decay, value);
  if (hasKey(currKey,"BR")) {
    return valuesBR[currKey];
  }
  else {
    std::cout << "BRXSReader Error! No match for " << mass << " and " 
	      << decay << " found." << std::endl;
    return 0;
  }
}

/**
   Returns the cross-section, or related uncertainties, for a production 
   process at a particular mass.
*/
float BRXSReader::getXS(double mass, TString production, TString value) {
  TString currKey = getMapKey(mass, production, value);
  if (hasKey(currKey,"XS")) {
    return valuesXS[currKey];
  }
  else {
    std::cout << "BRXSReader Error! No match for " << mass << " and " 
	      << production << " found." << std::endl;
    return 0;
  }
}

/**
   Print the specified branching-ratio values.
*/
void BRXSReader::printBR(double mass, TString decay) {
  if (hasKey(getMapKey(mass, decay, "BR"),"BR")) {
    std::cout << "Branching-ratio data for " << decay << std::endl;
    std::cout << "  BR=" << valuesXS[getMapKey(mass, decay, "BR")];
    std::cout << "  +ERR=" << valuesXS[getMapKey(mass, decay, "+ERR")];
    std::cout << "  -ERR=" << valuesXS[getMapKey(mass, decay, "-ERR")];
    std::cout << std::endl;
  }
  else {
    std::cout <<"BRXSReader: Error! No matching entry." << std::endl;
  }
}

/**
   Print the specified cross-section values.
*/
void BRXSReader::printXS(double mass, TString production) {
  if (hasKey(getMapKey(mass, production, "XS"),"XS")) {
    std::cout << "Cross-section data for " << production << std::endl;
    std::cout << "  XS=" << valuesXS[getMapKey(mass, production, "XS")];
    std::cout << "  +QCD=" << valuesXS[getMapKey(mass, production, "+QCD")];
    std::cout << "  -QCD=" << valuesXS[getMapKey(mass, production, "-QCD")];
    std::cout << "  +PDF=" << valuesXS[getMapKey(mass, production, "+PDF")];
    std::cout << "  -PDF=" << valuesXS[getMapKey(mass, production, "-PDF")];
    std::cout << std::endl;
  }
  else {
    std::cout <<"BRXSReader: Error! No matching entry." << std::endl;
  }
}

/**
   Load the branching ratios from input files.
*/
void BRXSReader::loadBR(TString decayClass) {
  
  TString fileName = Form("%s/BR_%s.txt", directory.Data(), decayClass.Data());
  ifstream currFile(fileName);
  if (currFile.is_open()) {
    
    double currMass;
    float currIn[19];
    
    while (!currFile.eof()) {
      currFile >> currMass >> currIn[1] >> currIn[2] >> currIn[3] >> currIn[4]
	       >> currIn[5] >> currIn[6] >> currIn[7] >> currIn[8] >> currIn[9]
	       >> currIn[10] >> currIn[11] >> currIn[12] >> currIn[13]
	       >> currIn[14] >> currIn[15] >> currIn[16] >> currIn[17]
	       >> currIn[18];
           
      if (decayClass.Contains("2bosons")) {
	valuesBR[getMapKey(currMass, "gg", "BR")] = currIn[1];
	valuesBR[getMapKey(currMass, "gg", "+ERR")] = currIn[2];
	valuesBR[getMapKey(currMass, "gg", "-ERR")] = currIn[3];
	valuesBR[getMapKey(currMass, "gammagamma", "BR")] = currIn[4];
	valuesBR[getMapKey(currMass, "gammagamma", "+ERR")] = currIn[5];
	valuesBR[getMapKey(currMass, "gammagamma", "-ERR")] = currIn[6];
	valuesBR[getMapKey(currMass, "Zgamma", "BR")] = currIn[7];
	valuesBR[getMapKey(currMass, "Zgamma", "+ERR")] = currIn[8];
	valuesBR[getMapKey(currMass, "Zgamma", "-ERR")] = currIn[9];
	valuesBR[getMapKey(currMass, "WW", "BR")] = currIn[10];
	valuesBR[getMapKey(currMass, "WW", "+ERR")] = currIn[11];
	valuesBR[getMapKey(currMass, "WW", "-ERR")] = currIn[12];
	valuesBR[getMapKey(currMass, "ZZ", "BR")] = currIn[13];
	valuesBR[getMapKey(currMass, "ZZ", "+ERR")] = currIn[14];
	valuesBR[getMapKey(currMass, "ZZ", "-ERR")] = currIn[15];
      }
      
      else if (decayClass.Contains("2fermions")) {
	valuesBR[getMapKey(currMass, "bb", "BR")] = currIn[1];
	valuesBR[getMapKey(currMass, "bb", "+ERR")] = currIn[2];
	valuesBR[getMapKey(currMass, "bb", "-ERR")] = currIn[3];
	valuesBR[getMapKey(currMass, "tautau", "BR")] = currIn[4];
	valuesBR[getMapKey(currMass, "tautau", "+ERR")] = currIn[5];
	valuesBR[getMapKey(currMass, "tautau", "-ERR")] = currIn[6];
	valuesBR[getMapKey(currMass, "mumu", "BR")] = currIn[7];
	valuesBR[getMapKey(currMass, "mumu", "+ERR")] = currIn[8];
	valuesBR[getMapKey(currMass, "mumu", "-ERR")] = currIn[9];
	valuesBR[getMapKey(currMass, "cc", "BR")] = currIn[10];
	valuesBR[getMapKey(currMass, "cc", "+ERR")] = currIn[11];
	valuesBR[getMapKey(currMass, "cc", "-ERR")] = currIn[12];
	valuesBR[getMapKey(currMass, "ss", "BR")] = currIn[13];
	valuesBR[getMapKey(currMass, "ss", "+ERR")] = currIn[14];
	valuesBR[getMapKey(currMass, "ss", "-ERR")] = currIn[15];
	valuesBR[getMapKey(currMass, "tt", "BR")] = currIn[16];
	valuesBR[getMapKey(currMass, "tt", "+ERR")] = currIn[17];
	valuesBR[getMapKey(currMass, "tt", "-ERR")] = currIn[18];
      }
      
      else {
	std::cout << "BRXSReader: No decay class provided!" << std::endl;
      }
    }
  }
  currFile.close();
}

/**
   Load the cross-sections from input files.
*/
void BRXSReader::loadXS(TString production) {
  
  TString fileName = Form("%s/XS_%s.txt", directory.Data(), production.Data());
  ifstream currFile(fileName);
  if (currFile.is_open()) {
    
    double currMass;
    float currIn[6];
    
    while (!currFile.eof()) {
      currFile >> currMass >> currIn[1] >> currIn[2] >> currIn[3] >> currIn[4]
	       >> currIn[5];
      
      valuesXS[getMapKey(currMass, production, "XS")] = currIn[1];
      valuesXS[getMapKey(currMass, production, "+QCD")] = currIn[2];
      valuesXS[getMapKey(currMass, production, "-QCD")] = currIn[3];
      valuesXS[getMapKey(currMass, production, "+PDF")] = currIn[4];
      valuesXS[getMapKey(currMass, production, "-PDF")] = currIn[5];
    }
  }
  currFile.close();
}

/**
   Convert variables to the key for the map.
*/
TString BRXSReader::getMapKey(double mass, TString type, TString value) {
  int massInt = (int)(100*mass);
  TString key = Form("%d_%s_%s", massInt, type.Data(), value.Data());
  return key;
}

/**
   Check whether the specified key is contained in the specified map.
*/
bool BRXSReader::hasKey(TString key, TString mapType) {
  if (mapType.EqualTo("XS")) {
    return (valuesXS.count(key) > 0);
  }
  else if (mapType.EqualTo("BR")) {
    return (valuesBR.count(key) > 0);
  }
  else {
    std::cout << "BRXSReader: Error! Improper mapType argument." << std::endl;
    return false;
  }
}
