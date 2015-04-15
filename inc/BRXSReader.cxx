////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BRXSReader.cxx                                                      //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/03/2015                                                          //
//                                                                            //
//  This class is designed to load the cross-sections and branching-ratios    //
//  for the SM Higgs boson at a variety of masses, for sqrt(s) = 13 TeV. All  //
//  cross sections are in units of pb.                                        //
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
//  double getSMXS(int mass, TString production, TString value);              //
//    production = "ggH", "VBF", "WH", "ZH", "ttH", "bbH"                     //
//    type = "XS" for nominal cross-section value                             //
//           "+QCD"   for the XS + QCD uncertainty in %                       //
//           "-QCD"   for the XS - QCD uncertainty in %                       //
//           "+PDF"   for the XS + PDF uncertainty in %                       //
//           "-PDF"   for the XS - PDF uncertainty in %                       //
//                                                                            //
//  double getSMBR(int mass, TString decay, TString value);                   //
//    decay = "gg", "gammagamma", "Zgamma", "WW", "ZZ",                       //
//            "bb", "tauttau", "mumu", "cc", "ss", "tt"                       //
//    type = "+ERR"   for the value of the BR + total error in %              //
//           "-ERR"   for the value of the BR - total error in %              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "BRXSReader.h"

/**
   Initializes the tool and loads XS, BR values from files. 
   @param inputDirectory - The location of the XS and BR tables from the LHCWG
   @returns void.
*/
BRXSReader::BRXSReader(TString inputDirectory) {
  
  directory = inputDirectory;
  valuesXS.clear();
  valuesBR.clear();
  
  // Open SM cross-section files and store values.
  loadSMXS("ggH");
  loadSMXS("VBF");
  loadSMXS("WH");
  loadSMXS("ZH");
  loadSMXS("ttH");
  loadSMXS("bbH");
  
  // Open the DM cross-section files and store the values.
  loadDMXS();

  // Open BR files and store values.
  loadSMBR("2bosons");
  loadSMBR("2fermions");
}

/**
   Returns the Standard Model branching-ratio, or related uncertainties, for a
   decay mode at a particular Higgs boson mass. 
   @param mass - The signal mass of interest.
   @param decay - The decay name. 
   @param value - The value (BR, +ERR, -ERR).
   @returns - The value of the branching-ratio.
*/
float BRXSReader::getSMBR(double mass, TString decay, TString value) {
  TString currKey = getSMMapKey(mass, decay, value);
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
   Returns the cross-section, or related uncertainties, for a Dark Matter
   production process at a particular intermediate and fermion mass.
   @param massIntermediate - The scalar or Zprime mass.
   @param massFermion - The dark matter particle mass.
   @param type - The production type (shxx_gg, zphxx_gg).
   @param value - The value (XS, +ERR, -ERR).
   @returns - The value of the production cross-section.
*/
float BRXSReader::getDMXSBR(int massIntermediate, int massFermion, TString type,
			    TString value) {
  TString currKey = getDMMapKey(massIntermediate, massFermion, type, value);
  if (hasKey(currKey,"XS")) {
    return valuesXS[currKey];
  }
  else {
    std::cout << "BRXSReader Error! No match for " << type << " "
	      << massIntermediate << " and " << massFermion
	      << " found." << std::endl;
    return 0;
  }
}

/**
   Returns the cross-section, or related uncertainties, for a Standard Model
   production process at a particular Higgs boson mass.
   @param mass - The signal mass of interest.
   @param production - The production mode name. 
   @param value - The value (XS, +QCD, -QCD, +PDF, -PDF).
   @returns - The value of the production cross-section.
*/
float BRXSReader::getSMXS(double mass, TString production, TString value) {
  TString currKey = getSMMapKey(mass, production, value);
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
   @param mass - The signal mass of interest.
   @param decay - The decay name. 
   @returns - void after printing BR data.
*/
void BRXSReader::printSMBR(double mass, TString decay) {
  if (hasKey(getSMMapKey(mass, decay, "BR"),"BR")) {
    std::cout << "Branching-ratio data for " << decay << std::endl;
    std::cout << "  BR=" << valuesXS[getSMMapKey(mass, decay, "BR")];
    std::cout << "  +ERR=" << valuesXS[getSMMapKey(mass, decay, "+ERR")];
    std::cout << "  -ERR=" << valuesXS[getSMMapKey(mass, decay, "-ERR")];
    std::cout << std::endl;
  }
  else {
    std::cout <<"BRXSReader: Error! No matching entry." << std::endl;
  }
}

/**
   Print the specified cross-section values.
   @param mass - The signal mass of interest.
   @param production - The production mode name. 
   @returns - void after printing XS data.
*/
void BRXSReader::printSMXS(double mass, TString production) {
  if (hasKey(getSMMapKey(mass, production, "XS"),"XS")) {
    std::cout << "Cross-section data for " << production << std::endl;
    std::cout << "  XS=" << valuesXS[getSMMapKey(mass, production, "XS")];
    std::cout << "  +QCD=" << valuesXS[getSMMapKey(mass, production, "+QCD")];
    std::cout << "  -QCD=" << valuesXS[getSMMapKey(mass, production, "-QCD")];
    std::cout << "  +PDF=" << valuesXS[getSMMapKey(mass, production, "+PDF")];
    std::cout << "  -PDF=" << valuesXS[getSMMapKey(mass, production, "-PDF")];
    std::cout << std::endl;
  }
  else {
    std::cout <<"BRXSReader: Error! No matching entry." << std::endl;
  }
}

/**
   Load the Standard Model branching ratios from input files.
   @param decayClass - The class of decays ("2bosons", "2fermions").
   @returns void
*/
void BRXSReader::loadSMBR(TString decayClass) {
  
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
	valuesBR[getSMMapKey(currMass, "gg", "BR")] = currIn[1];
	valuesBR[getSMMapKey(currMass, "gg", "+ERR")] = currIn[2];
	valuesBR[getSMMapKey(currMass, "gg", "-ERR")] = currIn[3];
	valuesBR[getSMMapKey(currMass, "gammagamma", "BR")] = currIn[4];
	valuesBR[getSMMapKey(currMass, "gammagamma", "+ERR")] = currIn[5];
	valuesBR[getSMMapKey(currMass, "gammagamma", "-ERR")] = currIn[6];
	valuesBR[getSMMapKey(currMass, "Zgamma", "BR")] = currIn[7];
	valuesBR[getSMMapKey(currMass, "Zgamma", "+ERR")] = currIn[8];
	valuesBR[getSMMapKey(currMass, "Zgamma", "-ERR")] = currIn[9];
	valuesBR[getSMMapKey(currMass, "WW", "BR")] = currIn[10];
	valuesBR[getSMMapKey(currMass, "WW", "+ERR")] = currIn[11];
	valuesBR[getSMMapKey(currMass, "WW", "-ERR")] = currIn[12];
	valuesBR[getSMMapKey(currMass, "ZZ", "BR")] = currIn[13];
	valuesBR[getSMMapKey(currMass, "ZZ", "+ERR")] = currIn[14];
	valuesBR[getSMMapKey(currMass, "ZZ", "-ERR")] = currIn[15];
      }
      
      else if (decayClass.Contains("2fermions")) {
	valuesBR[getSMMapKey(currMass, "bb", "BR")] = currIn[1];
	valuesBR[getSMMapKey(currMass, "bb", "+ERR")] = currIn[2];
	valuesBR[getSMMapKey(currMass, "bb", "-ERR")] = currIn[3];
	valuesBR[getSMMapKey(currMass, "tautau", "BR")] = currIn[4];
	valuesBR[getSMMapKey(currMass, "tautau", "+ERR")] = currIn[5];
	valuesBR[getSMMapKey(currMass, "tautau", "-ERR")] = currIn[6];
	valuesBR[getSMMapKey(currMass, "mumu", "BR")] = currIn[7];
	valuesBR[getSMMapKey(currMass, "mumu", "+ERR")] = currIn[8];
	valuesBR[getSMMapKey(currMass, "mumu", "-ERR")] = currIn[9];
	valuesBR[getSMMapKey(currMass, "cc", "BR")] = currIn[10];
	valuesBR[getSMMapKey(currMass, "cc", "+ERR")] = currIn[11];
	valuesBR[getSMMapKey(currMass, "cc", "-ERR")] = currIn[12];
	valuesBR[getSMMapKey(currMass, "ss", "BR")] = currIn[13];
	valuesBR[getSMMapKey(currMass, "ss", "+ERR")] = currIn[14];
	valuesBR[getSMMapKey(currMass, "ss", "-ERR")] = currIn[15];
	valuesBR[getSMMapKey(currMass, "tt", "BR")] = currIn[16];
	valuesBR[getSMMapKey(currMass, "tt", "+ERR")] = currIn[17];
	valuesBR[getSMMapKey(currMass, "tt", "-ERR")] = currIn[18];
      }
      
      else {
	std::cout << "BRXSReader: No decay class provided!" << std::endl;
      }
    }
  }
  currFile.close();
}

/**
   Load the Dark Matter model cross-sections from input files. 
   @returns - void.
*/
void BRXSReader::loadDMXS() {

  ifstream currFile(Form("%s/XS_DM.txt", directory.Data()));
  if (currFile.is_open()) {
    
    TString intermediateName;
    int intermediateMass;
    int fermionMass;
    float currIn[3];
    
    while (!currFile.eof()) {
      currFile >> intermediateName >> intermediateMass >> fermionMass
	       >> currIn[0] >> currIn[1] >> currIn[5];
      
      valuesXS[getDMMapKey(intermediateMass, fermionMass,
			   intermediateName, "XS")] = currIn[0];
      valuesXS[getDMMapKey(intermediateMass, fermionMass,
			   intermediateName, "+ERR")] = currIn[1];
      valuesXS[getDMMapKey(intermediateMass, fermionMass,
			   intermediateName, "-ERR")] = currIn[2];
    }
  }
  currFile.close();
}

/**
   Load the Standard Model cross-sections from input files.
   @param production - The production mode name. 
   @returns - void.
*/
void BRXSReader::loadSMXS(TString production) {
  
  TString fileName = Form("%s/XS_%s.txt", directory.Data(), production.Data());
  ifstream currFile(fileName);
  if (currFile.is_open()) {
    
    double currMass;
    float currIn[6];
    
    while (!currFile.eof()) {
      currFile >> currMass >> currIn[1] >> currIn[2] >> currIn[3] >> currIn[4]
	       >> currIn[5];
      
      valuesXS[getSMMapKey(currMass, production, "XS")] = currIn[1];
      valuesXS[getSMMapKey(currMass, production, "+QCD")] = currIn[2];
      valuesXS[getSMMapKey(currMass, production, "-QCD")] = currIn[3];
      valuesXS[getSMMapKey(currMass, production, "+PDF")] = currIn[4];
      valuesXS[getSMMapKey(currMass, production, "-PDF")] = currIn[5];
    }
  }
  currFile.close();
}

/**
   Convert variables to the key for the Dark Matter map.
   @param massIntermediate - The mass of the scalar or Zprime intermediary. 
   @param massFermion - The mass of the fermionic dark matter particle.
   @param type - The intermediate type (shxx_gg, zphxx_gg).
   @param value - The value (XS, +ERR, -ERR).
   @returns - the map key.
*/
TString BRXSReader::getDMMapKey(int massIntermediate, int massFermion,
				TString type, TString value) {
  TString key = Form("%s_%d_%d_%s", type.Data(), massIntermediate,
		     massFermion, value.Data());
  return key;
}

/**
   Convert variables to the key for the Standard Model map.
   @param mass - The signal mass of interest.
   @param type - The production or decay mode name. 
   @param value - The value (XS, +QCD, -QCD, +PDF, -PDF, BR, +ERR, -ERR).
   @returns - the map key.
*/
TString BRXSReader::getSMMapKey(double mass, TString type, TString value) {
  int massInt = (int)(100*mass);
  TString key = Form("mH%dGeV_%s_%s", massInt, type.Data(), value.Data());
  return key;
}

/**
   Check whether the specified key is contained in the specified map.
   @param key - The key for either the valuesXS or valuesBR map.
   @param mapType - The map type (either "XS", "BR").
   @returns - true if the specified key is in the map. 
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
