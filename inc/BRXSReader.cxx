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
////////////////////////////////////////////////////////////////////////////////

#include "BRXSReader.h"

/**
   Initializes the tool and loads XS, BR values from files. 
  */
BRXSReader::BRXSReader(std::string inputDirectory) {
  loadProductionXS("ggH");
  loadProductionXS("VBF");
  loadProductionXS("WH");
  loadProductionXS("ZH");
  loadProductionXS("ttH");
  loadProductionXS("bbH");
  
}

/**
   Returns the cross-section, or related uncertainties, for a production 
   process at a particular mass.
*/
double BRXSReader::getXSection(int mass, std::string production, 
			       std::string value) {
  
}

/**
   Returns the branching-ratio, or related uncertainties, for a decay mode
   at a particular Higgs mass. 
 */
double BRXSReader::getBranchingRatio(int mass, std::string decay, 
				     std::string value) {
  
}

/**
   Load the cross sections from input files.
 */
void BRXSReader::loadProductionXS(std::string production) {
  char fileName[100];
  int n = sprintf(buffer, "%s/xsection_%s.txt", directory, inputDirectory);
  ifstream currFile(fileName);
  if (currFile.is_open()) {
    
    double inMass; double inXS; double inPQCD; double inMQCD; double inPPDF; double inMPDF;
    
    while (!currFile.eof()) {
      currFile >> inMass >> inXS inPQCD >> inMQCD >> inPPDF >> inMPDF;
      int saveMass = (int)(10*inMass);
      xSections[saveMass][production] = inXS;
      plusQCD[saveMass][production] = inPQCD/100.0;
      minusQCD[saveMass][production] = inMQCD/100.0;
      plusPDF[saveMass][production] = inPPDF/100.0;
      minusPDF[saveMass][production] = inMPDF/100.0;
      
    }
  }
  currFile.close();
}
