////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: ESSReader.cxx                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
//  This code reads the ESS table and provides the values for the workspace.  //
//                                                                            //
//  WARNING! The category indexing is a bit unclear.                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ESSReader.h"

/**
   Initialize the ESSReader class.
   @param inputFileName - the file location containing ESS values.
   @param nCategories - the number of analysis categories.
   @returns - void.
*/
ESSReader::ESSReader(TString inputFileName, int nCategories) {
  
  std::cout << "\nInitializing the ESSReader class" << std::endl;
  
  nAnalysisCategories = nCategories;
  nESSParams = 0;
  
  nameListESS.clear();
  
  TString tempSourceName = "";
  double tempValues[20] = {0.0};
  
  ifstream sysFile(inputFileName);
  while( !sysFile.eof() )
  {
    sysFile >> tempSourceName;
    for (int i_c = 0; i_c <= nCategories; i_c++) {
      sysFile >> tempValues[i_c];
      valuesESS[nESSParams][i_c] = tempValues[i_c];
    }
    nameListESS.push_back(tempSourceName);
    nESSParams++;
  }
  sysFile.close();
  std::cout << "Finished loading data for the ESSReader class." << std::endl;
}

/**
   Finds the index associated with a particular ESS source.
   @param name - the name of the ESS.
   @returns - the index of the ESS.
*/
int ESSReader::sourceNameToIndex(TString name) {
  for (int i_n = 0; i_n < (int)nameListESS.size(); i_n++) {
    TString current_name = nameListESS[i_n];
    if (name.Contains(current_name)) {
      return i_n;
    }
  }
  return -1;
}

/**
   Get the value of a particular ESS in a given category.
   @param name - the name of the ESS.
   @param cateIndex - the index of the category.
   @returns - the value of the ESS in the category.
*/
double ESSReader::getValue(TString name, int cateIndex) {
  return fabs(valuesESS[sourceNameToIndex(name)][cateIndex]);
}

/**
   Returns the sign of the ESS.
   @param name - the name of the ESS.
   @param cateIndex - the index of the category.
   @returns - the sign of the ESS in the category.
*/
int ESSReader::getSign(TString name, int cateIndex) {
  double value = valuesESS[sourceNameToIndex(name)][cateIndex];
  // make it absolute value (sign comes from GetSign():
  double sign = value / fabs(value);
  int result = (sign >= 0) ? 1 : -1;
  return result;
}

/**
   Returns the number of ESS sources that have been defined.
   @returns - the number of ESS sources.
*/
int ESSReader::getNumberOfSources() {
  return nESSParams;
}

/**
   Retrieve the name of a source from its index.
   @param sourceIndex - the index of the ESS.
   @returns - the name of the ESS.
*/
TString ESSReader::getNameOfSource(int indexESS) {
  return nameListESS[indexESS];
}
