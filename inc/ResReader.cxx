////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: ResReader.cxx                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 15/04/2015                                                          //
//                                                                            //
//  This code reads the resolution systematic uncertainty table and provides  //
//  the values for the workspace code.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ResReader.h"

/**
   Initialize the ResReader class.
   @param inputFileName - the file location containing Res sys values.
   @param nCategories - the number of analysis categories.
   @returns - void.
*/
ResReader::ResReader(TString inputFileName, int nCategories) {
  
  std::cout << "\nInitializing the ResReader class" << std::endl;
  
  nResParams = 0;
  nameListRes.clear();
  
  TString tempSourceName = "";
  double tempValues[20] = {0.0};
  
  ifstream sysFile(inputFileName);
  while (!sysFile.eof()) {
    
    sysFile >> tempSourceName;
    for (int i_c = 0; i_c <= nCategories; i_c++) {
      sysFile >> tempValues[i_c];
      valuesRes[nResParams][i_c] = tempValues[i_c];
    }
    nameListRes.push_back(tempSourceName);
    nResParams++;
  }
  sysFile.close();
  std::cout << "Finished loading data for the ResReader class." << std::endl;
}

/**
   Finds the index associated with a particular Res sys source.
   @param name - the name of the Res sys.
   @returns - the index of the Res sys.
*/
int ResReader::sourceNameToIndex(TString name) {
  for (int i_n = 0; i_n < (int)nameListRes.size(); i_n++) {
    if (name.Contains(nameListRes[i_n])) {
      return i_n;
    }
  }
  return -1;
}

/**
   Get the value of a particular Res sys. in a given category.
   @param name - the name of the Res sys.
   @param cateIndex - the index of the category.
   @returns - the value of the Res sys. in the category.
*/
double ResReader::getValue(TString name, int cateIndex) {
  return fabs(valuesRes[sourceNameToIndex(name)][cateIndex]);
}

/**
   Returns the sign of the Res sys.
   @param name - the name of the Res sys.
   @param cateIndex - the index of the category.
   @returns - the sign of the Res sys. in the category.
*/
int ResReader::getSign(TString name, int cateIndex) {
  double value = valuesRes[sourceNameToIndex(name)][cateIndex];
  // make it absolute value (sign comes from GetSign():
  double sign = value / fabs(value);
  int result = ( sign >= 0 ) ? 1 : -1;
  return result;
}

/**
   Returns the number of Res sys. sources that have been defined.
   @returns - the number of Res sys. sources.
*/
int ResReader::getNumberOfSources() {
  return nResParams;
}

/**
   Retrieve the name of a source from its index.
   @param sourceIndex - the index of the Res sys.
   @returns - the name of the Res sys.
*/
TString ResReader::getNameOfSource(int resIndex) {
  return nameListRes[resIndex];
}
