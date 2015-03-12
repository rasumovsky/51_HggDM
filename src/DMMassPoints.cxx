////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMassPoints.cxx                                                    //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
//  This class uses a ROOT file with a TTree to implement the analysis event  //
//  categorization and produce text files or RooDataSets.                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMMassPoints.h"

/**
   Initialize the MassPoint class.
*/
DMMassPoints::DMMassPoints(TString newJobname, TString newOption) {
  std::cout << std::endl << "DMMassPoints::Initializing..." << std::endl;
  
  jobname = newJobname;
  option = newOption;
    
  return;
}

/**
 */


/**
   Return the number of events in a given category, or inclusively. Use 
   cateType="inclusive" to get the inclusive number of events.
*/
int getNEvents(TString cateType, int cate) {
  
}

/**
   Print the cutflow.
*/
void DMMassPoints::printCutflow() {
  
  std::cout << "Printing the cutflow:" << std::endl;
  for (int i_c = 0; i_c < 20; i_c++) {
    
    if (nPassing[i_c] == 0) {
      continue;
    }
    std::cout << "\t" << cutNames[i_c] << " \t" << nPassing[i_c] << "\t" << nPassing_weighted[20] << std::endl;
  }
  
  std::cout << "Printing the category results:" << std::endl;
  for (int i_b = 0; i_b < 20; i_b++) {
    
    if (
  
}
