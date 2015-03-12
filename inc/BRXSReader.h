////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: BRXSReader.h                                                        //
//  Class: BRXSReader.cxx                                                     //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BRXSReader_h
#define BRXSReader_h

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <map>

class BRXSReader 
{
  
 public:
  
  BRXSReader(std::string inputDirectory);
  ~BRXSReader();
  
  // Accessors:
  double getXSection(int mass, std::string production, std::string value);
  double getBranchingRatio(int mass, std::string decay, std::string value);
    
 private:
  
  void loadProductionXS(std::string production);
  
  // Member objects:
  std::string directory;
  
  std::map<std::string,double> xSectionValues[3];
  
  
  std::map<int, std::map<std::string,double> > branchingRatios;

};

#endif
