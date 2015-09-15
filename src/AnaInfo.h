////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  AnaInfo.h                                                                 //
//  Class: AnaInfo.cxx                                                        //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 13/09/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef AnaInfo_h
#define AnaInfo_h

#include "CommonHead.h"
//#include "CommonFunc.h"

class AnaInfo {

 public:
  
  // Class constructor and destructor:
  AnaInfo(int index);
  virtual ~AnaInfo() {};
  
  // Public accessors:
  std::vector<TString> getCutList();
  double getAnaP0(TString signal, bool observed);
  double getAnaCL(TString signal, bool observed);
  double getCutVal(TString cutName);
  int getIndex();
  std::vector<TString> getSignalList();
  double getStatVal(TString signal, TString statistic);
  bool isGood();
  void printAna(TString signal);
  
  // Public mutators:
  void setCutVal(TString cutName, double cutVal);
  void setGood(bool isGood);
  void setIndex(int index);
  void setStatVal(TString signal, TString statistic, double value);
  
 private:

  // Private mutators:
  void addCut(TString cutName);
  void addSignal(TString signalName);
  
  // Private member variables:
  int m_index; // Job index.
  bool m_isGood; // Data loaded successfully?
  std::map<TString,double> m_cutValues; // Cut names and their values.
  std::vector<TString> m_signals; // A list of signal names
  std::vector<TString> m_cutNames; // A list of cut names
  std::map<TString,double> m_statValues; // statistics for ALL signals.
  
};

#endif

