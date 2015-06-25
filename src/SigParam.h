////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SigParam.h                                                          //
//  Class: SigParam.cxx                                                       //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/06/2015                                                          //
//                                                                            //
//  Accessors access class data without modifying the member objects, while   //
//  mutators modify the state of the class (and also sometimes return data.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SigParam_h
#define SigParam_h

// C++ includes:
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

// ROOT includes:
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

// RooFit headers:
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsRealLValue.h>
#include <RooAddition.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooBernstein.h>
#include <RooCategory.h>
#include <RooCBShape.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooFormula.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooNumIntConfig.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooProfileLL.h>
#include <RooRealSumPdf.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

class SigParam {
  
 public:
  
  // Constructor and destructor:
  SigParam(TString options);
  virtual ~SigParam() {};
  
  // Public Accessors:
  void addSigToWS(RooWorkspace *&workspace, RooArgSet *&nuisParams, 
		  int cateIndex);
  void addSigToWS(RooWorkspace *&workspace, RooArgSet *&nuisParams,
		  double resonanceMass, int cateIndex);
  double getParameterError(TString paramName, double resonanceMass,
			   int cateIndex);
  double getParameterError(TString paramName, int cateIndex);
  double getParameterValue(TString paramName, double resonanceMass, 
			   int cateIndex);
  double getParameterValue(TString paramName, int cateIndex);
  RooAbsPdf* getSingleResonance(double resonanceMass, int cateIndex);
  RooWorkspace* getWorkspace();
  double getYieldInCategory(double resonanceMass, int cateIndex);
  double getYieldTotal(double resonanceMass);
  
  // Public Mutators:
  void addMResSystematic(TString nameMResSys);
  void addMResSystematics(std::vector<TString> nameMResSys);
  void addMScaleSystematic(TString nameMScaleSys);
  void addMScaleSystematics(std::vector<TString> nameMScaleSys);
  //void addCombData(RooDataSet* dataSet, TString observableName,/// !
  //		   TString weightName, TString cateName, TString resMassName);
  void addDataSet(double resonanceMass, int cateIndex, RooDataSet* dataSet,
		  TString observableName);
  void addDataTree(double resonanceMass, int cateIndex, TTree *dataTree,
		   TString massBranchName, TString weightBranchName);
  void addMassPoint(double resonanceMass, int cateIndex, double diphotonMass,
		    double eventWeight);
  bool loadParameterization(TString fileName);
  bool makeAllParameterizations(TString function);
  bool makeCategoryParameterization(int cateIndex, TString function);
  bool makeSingleResonance(double resonanceMass, int cateIndex,
			   TString function);
  //void mergeParameterizations(SigParam *otherParam);/// !
  void plotCategoryResonances(int cateIndex, TString fileName);
  void plotSingleResonance(double resonanceMass, int cateIndex,
			   TString fileName);
  void saveParameterization(TString fileName);
  void setVarParameterization(TString varName, TString function);
  
 private:
  
  // Private Accessors:
  std::vector<int> categoriesForMass(double resonanceMass);
  bool dataExists(double resonanceMass, int cateIndex);
  bool equalMasses(double massValue1, double massValue2);
  TString getKey(double resonanceMass, int cateIndex);
  double massIntToDouble(int massInteger);
  int massDoubleToInt(double resonanceMass);
  std::vector<double> massPointsForCategory(int cateIndex);
  
  // Private Mutators:
  //RooFormulaVar* constructFormulaVar(TString varName, int polyOrder);
  RooFitResult* fitResult(int cateIndex);
  RooFitResult* fitResult(double resonanceMass, int cateIndex);
  int getNCategories();
  double regularizedMass(double resonanceMass);
  void resonanceCreator(double resonanceMass, int cateIndex, TString function);
  void setParamsConstant(RooAbsPdf* pdf, bool isConstant);
  
  // Member variables:
  std::vector<std::pair<double,int> > m_massCatePairs;
  std::map<TString,TString> m_funcList;
  
  int m_nCategories;
  TString m_options;
  
  // Objects for fitting:
  RooRealVar *m_yy;
  RooRealVar *m_wt;
  RooWorkspace *m_ws;
  RooCategory *m_cat;
  
  // Shape information:
  TString m_listMRS;
  TString m_listMSS;
  
};

#endif
