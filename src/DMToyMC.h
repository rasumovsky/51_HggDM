////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMToyMC.cxx                                                         //
//  Class: DMToyMC.cxx                                                        //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 17/04/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef TOY_H
#define TOY_H
#include "CommonHead.h"
#include "CommonFunc.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.hh"
#include "RooBernsteinM.h"

#include "DMAnalysis.h"

// INPUTS:
TString jobName;
TString option;
int seed;
int nToysPerJob;
int inputMuDM;
int input_lambda;
int input_lifetime;

TString outputDir;

vector<double> nToyEvents;
vector<double> valueNPMu1, valueNPMu0, valueNPMufree;
vector<double> value_globs_mu1, value_globs_mu0, value_globs_mufree;
vector<string> nameNP;
vector<string> nameGlobs;

bool firstPlot = true;
bool isFirstToy = true;
TCanvas *can;
