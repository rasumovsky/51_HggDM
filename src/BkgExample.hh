//////////////////////////////////////////
//                                      //
//  SimultaneousBkg.hh                  //
//                                      //
//  Created: Andrew Hard 07/04/2013     //
//                                      //
//  Used by:                            //
//    SimultaneousBkg.cc                //
//                                      //
//////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "RooBernsteinM.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace CommonFunc;

TString output_directory;

// fit ranges:
double signal_band1 = 122.0;
double signal_band2 = 130.0;

int const NBINS = 10;

double luminosity = 20.7;

// signal parameters:
double param_values[NBINS][102][9];

//map<TString,double> map_fitresultSR;
vector<TString> vector_fitnames;
vector<vector<double> > vector_background;
vector<vector<double> > vector_coeff_Background;
vector<vector<double> > vector_coeff_Signal;
vector<vector<double> > vector_nbkg;

int FIT_PrintLevel = 0;
bool FIT_Verbose = true;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// LoadSignalParameterization:

void LoadSignalParameterization( TString param_directory )
{
  cout << "Loading Camila's signal parameterization text files." << endl;
  // loop over the bins:
  for( int i_b = 0; i_b < NBINS; i_b++ )
  {
    ifstream file_param(Form("%s/all_CosThetaStar%i_mc12a_commonPDF.txt",param_directory.Data(),i_b+1));
    int mass_index = 0;
    while( !file_param.eof() )
    {
      
      file_param >> param_values[i_b][mass_index][0] >> param_values[i_b][mass_index][1] >> param_values[i_b][mass_index][2] >> param_values[i_b][mass_index][3] >> param_values[i_b][mass_index][4] >> param_values[i_b][mass_index][5] >> param_values[i_b][mass_index][6] >> param_values[i_b][mass_index][7] >> param_values[i_b][mass_index][8];
      
      //cout << param_values[i_b][mass_index][0] << " " << param_values[i_b][mass_index][1] << " " << param_values[i_b][mass_index][2] << " " << param_values[i_b][mass_index][3] << " " << param_values[i_b][mass_index][4] << " " << param_values[i_b][mass_index][5] << " " << param_values[i_b][mass_index][6] << " " << param_values[i_b][mass_index][7] << " " << param_values[i_b][mass_index][8] << endl;
      
      mass_index++;
    }
    file_param.close();
  }
  cout << "Finished loading Camila's signal parameterization." << endl;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// MassToIndex:

int MassToIndex( double mass )
{
  int result = (int)(2 * ( mass - 100 ) );
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// PrintProgressBar:

void PrintProgressBar( int index, int total )
{
  if( index%10000 == 0 )
  {
    TString print_bar = " [";
    for( int bar = 0; bar < 20; bar++ )
    {
      double current_fraction = double(bar) / 20.0;
      if( double(index)/double(total) > current_fraction )
	print_bar.Append("/");
      else
	print_bar.Append(".");
    }
    print_bar.Append("] ");
    std::cout << print_bar << 100.*(double(index)/double(total)) << "%\r" << std::flush; 
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// FitOrder:

int FitOrder( TString fitname )
{
  int order = 0;
  for( int i_n = 1; i_n <= 8; i_n++ )
  {
    TString test = Form("ORDER%i",i_n);
    if( fitname.Contains(test) )
    {
      order = i_n;
      break;
    }
  }
  return order;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// NewFunctionName:

TString NewFunctionName( TString old_fitname, int param_change )
{
  int old_order = FitOrder( old_fitname );
  TString new_fitname = old_fitname;
  new_fitname.ReplaceAll(Form("ORDER%i",old_order),Form("ORDER%i",old_order+param_change));
  return new_fitname;
}

