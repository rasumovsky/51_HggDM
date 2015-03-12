////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: DMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 10/03/2015                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DMHeader.h"

using namespace std;

void MakeExe(TString exename);

void SubmitWSViaBsub(TString executable_name, TString executable_jobname, TString executable_option, int executable_lambda, int executable_lifetime);

void SubmitToysViaBsub(TString executable_name, TString executable_jobname, TString executable_option, int executable_seed, int executable_toys_per_job, int chosen_lambda, int chosen_lifetime);
