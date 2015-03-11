/////////////////////////////////////////
//                                     //
//  spinv1_masspoints.hh               //
//                                     //
//  Author: Andrew Hard                //
//  Date: 21/01/2014                   //
//                                     //
/////////////////////////////////////////

#include "../inc/CommonHead.h"
#include "../inc/CommonFunc.h"
#include "../inc/RooFitHead.h"
#include "../inc/statistics.hh"
#include "../inc/ntup_SpinPub.hh"

#include "../src/spinv1_Master.hh"

using namespace std;
using namespace CommonFunc;

int energy;// 7 or 8

// ntuple class:
ntup_SpinPub* p;

TString output_directory;
TString jobname;

TH1F *h_relpt1;
TH1F *h_relpt2;
TH2D *h2_relpt;
TH1F *h_sinhdeta;

TH1F *h_mass;
TH1F *h_cts;
TH1F *h_ptt;
TH1F *h_pt_gg;
TH1F *h_pt_g1;
TH1F *h_pt_g2;
TH1F *h_eta_g1;
TH1F *h_eta_g2;
TH1F *h_sinh_deta;
TH2D *h2_cts_ptt;
TH2D *h2_cts_pt_gg;
TH2D *h2_cts_sinh_deta;
TH2D *h2_sinh_ptgg;

TCanvas *can;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// get_event:

void get_event( int i ) 
{
  if ( p->LoadTree(i) < 0) 
  { 
    cout << "\nProblem in LoadTree." << "\nEntry: " <<i<<endl;
    exit(0);
  }
  p->fChain->GetEntry(i);
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
      if( double(index)/double(total) > current_fraction ) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.*(double(index)/double(total));
    TString text = Form("%s %2.2f ",print_bar.Data(),percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetCTSPTBin:

int GetCTSPTBin( int nbins, double cts, double pt )
{
  int result = -1;
  for( int b = 0; b < 10; b++ )
  {
    double bin_lo = ( (double) b ) / 10;
    double bin_hi = ( (double) b+1 ) / 10;
    if( cts > bin_lo && cts < bin_hi ){ result = b; break; }
  }
  if( pt > 125000 ) result = 10;
  if( result == -1 ){ cout << "Problem! No bin assignment." << endl; exit(0); }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// GetCTSBin:

int GetCTSBin( int nbins, double cts )
{
  int result = -1;
  for( int b = 0; b < nbins; b++ )
  {
    double bin_lo = ( (double) b ) / ( (double) nbins );
    double bin_hi = ( (double) b+1 ) / ( (double) nbins );
    if( cts > bin_lo && cts < bin_hi ){ result = b; break; }
  }
  if( result == -1 ){ cout << "Problem! No bin assignment." << endl; exit(0); }
  return result;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// DrawRelPtPlots:

void DrawRelPtPlots()
{
  can->cd();
  can->SetRightMargin(0.1);
  h_relpt1->GetXaxis()->SetTitle("p_{T}(#gamma1) / m_{#gamma#gamma}");
  h_relpt1->GetYaxis()->SetTitle("Entries");
  h_relpt1->SetLineColor(kBlue);
  h_relpt1->Draw("");
  can->Print(Form("%s/plots/plot_RelPt_g1.eps",output_directory.Data()));
  can->Clear();
  
  h_relpt2->GetXaxis()->SetTitle("p_{T}(#gamma2) / m_{#gamma#gamma}");
  h_relpt2->GetYaxis()->SetTitle("Entries");
  h_relpt2->SetLineColor(kBlue);
  h_relpt2->Draw("");
  can->Print(Form("%s/plots/plot_RelPt_g2.eps",output_directory.Data()));
  can->Clear();
  
  h_sinhdeta->GetXaxis()->SetTitle("sinh(#Delta#eta_{#gamma#gamma})");
  h_sinhdeta->GetYaxis()->SetTitle("Entries");
  h_sinhdeta->SetLineColor(kBlue);
  h_sinhdeta->Draw("");
  can->Print(Form("%s/plots/plot_sinhdeta.eps",output_directory.Data()));
  can->Clear();
  
  gStyle->SetPalette(1);
  gPad->SetLogz();
  h2_relpt->GetXaxis()->SetTitle("p_{T}(#gamma1) / m_{#gamma#gamma}");
  h2_relpt->GetYaxis()->SetTitle("p_{T}(#gamma2) / m_{#gamma#gamma}");
  h2_relpt->GetZaxis()->SetTitle("Entries");
  h2_relpt->Draw("COLZ");
  can->Print(Form("%s/plots/plot_RelPt_2D.eps",output_directory.Data()));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// InitializeKinematicsHists:

void InitializeKinematicHists()
{
  h_mass = new TH1F("h_mass","h_mass",55,105,160);
  h_cts = new TH1F("h_cts","h_cts",10,0,1);
  h_ptt = new TH1F("h_ptt","h_ptt",40,0,120);
  h_pt_gg = new TH1F("h_pt_gg","h_pt_gg",40,0,120);
  h_pt_g1 = new TH1F("h_pt_g1","h_pt_g1",40,35,135);
  h_pt_g2 = new TH1F("h_pt_g2","h_pt_g2",40,25,125);
  h_eta_g1 = new TH1F("h_eta_g1","h_eta_g1",40,-2.4,2.5);
  h_eta_g2 = new TH1F("h_eta_g2","h_eta_g2",40,-2.4,2.5);
  h_sinh_deta = new TH1F("h_sinh_deta","h_sinh_deta",40,-8,8);

  h_mass->Sumw2();
  h_cts->Sumw2();
  h_ptt->Sumw2();
  h_pt_gg->Sumw2();
  h_pt_g1->Sumw2();
  h_pt_g2->Sumw2();
  h_eta_g1->Sumw2();
  h_eta_g2->Sumw2();
  h_sinh_deta->Sumw2();
  
  h2_cts_ptt = new TH2D("h2_cts_ptt","h2_cts_ptt",10,0,1,80,0,120);
  h2_cts_pt_gg = new TH2D("h2_cts_pt_gg","h2_cts_pt_gg",10,0,1,80,0,120);
  h2_cts_sinh_deta = new TH2D("h2_cts_sinh_deta","h2_cts_sinh_deta",10,0,1,80,-8,8); 
  h2_sinh_ptgg = new TH2D("h2_sinh_ptgg","h2_sinh_ptgg",4,0,4,7,0,70);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// FillKinematicHists:

void FillKinematicHists( TLorentzVector g1, TLorentzVector g2, double cts, double ptt, double weight )
{
  TLorentzVector gg = g1 + g2;
  h_mass->Fill( gg.M(), weight );
  h_cts->Fill( cts, weight );
  h_ptt->Fill( ptt, weight );
  h_pt_gg->Fill( gg.Pt(), weight );
  h_pt_g1->Fill( g1.Pt(), weight );
  h_pt_g2->Fill( g2.Pt(), weight );
  h_eta_g1->Fill( g1.Eta(), weight );
  h_eta_g2->Fill( g2.Eta(), weight );
  h_sinh_deta->Fill( TMath::SinH( g1.Eta() - g2.Eta() ), weight );
  h2_cts_ptt->Fill( cts, ptt, weight );
  h2_cts_pt_gg->Fill( cts, gg.Pt(), weight );
  h2_cts_sinh_deta->Fill( cts, TMath::SinH( g1.Eta() - g2.Eta() ), weight );
  double test_sinh = TMath::SinH( g1.Eta() - g2.Eta() );
  if( test_sinh >= 4.0 ) test_sinh = 3.999;
  double test_ptgg = gg.Pt();
  if( test_ptgg >= 70.0 ) test_ptgg = 69.999;
  h2_sinh_ptgg->Fill( fabs(test_sinh), test_ptgg, weight );
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// SaveKinematicHists:

void SaveKinematicHists()
{
  TFile *output_file = new TFile(Form("%s/hists.root",output_directory.Data()),"RECREATE");
  
  h_mass->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
  h_mass->GetYaxis()->SetTitle("Entries / GeV");
  h_mass->Write();
  
  h_cts->GetXaxis()->SetTitle("cos(#theta*)_{CS}");
  h_cts->GetYaxis()->SetTitle("Entries / 0.1");
  h_cts->Write();
  
  h_ptt->GetXaxis()->SetTitle("p_{Tt} [GeV]");
  h_ptt->GetYaxis()->SetTitle("Entries");
  h_ptt->Write();
  
  h_pt_gg->GetXaxis()->SetTitle("p_{T}(#gamma#gamma) [GeV]");
  h_pt_gg->GetYaxis()->SetTitle("Entries");
  h_pt_gg->Write();
  
  h_pt_g1->GetXaxis()->SetTitle("p_{T}(#gamma1) [GeV]");
  h_pt_g1->GetYaxis()->SetTitle("Entries");
  h_pt_g1->Write();
  
  h_pt_g2->GetXaxis()->SetTitle("p_{T}(#gamma2) [GeV]");
  h_pt_g2->GetYaxis()->SetTitle("Entries");
  h_pt_g2->Write();
  
  h_eta_g1->GetXaxis()->SetTitle("#eta_{#gamma1}");
  h_eta_g1->GetYaxis()->SetTitle("Entries");
  h_eta_g1->Write();
  
  h_eta_g2->GetXaxis()->SetTitle("#eta_{#gamma2}");
  h_eta_g2->GetYaxis()->SetTitle("Entries");
  h_eta_g2->Write();
  
  h_sinh_deta->GetXaxis()->SetTitle("sinh(#Delta#eta_{#gamma#gamma})");
  h_sinh_deta->GetYaxis()->SetTitle("Entries");
  h_sinh_deta->Write();  
  
  h2_cts_ptt->GetXaxis()->SetTitle("cos(#theta*)");
  h2_cts_ptt->GetYaxis()->SetTitle("p_{Tt} [GeV]");
  h2_cts_ptt->Write();
  
  h2_cts_pt_gg->GetXaxis()->SetTitle("cos(#theta*)");
  h2_cts_pt_gg->GetYaxis()->SetTitle("p_{T}(#gamma#gamma) [GeV]");
  h2_cts_pt_gg->Write();
  
  h2_cts_sinh_deta->GetXaxis()->SetTitle("cos(#theta*)");
  h2_cts_sinh_deta->GetYaxis()->SetTitle("sinh(#Delta#eta_{#gamma#gamma})");
  h2_cts_sinh_deta->Write();
  
  h2_sinh_ptgg->GetXaxis()->SetTitle("sinh(#Delta#eta)");
  h2_sinh_ptgg->GetYaxis()->SetTitle("p_{T}(#gamma#gamma) [GeV]");
  h2_sinh_ptgg->Write();
  
  output_file->Write();
  output_file->Close();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////// CompareTwoFiles:

void CompareTwoFiles( TString file1, TString file2 )
{
  // 2nd file should have the extra events.
  
  ofstream differences;
  differences.open("Difference.txt");
  int ndifferences = 0;
  
  int Run2, Event2;
  ifstream list2(file2);
  cout << "Opened file 2." << endl;
  int counter = 0;
  while( !list2.eof() )
  {
    counter++;
    bool match = false;
    list2 >> Run2 >> Event2;
    if( counter%1000 == 0 ) cout << counter << endl;
    
    int Run1, Event1;
    bool zveto1;
    ifstream list1(file1);
    while( !list1.eof() )
    {
      list1 >> Run1 >> Event1;
      if( Event1 == Event2 && Run1 == Run2 )
      {
	match = true;
	break;
      }
    }
    list1.close();
    if( !match ) 
    {
      cout << "No match for run " << Run2 << " event: " << Event2 << endl;
      differences << Event2 << " " << endl;
      ndifferences++;
    }
  }
  cout << "NUMBER OF DIFFERENCES = " << ndifferences << endl;
}
