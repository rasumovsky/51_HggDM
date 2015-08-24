//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 24 13:38:52 2015 by ROOT version 5.34/09
// from TTree CollectionTree/xAOD event tree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef DMTree_h
#define DMTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

using namespace std;
using std::vector;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxMET_MyRefFinalAux = 1;
const Int_t kMaxHGamEventInfoAux = 1;

class DMTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //xAOD::MissingETAuxContainer_v1 *MET_MyRefFinalAux_;
 //xAOD::MissingETAuxContainer_v1 *MET_MyRefFinalAux_xAOD__AuxContainerBase;
   vector<double>  MET_MyRefFinalAux_mpx;
   vector<double>  MET_MyRefFinalAux_mpy;
   vector<double>  MET_MyRefFinalAux_sumet;
   vector<string>  MET_MyRefFinalAux_name;
 //vector<ULong64_t> MET_MyRefFinalAux_source;
   //xAOD::MissingETContainer_v1 *MET_MyRefFinal;
   //xAOD::EventInfo_v1 *HGamEventInfo;
   //xAOD::AuxInfoBase *HGamEventInfoAux_;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_m;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_pt;
   Float_t         HGamEventInfoAuxDyn_weightInitial;
   Float_t         HGamEventInfoAuxDyn_HighMet_y1_pt;
   Float_t         HGamEventInfoAuxDyn_HighMet_y2_pt;
   Float_t         HGamEventInfoAuxDyn_HighMet_y1_eta;
   Float_t         HGamEventInfoAuxDyn_HighMet_y2_eta;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_eta;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_deltaR;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_deltaPhi;
   Float_t         HGamEventInfoAuxDyn_HighMet_yy_deltaEta;
   Double_t        HGamEventInfoAuxDyn_HighMet_MET_reb_CST1;
   Double_t        HGamEventInfoAuxDyn_HighMet_MET_reb_TST1;// use this ETMiss
   Float_t         HGamEventInfoAuxDyn_HighMet_dPhi_yy_MET;
   Int_t           HGamEventInfoAuxDyn_HighMet_lep_n2;
   Int_t           HGamEventInfoAuxDyn_HighMet_jet_n;

   // List of branches
   TBranch        *b_MET_MyRefFinalAux_mpx;   //!
   TBranch        *b_MET_MyRefFinalAux_mpy;   //!
   TBranch        *b_MET_MyRefFinalAux_sumet;   //!
   TBranch        *b_MET_MyRefFinalAux_name;   //!
   TBranch        *b_MET_MyRefFinal;   //!
   TBranch        *b_HGamEventInfo;   //!
   TBranch        *b_HGamEventInfoAux_;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_m;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_pt;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weightInitial;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_y1_pt;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_y2_pt;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_y1_eta;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_y2_eta;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_eta;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_deltaR;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_deltaPhi;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_yy_deltaEta;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_MET_reb_CST1;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_MET_reb_TST1;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_dPhi_yy_MET;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_lep_n2;   //!
   TBranch        *b_HGamEventInfoAuxDyn_HighMet_jet_n;   //!

   DMTree(TTree *tree=0);
   virtual ~DMTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DMTree_cxx
DMTree::DMTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sample.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sample.root");
      }
      f->GetObject("CollectionTree",tree);

   }
   Init(tree);
}

DMTree::~DMTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DMTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DMTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DMTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //MET_MyRefFinal = 0;
   //HGamEventInfo = 0;
   //HGamEventInfoAux_ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("MET_MyRefFinalAux.mpx", &MET_MyRefFinalAux_mpx, &b_MET_MyRefFinalAux_mpx);
   fChain->SetBranchAddress("MET_MyRefFinalAux.mpy", &MET_MyRefFinalAux_mpy, &b_MET_MyRefFinalAux_mpy);
   fChain->SetBranchAddress("MET_MyRefFinalAux.sumet", &MET_MyRefFinalAux_sumet, &b_MET_MyRefFinalAux_sumet);
   fChain->SetBranchAddress("MET_MyRefFinalAux.name", &MET_MyRefFinalAux_name, &b_MET_MyRefFinalAux_name);
   //fChain->SetBranchAddress("MET_MyRefFinal", &MET_MyRefFinal, &b_MET_MyRefFinal);
   //fChain->SetBranchAddress("HGamEventInfo", &HGamEventInfo, &b_HGamEventInfo);
   //fChain->SetBranchAddress("HGamEventInfoAux.", &HGamEventInfoAux_, &b_HGamEventInfoAux_);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_m", &HGamEventInfoAuxDyn_HighMet_yy_m, &b_HGamEventInfoAuxDyn_HighMet_yy_m);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_pt", &HGamEventInfoAuxDyn_HighMet_yy_pt, &b_HGamEventInfoAuxDyn_HighMet_yy_pt);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightInitial", &HGamEventInfoAuxDyn_weightInitial, &b_HGamEventInfoAuxDyn_weightInitial);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_y1_pt", &HGamEventInfoAuxDyn_HighMet_y1_pt, &b_HGamEventInfoAuxDyn_HighMet_y1_pt);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_y2_pt", &HGamEventInfoAuxDyn_HighMet_y2_pt, &b_HGamEventInfoAuxDyn_HighMet_y2_pt);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_y1_eta", &HGamEventInfoAuxDyn_HighMet_y1_eta, &b_HGamEventInfoAuxDyn_HighMet_y1_eta);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_y2_eta", &HGamEventInfoAuxDyn_HighMet_y2_eta, &b_HGamEventInfoAuxDyn_HighMet_y2_eta);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_eta", &HGamEventInfoAuxDyn_HighMet_yy_eta, &b_HGamEventInfoAuxDyn_HighMet_yy_eta);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_deltaR", &HGamEventInfoAuxDyn_HighMet_yy_deltaR, &b_HGamEventInfoAuxDyn_HighMet_yy_deltaR);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_deltaPhi", &HGamEventInfoAuxDyn_HighMet_yy_deltaPhi, &b_HGamEventInfoAuxDyn_HighMet_yy_deltaPhi);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_yy_deltaEta", &HGamEventInfoAuxDyn_HighMet_yy_deltaEta, &b_HGamEventInfoAuxDyn_HighMet_yy_deltaEta);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_MET_reb_CST1", &HGamEventInfoAuxDyn_HighMet_MET_reb_CST1, &b_HGamEventInfoAuxDyn_HighMet_MET_reb_CST1);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_MET_reb_TST1", &HGamEventInfoAuxDyn_HighMet_MET_reb_TST1, &b_HGamEventInfoAuxDyn_HighMet_MET_reb_TST1);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_dPhi_yy_MET", &HGamEventInfoAuxDyn_HighMet_dPhi_yy_MET, &b_HGamEventInfoAuxDyn_HighMet_dPhi_yy_MET);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_lep_n2", &HGamEventInfoAuxDyn_HighMet_lep_n2, &b_HGamEventInfoAuxDyn_HighMet_lep_n2);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.HighMet_jet_n", &HGamEventInfoAuxDyn_HighMet_jet_n, &b_HGamEventInfoAuxDyn_HighMet_jet_n);
   Notify();
}

Bool_t DMTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DMTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DMTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DMTree_cxx
