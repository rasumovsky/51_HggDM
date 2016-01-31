//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 15:01:29 2015 by ROOT version 6.02/12
// from TTree CollectionTree/xAOD event tree
// found on file: PowhegPy8_VBF125.MxAOD.p2421.h008.root
//////////////////////////////////////////////////////////

#ifndef DMTree_h
#define DMTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>

using namespace std;
using std::vector;

class DMTree {
 
 public :
  
  TTree *fChain;  //!pointer to the analyzed TTree or TChain
  Int_t fCurrent; //!current Tree number in a TChain
  
  bool m_useSys;
  std::vector<TString> m_sysNames;
  
  // Declaration of nominal and systematics leaf types:
  std::map<TString, Float_t> HGamEventInfoAuxDyn_m_yy;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_pT_yy;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_cosTS_yy;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_Njets;
  std::map<TString, Int_t>   HGamEventInfoAuxDyn_cutFlow;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_weightInitial;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_weight;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_TST_met;
  std::map<TString, Float_t> HGamEventInfoAuxDyn_crossSectionBRfilterEff;
  
  std::vector<float> *HGamElectronsAuxDyn_pt;
  std::vector<float> *HGamMuonsAuxDyn_pt;
  std::vector<float> *HGamTruthHiggsBosonsAuxDyn_m;
  
  std::vector<float> *HGamPhotonsAuxDyn_pt;
  std::vector<float> *HGamPhotonsAuxDyn_eta;
  std::vector<float> *HGamPhotonsAuxDyn_phi;
  std::vector<float> *HGamPhotonsAuxDyn_m;
  
  std::vector<float> *HGamAntiKt4EMTopoJetsAuxDyn_pt;
  std::vector<float> *HGamAntiKt4EMTopoJetsAuxDyn_eta;
  std::vector<float> *HGamAntiKt4EMTopoJetsAuxDyn_phi;
  std::vector<float> *HGamAntiKt4EMTopoJetsAuxDyn_m;
  
  // List of branches
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_m_yy;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_pT_yy;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_cosTS_yy;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_Njets;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_cutFlow;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_weightInitial;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_weight;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_TST_met;
  std::map<TString, TBranch*> b_HGamEventInfoAuxDyn_crossSectionBRfilterEff;
  
  TBranch *b_HGamElectronsAuxDyn_pt;
  TBranch *b_HGamMuonsAuxDyn_pt;
  TBranch *b_HGamTruthHiggsBosonsAuxDyn_m;
  
  TBranch *b_HGamPhotonsAuxDyn_pt;
  TBranch *b_HGamPhotonsAuxDyn_eta;
  TBranch *b_HGamPhotonsAuxDyn_phi;
  TBranch *b_HGamPhotonsAuxDyn_m;
  
  TBranch *b_HGamAntiKt4EMTopoJetsAuxDyn_pt;
  TBranch *b_HGamAntiKt4EMTopoJetsAuxDyn_eta;
  TBranch *b_HGamAntiKt4EMTopoJetsAuxDyn_phi;
  TBranch *b_HGamAntiKt4EMTopoJetsAuxDyn_m;
  
  // Methods:
  DMTree(TTree *tree, std::vector<TString> sysNames);
  virtual ~DMTree();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
  
  void DefineSystematics(std::vector<TString> sysNames);

};

#endif

#ifdef DMTree_cxx
DMTree::DMTree(TTree *tree, std::vector<TString> sysNames) : fChain(0) {
  
  m_sysNames.clear();
  m_sysNames.push_back("Nominal");
  for (int i_s = 0; i_s < (int) sysNames.size(); i_s++) {
    m_sysNames.push_back(sysNames[i_s]);
  }
  Init(tree);
}

DMTree::~DMTree() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DMTree::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t DMTree::LoadTree(Long64_t entry) {
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

void DMTree::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
   // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
    
  // Set object pointer
  HGamElectronsAuxDyn_pt = 0;
  HGamMuonsAuxDyn_pt = 0;
  HGamTruthHiggsBosonsAuxDyn_m = 0;
  
  // Loop over systematics:
  for (int i_s = 0; i_s < (int)m_sysNames.size(); i_s++) {
    
    // Create a map tag (key for maps) and data tag (name for MxAOD branches).
    TString mTag = m_sysNames[i_s];
    TString dTag = (mTag.EqualTo("Nominal")) ?
      "" : Form("_%s", m_sysNames[i_s].Data());
  
    // Someone carelessly put a space in the branch name:
    if (dTag.Contains("FT_EFF_extrapolation_from_charm")) {
      dTag.ReplaceAll("FT_EFF_extrapolation_from_charm",
		      "FT_EFF_extrapolation from charm");
    }

    b_HGamEventInfoAuxDyn_m_yy[mTag] = NULL;
    b_HGamEventInfoAuxDyn_pT_yy[mTag] = NULL;
    b_HGamEventInfoAuxDyn_cosTS_yy[mTag] = NULL;
    b_HGamEventInfoAuxDyn_Njets[mTag] = NULL;
    b_HGamEventInfoAuxDyn_cutFlow[mTag] = NULL;
    b_HGamEventInfoAuxDyn_weightInitial[mTag] = NULL;
    b_HGamEventInfoAuxDyn_weight[mTag] = NULL;
    b_HGamEventInfoAuxDyn_TST_met[mTag] = NULL;
    b_HGamEventInfoAuxDyn_crossSectionBRfilterEff[mTag] = NULL;
        
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.m_yy",dTag.Data()),
			     &HGamEventInfoAuxDyn_m_yy[mTag],
			     &b_HGamEventInfoAuxDyn_m_yy[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.pT_yy",dTag.Data()),
			     &HGamEventInfoAuxDyn_pT_yy[mTag],
			     &b_HGamEventInfoAuxDyn_pT_yy[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.cosTS_yy",dTag.Data()),
			     &HGamEventInfoAuxDyn_cosTS_yy[mTag],
			     &b_HGamEventInfoAuxDyn_cosTS_yy[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.Njets",dTag.Data()),
			     &HGamEventInfoAuxDyn_Njets[mTag],
			     &b_HGamEventInfoAuxDyn_Njets[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.cutFlow",dTag.Data()), 
			     &HGamEventInfoAuxDyn_cutFlow[mTag],
			     &b_HGamEventInfoAuxDyn_cutFlow[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.weightInitial",
				  dTag.Data()),
			     &HGamEventInfoAuxDyn_weightInitial[mTag], 
			     &b_HGamEventInfoAuxDyn_weightInitial[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.weight",dTag.Data()),
			     &HGamEventInfoAuxDyn_weight[mTag],
			     &b_HGamEventInfoAuxDyn_weight[mTag]);
    fChain->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.TST_met",dTag.Data()),
			     &HGamEventInfoAuxDyn_TST_met[mTag],
			     &b_HGamEventInfoAuxDyn_TST_met[mTag]);
    fChain
      ->SetBranchAddress(Form("HGamEventInfo%sAuxDyn.crossSectionBRfilterEff",
			      dTag.Data()),
			 &HGamEventInfoAuxDyn_crossSectionBRfilterEff[mTag],
			 &b_HGamEventInfoAuxDyn_crossSectionBRfilterEff[mTag]);
  }
  
  // The branches that don't have systematic variations:
  fChain->SetBranchAddress("HGamElectronsAuxDyn.pt",
			   &HGamElectronsAuxDyn_pt,
			   &b_HGamElectronsAuxDyn_pt);
  fChain->SetBranchAddress("HGamMuonsAuxDyn.pt",
			   &HGamMuonsAuxDyn_pt,
			   &b_HGamMuonsAuxDyn_pt);
  fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.m",
			   &HGamTruthHiggsBosonsAuxDyn_m,
			   &b_HGamTruthHiggsBosonsAuxDyn_m);
  
  fChain->SetBranchAddress("HGamPhotonsAuxDyn_pt", &HGamPhotonsAuxDyn_pt, 
			   &b_HGamPhotonsAuxDyn_pt);
  fChain->SetBranchAddress("HGamPhotonsAuxDyn_eta", &HGamPhotonsAuxDyn_eta,
			   &b_HGamPhotonsAuxDyn_eta);
  fChain->SetBranchAddress("HGamPhotonsAuxDyn_phi", &HGamPhotonsAuxDyn_phi,
			   &b_HGamPhotonsAuxDyn_phi);
  fChain->SetBranchAddress("HGamPhotonsAuxDyn_m", &HGamPhotonsAuxDyn_m,
			   &b_HGamPhotonsAuxDyn_m);
  
  fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn_pt", 
			   &HGamAntiKt4EMTopoJetsAuxDyn_pt,
			   &b_HGamAntiKt4EMTopoJetsAuxDyn_pt);
  fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn_eta",
			   &HGamAntiKt4EMTopoJetsAuxDyn_eta,
			   &b_HGamAntiKt4EMTopoJetsAuxDyn_eta);
  fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn_phi",
			   &HGamAntiKt4EMTopoJetsAuxDyn_phi,
			   &b_HGamAntiKt4EMTopoJetsAuxDyn_phi);
  fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn_m",
			   &HGamAntiKt4EMTopoJetsAuxDyn_m,
			   &b_HGamAntiKt4EMTopoJetsAuxDyn_m);
  
  Notify();
}

Bool_t DMTree::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   return kTRUE;
}

void DMTree::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t DMTree::Cut(Long64_t entry) {
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef DMTree_cxx
