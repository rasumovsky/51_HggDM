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
  
  // Methods:
  DMTree(TTree *tree = 0);
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
DMTree::DMTree(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PowhegPy8_VBF125.MxAOD.p2421.h008.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("PowhegPy8_VBF125.MxAOD.p2421.h008.root");
    }
    f->GetObject("CollectionTree",tree);
    
  }
  Init(tree);
  m_sysNames.clear();
  m_sysNames.push_back("Nominal");
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
    
    TString tag = "";;
    if (!tag.EqualTo("Nominal")) {
      tag = Form("_%s", m_sysNames[i_s].Data());
    }
    
    b_HGamEventInfoAuxDyn_m_yy[tag] = NULL;
    b_HGamEventInfoAuxDyn_pT_yy[tag] = NULL;
    b_HGamEventInfoAuxDyn_cosTS_yy[tag] = NULL;
    b_HGamEventInfoAuxDyn_Njets[tag] = NULL;
    b_HGamEventInfoAuxDyn_cutFlow[tag] = NULL;
    b_HGamEventInfoAuxDyn_weightInitial[tag] = NULL;
    b_HGamEventInfoAuxDyn_weight[tag] = NULL;
    b_HGamEventInfoAuxDyn_TST_met[tag] = NULL;
    b_HGamEventInfoAuxDyn_crossSectionBRfilterEff[tag] = NULL;
        
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.m_yy",tag.Data()),
			     &HGamEventInfoAuxDyn_m_yy[tag],
			     &b_HGamEventInfoAuxDyn_m_yy[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.pT_yy",tag.Data()),
			     &HGamEventInfoAuxDyn_pT_yy[tag],
			     &b_HGamEventInfoAuxDyn_pT_yy[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.cosTS_yy",tag.Data()),
			     &HGamEventInfoAuxDyn_cosTS_yy[tag],
			     &b_HGamEventInfoAuxDyn_cosTS_yy[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.Njets",tag.Data()),
			     &HGamEventInfoAuxDyn_Njets[tag],
			     &b_HGamEventInfoAuxDyn_Njets[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.cutFlow",tag.Data()), 
			     &HGamEventInfoAuxDyn_cutFlow[tag],
			     &b_HGamEventInfoAuxDyn_cutFlow[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.weightInitial",
				  tag.Data()),
			     &HGamEventInfoAuxDyn_weightInitial[tag], 
			     &b_HGamEventInfoAuxDyn_weightInitial[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.weight",tag.Data()),
			     &HGamEventInfoAuxDyn_weight[tag],
			     &b_HGamEventInfoAuxDyn_weight[tag]);
    fChain->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.TST_met",tag.Data()),
			     &HGamEventInfoAuxDyn_TST_met[tag],
			     &b_HGamEventInfoAuxDyn_TST_met[tag]);
    fChain
      ->SetBranchAddress(Form("HGamEventInfoAuxDyn%s.crossSectionBRfilterEff",
			      tag.Data()),
			 &HGamEventInfoAuxDyn_crossSectionBRfilterEff[tag],
			 &b_HGamEventInfoAuxDyn_crossSectionBRfilterEff[tag]);
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

void DMTree::DefineSystematics(std::vector<TString> sysNames) {
  m_sysNames = sysNames;
}

#endif // #ifdef DMTree_cxx
