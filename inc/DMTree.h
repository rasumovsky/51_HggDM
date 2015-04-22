//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 21 16:41:53 2015 by ROOT version 5.34/05
// from TTree CollectionTree/xAOD event tree
// found on file: sample.root
//////////////////////////////////////////////////////////

#ifndef DMTree_h
#define DMTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;
using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class DMTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   //xAOD::EventInfo_v1 *EventInfo;
   UInt_t          EventInfoAuxDyn_eventTypeBitmask;
   vector<string>  *EventInfoAuxDyn_streamTagNames;
   UInt_t          EventInfoAuxDyn_statusElement;
   vector<string>  *EventInfoAuxDyn_streamTagTypes;
   UInt_t          EventInfoAuxDyn_extendedLevel1ID;
   vector<char>    *EventInfoAuxDyn_streamTagObeysLumiblock;
   UShort_t        EventInfoAuxDyn_level1TriggerType;
   vector<unsigned short> *EventInfoAuxDyn_subEventTime;
   Float_t         EventInfoAuxDyn_actualInteractionsPerCrossing;
   //vector<ElementLink<DataVector<xAOD::EventInfo_v1> > > *EventInfoAuxDyn_subEventLink;
   Float_t         EventInfoAuxDyn_averageInteractionsPerCrossing;
   vector<unsigned short> *EventInfoAuxDyn_subEventType;
   UInt_t          EventInfoAuxDyn_pixelFlags;
   UInt_t          EventInfoAuxDyn_sctFlags;
   UInt_t          EventInfoAuxDyn_trtFlags;
   UInt_t          EventInfoAuxDyn_larFlags;
   UInt_t          EventInfoAuxDyn_tileFlags;
   UInt_t          EventInfoAuxDyn_muonFlags;
   UInt_t          EventInfoAuxDyn_forwardDetFlags;
   UInt_t          EventInfoAuxDyn_coreFlags;
   UInt_t          EventInfoAuxDyn_backgroundFlags;
   UInt_t          EventInfoAuxDyn_lumiFlags;
   Float_t         EventInfoAuxDyn_beamPosX;
   Float_t         EventInfoAuxDyn_beamPosY;
   Float_t         EventInfoAuxDyn_beamPosZ;
   Float_t         EventInfoAuxDyn_beamPosSigmaX;
   Float_t         EventInfoAuxDyn_beamPosSigmaY;
   Float_t         EventInfoAuxDyn_beamPosSigmaZ;
   Float_t         EventInfoAuxDyn_beamPosSigmaXY;
   Float_t         EventInfoAuxDyn_beamTiltXZ;
   Float_t         EventInfoAuxDyn_beamTiltYZ;
   UInt_t          EventInfoAuxDyn_beamStatus;
   vector<float>   *EventInfoAuxDyn_mcEventWeights;
   UInt_t          EventInfoAuxDyn_mcChannelNumber;
   ULong64_t       EventInfoAuxDyn_mcEventNumber;
   Float_t         EventInfoAuxDyn_m_yy;
   Float_t         EventInfoAuxDyn_pt_yy;
   Float_t         EventInfoAuxDyn_y1_pt;
   Float_t         EventInfoAuxDyn_y1_eta;
   Float_t         EventInfoAuxDyn_y2_pt;
   Float_t         EventInfoAuxDyn_y2_eta;
   Float_t         EventInfoAuxDyn_y1_ID;
   Float_t         EventInfoAuxDyn_y2_ID;
   Float_t         EventInfoAuxDyn_y1_track_iso;
   Float_t         EventInfoAuxDyn_y2_track_iso;
   Float_t         EventInfoAuxDyn_metref_final;
   Float_t         EventInfoAuxDyn_PileupWeight;
   UInt_t          EventInfoAuxDyn_runNumber;
   ULong64_t       EventInfoAuxDyn_eventNumber;
   UInt_t          EventInfoAuxDyn_lumiBlock;
   UInt_t          EventInfoAuxDyn_timeStamp;
   UInt_t          EventInfoAuxDyn_timeStampNSOffset;
   UInt_t          EventInfoAuxDyn_bcid;
   UInt_t          EventInfoAuxDyn_detectorMask0;
   UInt_t          EventInfoAuxDyn_detectorMask1;
   //vector<pair<string,string> > *EventInfoAuxDyn_detDescrTags;

   // List of branches
   TBranch        *b_EventInfo;   //!
   TBranch        *b_EventInfoAuxDyn_eventTypeBitmask;   //!
   TBranch        *b_EventInfoAuxDyn_streamTagNames;   //!
   TBranch        *b_EventInfoAuxDyn_statusElement;   //!
   TBranch        *b_EventInfoAuxDyn_streamTagTypes;   //!
   TBranch        *b_EventInfoAuxDyn_extendedLevel1ID;   //!
   TBranch        *b_EventInfoAuxDyn_streamTagObeysLumiblock;   //!
   TBranch        *b_EventInfoAuxDyn_level1TriggerType;   //!
   TBranch        *b_EventInfoAuxDyn_subEventTime;   //!
   TBranch        *b_EventInfoAuxDyn_actualInteractionsPerCrossing;   //!
   TBranch        *b_EventInfoAuxDyn_subEventLink;   //!
   TBranch        *b_EventInfoAuxDyn_averageInteractionsPerCrossing;   //!
   TBranch        *b_EventInfoAuxDyn_subEventType;   //!
   TBranch        *b_EventInfoAuxDyn_pixelFlags;   //!
   TBranch        *b_EventInfoAuxDyn_sctFlags;   //!
   TBranch        *b_EventInfoAuxDyn_trtFlags;   //!
   TBranch        *b_EventInfoAuxDyn_larFlags;   //!
   TBranch        *b_EventInfoAuxDyn_tileFlags;   //!
   TBranch        *b_EventInfoAuxDyn_muonFlags;   //!
   TBranch        *b_EventInfoAuxDyn_forwardDetFlags;   //!
   TBranch        *b_EventInfoAuxDyn_coreFlags;   //!
   TBranch        *b_EventInfoAuxDyn_backgroundFlags;   //!
   TBranch        *b_EventInfoAuxDyn_lumiFlags;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosX;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosY;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosZ;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosSigmaX;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosSigmaY;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosSigmaZ;   //!
   TBranch        *b_EventInfoAuxDyn_beamPosSigmaXY;   //!
   TBranch        *b_EventInfoAuxDyn_beamTiltXZ;   //!
   TBranch        *b_EventInfoAuxDyn_beamTiltYZ;   //!
   TBranch        *b_EventInfoAuxDyn_beamStatus;   //!
   TBranch        *b_EventInfoAuxDyn_mcEventWeights;   //!
   TBranch        *b_EventInfoAuxDyn_mcChannelNumber;   //!
   TBranch        *b_EventInfoAuxDyn_mcEventNumber;   //!
   TBranch        *b_EventInfoAuxDyn_m_yy;   //!
   TBranch        *b_EventInfoAuxDyn_pt_yy;   //!
   TBranch        *b_EventInfoAuxDyn_y1_pt;   //!
   TBranch        *b_EventInfoAuxDyn_y1_eta;   //!
   TBranch        *b_EventInfoAuxDyn_y2_pt;   //!
   TBranch        *b_EventInfoAuxDyn_y2_eta;   //!
   TBranch        *b_EventInfoAuxDyn_y1_ID;   //!
   TBranch        *b_EventInfoAuxDyn_y2_ID;   //!
   TBranch        *b_EventInfoAuxDyn_y1_track_iso;   //!
   TBranch        *b_EventInfoAuxDyn_y2_track_iso;   //!
   TBranch        *b_EventInfoAuxDyn_metref_final;   //!
   TBranch        *b_EventInfoAuxDyn_PileupWeight;   //!
   TBranch        *b_EventInfoAuxDyn_runNumber;   //!
   TBranch        *b_EventInfoAuxDyn_eventNumber;   //!
   TBranch        *b_EventInfoAuxDyn_lumiBlock;   //!
   TBranch        *b_EventInfoAuxDyn_timeStamp;   //!
   TBranch        *b_EventInfoAuxDyn_timeStampNSOffset;   //!
   TBranch        *b_EventInfoAuxDyn_bcid;   //!
   TBranch        *b_EventInfoAuxDyn_detectorMask0;   //!
   TBranch        *b_EventInfoAuxDyn_detectorMask1;   //!
   //TBranch        *b_EventInfoAuxDyn_detDescrTags;   //!

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
   //EventInfo = 0;
   EventInfoAuxDyn_streamTagNames = 0;
   EventInfoAuxDyn_streamTagTypes = 0;
   EventInfoAuxDyn_streamTagObeysLumiblock = 0;
   EventInfoAuxDyn_subEventTime = 0;
   //EventInfoAuxDyn_subEventLink = 0;
   EventInfoAuxDyn_subEventType = 0;
   EventInfoAuxDyn_mcEventWeights = 0;
   //EventInfoAuxDyn_detDescrTags = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("EventInfo", &EventInfo, &b_EventInfo);
   fChain->SetBranchAddress("EventInfoAuxDyn.eventTypeBitmask", &EventInfoAuxDyn_eventTypeBitmask, &b_EventInfoAuxDyn_eventTypeBitmask);
   fChain->SetBranchAddress("EventInfoAuxDyn.streamTagNames", &EventInfoAuxDyn_streamTagNames, &b_EventInfoAuxDyn_streamTagNames);
   fChain->SetBranchAddress("EventInfoAuxDyn.statusElement", &EventInfoAuxDyn_statusElement, &b_EventInfoAuxDyn_statusElement);
   fChain->SetBranchAddress("EventInfoAuxDyn.streamTagTypes", &EventInfoAuxDyn_streamTagTypes, &b_EventInfoAuxDyn_streamTagTypes);
   fChain->SetBranchAddress("EventInfoAuxDyn.extendedLevel1ID", &EventInfoAuxDyn_extendedLevel1ID, &b_EventInfoAuxDyn_extendedLevel1ID);
   fChain->SetBranchAddress("EventInfoAuxDyn.streamTagObeysLumiblock", &EventInfoAuxDyn_streamTagObeysLumiblock, &b_EventInfoAuxDyn_streamTagObeysLumiblock);
   fChain->SetBranchAddress("EventInfoAuxDyn.level1TriggerType", &EventInfoAuxDyn_level1TriggerType, &b_EventInfoAuxDyn_level1TriggerType);
   fChain->SetBranchAddress("EventInfoAuxDyn.subEventTime", &EventInfoAuxDyn_subEventTime, &b_EventInfoAuxDyn_subEventTime);
   fChain->SetBranchAddress("EventInfoAuxDyn.actualInteractionsPerCrossing", &EventInfoAuxDyn_actualInteractionsPerCrossing, &b_EventInfoAuxDyn_actualInteractionsPerCrossing);
   //fChain->SetBranchAddress("EventInfoAuxDyn.subEventLink", &EventInfoAuxDyn_subEventLink, &b_EventInfoAuxDyn_subEventLink);
   fChain->SetBranchAddress("EventInfoAuxDyn.averageInteractionsPerCrossing", &EventInfoAuxDyn_averageInteractionsPerCrossing, &b_EventInfoAuxDyn_averageInteractionsPerCrossing);
   fChain->SetBranchAddress("EventInfoAuxDyn.subEventType", &EventInfoAuxDyn_subEventType, &b_EventInfoAuxDyn_subEventType);
   fChain->SetBranchAddress("EventInfoAuxDyn.pixelFlags", &EventInfoAuxDyn_pixelFlags, &b_EventInfoAuxDyn_pixelFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.sctFlags", &EventInfoAuxDyn_sctFlags, &b_EventInfoAuxDyn_sctFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.trtFlags", &EventInfoAuxDyn_trtFlags, &b_EventInfoAuxDyn_trtFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.larFlags", &EventInfoAuxDyn_larFlags, &b_EventInfoAuxDyn_larFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.tileFlags", &EventInfoAuxDyn_tileFlags, &b_EventInfoAuxDyn_tileFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.muonFlags", &EventInfoAuxDyn_muonFlags, &b_EventInfoAuxDyn_muonFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.forwardDetFlags", &EventInfoAuxDyn_forwardDetFlags, &b_EventInfoAuxDyn_forwardDetFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.coreFlags", &EventInfoAuxDyn_coreFlags, &b_EventInfoAuxDyn_coreFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.backgroundFlags", &EventInfoAuxDyn_backgroundFlags, &b_EventInfoAuxDyn_backgroundFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.lumiFlags", &EventInfoAuxDyn_lumiFlags, &b_EventInfoAuxDyn_lumiFlags);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosX", &EventInfoAuxDyn_beamPosX, &b_EventInfoAuxDyn_beamPosX);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosY", &EventInfoAuxDyn_beamPosY, &b_EventInfoAuxDyn_beamPosY);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosZ", &EventInfoAuxDyn_beamPosZ, &b_EventInfoAuxDyn_beamPosZ);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosSigmaX", &EventInfoAuxDyn_beamPosSigmaX, &b_EventInfoAuxDyn_beamPosSigmaX);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosSigmaY", &EventInfoAuxDyn_beamPosSigmaY, &b_EventInfoAuxDyn_beamPosSigmaY);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosSigmaZ", &EventInfoAuxDyn_beamPosSigmaZ, &b_EventInfoAuxDyn_beamPosSigmaZ);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamPosSigmaXY", &EventInfoAuxDyn_beamPosSigmaXY, &b_EventInfoAuxDyn_beamPosSigmaXY);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamTiltXZ", &EventInfoAuxDyn_beamTiltXZ, &b_EventInfoAuxDyn_beamTiltXZ);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamTiltYZ", &EventInfoAuxDyn_beamTiltYZ, &b_EventInfoAuxDyn_beamTiltYZ);
   fChain->SetBranchAddress("EventInfoAuxDyn.beamStatus", &EventInfoAuxDyn_beamStatus, &b_EventInfoAuxDyn_beamStatus);
   fChain->SetBranchAddress("EventInfoAuxDyn.mcEventWeights", &EventInfoAuxDyn_mcEventWeights, &b_EventInfoAuxDyn_mcEventWeights);
   fChain->SetBranchAddress("EventInfoAuxDyn.mcChannelNumber", &EventInfoAuxDyn_mcChannelNumber, &b_EventInfoAuxDyn_mcChannelNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.mcEventNumber", &EventInfoAuxDyn_mcEventNumber, &b_EventInfoAuxDyn_mcEventNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.m_yy", &EventInfoAuxDyn_m_yy, &b_EventInfoAuxDyn_m_yy);
   fChain->SetBranchAddress("EventInfoAuxDyn.pt_yy", &EventInfoAuxDyn_pt_yy, &b_EventInfoAuxDyn_pt_yy);
   fChain->SetBranchAddress("EventInfoAuxDyn.y1_pt", &EventInfoAuxDyn_y1_pt, &b_EventInfoAuxDyn_y1_pt);
   fChain->SetBranchAddress("EventInfoAuxDyn.y1_eta", &EventInfoAuxDyn_y1_eta, &b_EventInfoAuxDyn_y1_eta);
   fChain->SetBranchAddress("EventInfoAuxDyn.y2_pt", &EventInfoAuxDyn_y2_pt, &b_EventInfoAuxDyn_y2_pt);
   fChain->SetBranchAddress("EventInfoAuxDyn.y2_eta", &EventInfoAuxDyn_y2_eta, &b_EventInfoAuxDyn_y2_eta);
   fChain->SetBranchAddress("EventInfoAuxDyn.y1_ID", &EventInfoAuxDyn_y1_ID, &b_EventInfoAuxDyn_y1_ID);
   fChain->SetBranchAddress("EventInfoAuxDyn.y2_ID", &EventInfoAuxDyn_y2_ID, &b_EventInfoAuxDyn_y2_ID);
   fChain->SetBranchAddress("EventInfoAuxDyn.y1_track_iso", &EventInfoAuxDyn_y1_track_iso, &b_EventInfoAuxDyn_y1_track_iso);
   fChain->SetBranchAddress("EventInfoAuxDyn.y2_track_iso", &EventInfoAuxDyn_y2_track_iso, &b_EventInfoAuxDyn_y2_track_iso);
   fChain->SetBranchAddress("EventInfoAuxDyn.metref_final", &EventInfoAuxDyn_metref_final, &b_EventInfoAuxDyn_metref_final);
   fChain->SetBranchAddress("EventInfoAuxDyn.PileupWeight", &EventInfoAuxDyn_PileupWeight, &b_EventInfoAuxDyn_PileupWeight);
   fChain->SetBranchAddress("EventInfoAuxDyn.runNumber", &EventInfoAuxDyn_runNumber, &b_EventInfoAuxDyn_runNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.eventNumber", &EventInfoAuxDyn_eventNumber, &b_EventInfoAuxDyn_eventNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.lumiBlock", &EventInfoAuxDyn_lumiBlock, &b_EventInfoAuxDyn_lumiBlock);
   fChain->SetBranchAddress("EventInfoAuxDyn.timeStamp", &EventInfoAuxDyn_timeStamp, &b_EventInfoAuxDyn_timeStamp);
   fChain->SetBranchAddress("EventInfoAuxDyn.timeStampNSOffset", &EventInfoAuxDyn_timeStampNSOffset, &b_EventInfoAuxDyn_timeStampNSOffset);
   fChain->SetBranchAddress("EventInfoAuxDyn.bcid", &EventInfoAuxDyn_bcid, &b_EventInfoAuxDyn_bcid);
   fChain->SetBranchAddress("EventInfoAuxDyn.detectorMask0", &EventInfoAuxDyn_detectorMask0, &b_EventInfoAuxDyn_detectorMask0);
   fChain->SetBranchAddress("EventInfoAuxDyn.detectorMask1", &EventInfoAuxDyn_detectorMask1, &b_EventInfoAuxDyn_detectorMask1);
   //fChain->SetBranchAddress("EventInfoAuxDyn.detDescrTags", &EventInfoAuxDyn_detDescrTags, &b_EventInfoAuxDyn_detDescrTags);
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
