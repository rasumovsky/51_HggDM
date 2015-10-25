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

// Header file for the classes stored in the TTree if any.
//#include "xAODEventInfo/versions/EventInfo_v1.h"
//#include "vector"
//#include "AthContainers/DataVector.h"
//#include "xAODMissingET/versions/MissingETContainer_v1.h"


class DMTree {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types
  //xAOD::EventInfo_v1 *EventInfo;
  /*
  Char_t          EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose;
  Char_t          EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium;
  UInt_t          EventInfoAuxDyn_runNumber;
  UInt_t          EventInfoAuxDyn_eventTypeBitmask;
  ULong64_t       EventInfoAuxDyn_eventNumber;
  UInt_t          EventInfoAuxDyn_lumiBlock;
  UInt_t          EventInfoAuxDyn_mcChannelNumber;
  Float_t         EventInfoAuxDyn_averageInteractionsPerCrossing;
  vector<float>   *EventInfoAuxDyn_mcEventWeights;
  Char_t          EventInfoAuxDyn_passTrig_HLT_2g20_tight;
  //DataVector<xAOD::Jet_v1> *HGamAntiKt4EMTopoJets;
  vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60_Eff;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_60;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_pt;
  vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_eta;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70_Eff;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_phi;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_70;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_m;
  vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77_Eff;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_77;
  vector<char>    *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85_Eff;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_85;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_Jvt;
  vector<float>   *HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta;
  //DataVector<xAOD::Jet_v1> *HGamAntiKt4TruthJets;
  vector<float>   *HGamAntiKt4TruthJetsAuxDyn_pt;
  vector<float>   *HGamAntiKt4TruthJetsAuxDyn_eta;
  vector<float>   *HGamAntiKt4TruthJetsAuxDyn_phi;
  vector<float>   *HGamAntiKt4TruthJetsAuxDyn_m;
  //DataVector<xAOD::Electron_v1> *HGamElectrons;
  */
  vector<float>   *HGamElectronsAuxDyn_pt;
  /*
  vector<float>   *HGamElectronsAuxDyn_eta;
  vector<float>   *HGamElectronsAuxDyn_phi;
  vector<float>   *HGamElectronsAuxDyn_ptvarcone20;
  vector<float>   *HGamElectronsAuxDyn_topoetcone20;
  vector<float>   *HGamElectronsAuxDyn_m;
  //xAOD::EventInfo_v1 *HGamEventInfo;
  Float_t         HGamEventInfoAuxDyn_pTt_yy;
  Float_t         HGamEventInfoAuxDyn_m_jj;
  Float_t         HGamEventInfoAuxDyn_Dy_j_j;
  Float_t         HGamEventInfoAuxDyn_Dphi_yy_jj;
  Char_t          HGamEventInfoAuxDyn_isPassedBasic;
  Char_t          HGamEventInfoAuxDyn_isPassed;
  Char_t          HGamEventInfoAuxDyn_isPassedEventClean;
  Char_t          HGamEventInfoAuxDyn_isDalitz;
  Float_t         HGamEventInfoAuxDyn_m_yy_resolution;
  Int_t           HGamEventInfoAuxDyn_NLoosePhotons;
  */
  Float_t         HGamEventInfoAuxDyn_m_yy;
  Float_t         HGamEventInfoAuxDyn_pT_yy;
  Float_t         HGamEventInfoAuxDyn_cosTS_yy;
  Float_t         HGamEventInfoAuxDyn_Njets;
  Int_t           HGamEventInfoAuxDyn_cutFlow;
  Float_t         HGamEventInfoAuxDyn_weightInitial;
  Float_t         HGamEventInfoAuxDyn_weight;
  Float_t         HGamEventInfoAuxDyn_TST_met;
  /*
  Float_t         HGamEventInfoAuxDyn_weightCategory;
  Int_t           HGamEventInfoAuxDyn_category;
  Float_t         HGamEventInfoAuxDyn_TST_sumet;
  Int_t           HGamEventInfoAuxDyn_numberOfPrimaryVertices;
  Float_t         HGamEventInfoAuxDyn_CST_met;
  Float_t         HGamEventInfoAuxDyn_selectedVertexZ;
  Float_t         HGamEventInfoAuxDyn_CST_sumet;
  Float_t         HGamEventInfoAuxDyn_hardestVertexZ;
  Char_t          HGamEventInfoAuxDyn_isPassedPreselection;
  Float_t         HGamEventInfoAuxDyn_truthVertexZ;
  Char_t          HGamEventInfoAuxDyn_isPassedPID;
  Float_t         HGamEventInfoAuxDyn_eventShapeDensity;
  Char_t          HGamEventInfoAuxDyn_isPassedIsolation;
  Float_t         HGamEventInfoAuxDyn_mu;
  Char_t          HGamEventInfoAuxDyn_isPassedRelPtCuts;
  Char_t          HGamEventInfoAuxDyn_isPassedMassCut;
  */
  Float_t         HGamEventInfoAuxDyn_crossSectionBRfilterEff;
  /*
  Float_t         HGamEventInfoAuxDyn_yAbs_yy;
  //xAOD::MissingETContainer_v1 *HGamMET_Reference_AntiKt4EMTopo;
  vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy;
  vector<ULong64_t> *HGamMET_Reference_AntiKt4EMTopoAuxDyn_source;
  vector<string>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_name;
  vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet;
  vector<double>  *HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx;
  //DataVector<xAOD::MissingET_v1> *HGamMET_Truth;
  vector<double>  *HGamMET_TruthAuxDyn_mpy;
  vector<ULong64_t> *HGamMET_TruthAuxDyn_source;
  vector<string>  *HGamMET_TruthAuxDyn_name;
  vector<double>  *HGamMET_TruthAuxDyn_sumet;
  vector<double>  *HGamMET_TruthAuxDyn_mpx;
  //DataVector<xAOD::Muon_v1> *HGamMuons;
  */
  vector<float>   *HGamMuonsAuxDyn_pt;
  /*
  vector<char>    *HGamMuonsAuxDyn_passIPCut;
  vector<float>   *HGamMuonsAuxDyn_eta;
  vector<float>   *HGamMuonsAuxDyn_phi;
  vector<float>   *HGamMuonsAuxDyn_ptvarcone20;
  vector<float>   *HGamMuonsAuxDyn_topoetcone20;
  //DataVector<xAOD::Muon_v1> *HGamMuonsInJets;
  vector<float>   *HGamMuonsInJetsAuxDyn_pt;
  vector<char>    *HGamMuonsInJetsAuxDyn_passIPCut;
  vector<float>   *HGamMuonsInJetsAuxDyn_eta;
  vector<float>   *HGamMuonsInJetsAuxDyn_phi;
  vector<float>   *HGamMuonsInJetsAuxDyn_ptvarcone20;
  vector<float>   *HGamMuonsInJetsAuxDyn_topoetcone20;
  //DataVector<xAOD::Photon_v1> *HGamPhotons;
  vector<char>    *HGamPhotonsAuxDyn_isIsoCone40;
  vector<float>   *HGamPhotonsAuxDyn_topoetcone20;
  vector<unsigned int> *HGamPhotonsAuxDyn_isEMTight;
  vector<char>    *HGamPhotonsAuxDyn_isTight;
  vector<char>    *HGamPhotonsAuxDyn_isIsoCone40CaloOnly;
  vector<char>    *HGamPhotonsAuxDyn_isConv;
  vector<int>     *HGamPhotonsAuxDyn_truthOrigin;
  vector<int>     *HGamPhotonsAuxDyn_truthType;
  vector<char>    *HGamPhotonsAuxDyn_isTight_nofudge;
  vector<unsigned int> *HGamPhotonsAuxDyn_isEMTight_nofudge;
  vector<float>   *HGamPhotonsAuxDyn_ptcone20;
  vector<float>   *HGamPhotonsAuxDyn_eta_s2;
  vector<float>   *HGamPhotonsAuxDyn_pt;
  vector<float>   *HGamPhotonsAuxDyn_eta;
  vector<float>   *HGamPhotonsAuxDyn_phi;
  vector<float>   *HGamPhotonsAuxDyn_m;
  vector<float>   *HGamPhotonsAuxDyn_relEreso;
  vector<char>    *HGamPhotonsAuxDyn_isIsoCone20;
  vector<char>    *HGamPhotonsAuxDyn_isIsoCone20Higgs;
  //DataVector<xAOD::TruthParticle_v1> *HGamTruthElectrons;
  vector<float>   *HGamTruthElectronsAuxDyn_pt;
  vector<float>   *HGamTruthElectronsAuxDyn_eta;
  vector<float>   *HGamTruthElectronsAuxDyn_m;
  vector<float>   *HGamTruthElectronsAuxDyn_px;
  vector<float>   *HGamTruthElectronsAuxDyn_py;
  vector<float>   *HGamTruthElectronsAuxDyn_pz;
  vector<float>   *HGamTruthElectronsAuxDyn_e;
  //xAOD::EventInfo_v1 *HGamTruthEventInfo;
   Float_t         HGamTruthEventInfoAuxDyn_m_h1;
   Float_t         HGamTruthEventInfoAuxDyn_yAbs_yy;
   Float_t         HGamTruthEventInfoAuxDyn_m_h2;
   Float_t         HGamTruthEventInfoAuxDyn_pTt_yy;
   Float_t         HGamTruthEventInfoAuxDyn_m_yy;
   Float_t         HGamTruthEventInfoAuxDyn_pT_yy;
   Float_t         HGamTruthEventInfoAuxDyn_cosTS_yy;
   Float_t         HGamTruthEventInfoAuxDyn_Njets;
   Float_t         HGamTruthEventInfoAuxDyn_m_jj;
   Float_t         HGamTruthEventInfoAuxDyn_Dy_j_j;
   Float_t         HGamTruthEventInfoAuxDyn_Dphi_yy_jj;
   Char_t          HGamTruthEventInfoAuxDyn_isFiducial;
   Char_t          HGamTruthEventInfoAuxDyn_isFiducialKinOnly;
   Float_t         HGamTruthEventInfoAuxDyn_TruthNonInt_met;
   Float_t         HGamTruthEventInfoAuxDyn_TruthInt_sumet;
   Float_t         HGamTruthEventInfoAuxDyn_pT_h1;
   Float_t         HGamTruthEventInfoAuxDyn_pT_h2;
   Float_t         HGamTruthEventInfoAuxDyn_y_h1;
   Float_t         HGamTruthEventInfoAuxDyn_y_h2;
   //DataVector<xAOD::TruthParticle_v1> *HGamTruthHiggsBosons;
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_pt;
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_eta;
  */
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_m;
   /*
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_px;
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_py;
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_pz;
   vector<float>   *HGamTruthHiggsBosonsAuxDyn_e;
   //DataVector<xAOD::TruthParticle_v1> *HGamTruthMuons;
   vector<float>   *HGamTruthMuonsAuxDyn_pt;
   vector<float>   *HGamTruthMuonsAuxDyn_eta;
   vector<float>   *HGamTruthMuonsAuxDyn_m;
   vector<float>   *HGamTruthMuonsAuxDyn_px;
   vector<float>   *HGamTruthMuonsAuxDyn_py;
   vector<float>   *HGamTruthMuonsAuxDyn_pz;
   vector<float>   *HGamTruthMuonsAuxDyn_e;
   //DataVector<xAOD::TruthParticle_v1> *HGamTruthPhotons;
   vector<char>    *HGamTruthPhotonsAuxDyn_isIsolated;
   vector<float>   *HGamTruthPhotonsAuxDyn_etcone20;
   vector<float>   *HGamTruthPhotonsAuxDyn_etcone40;
   vector<int>     *HGamTruthPhotonsAuxDyn_truthOrigin;
   vector<int>     *HGamTruthPhotonsAuxDyn_truthType;
   vector<float>   *HGamTruthPhotonsAuxDyn_pt;
   vector<float>   *HGamTruthPhotonsAuxDyn_px;
   vector<float>   *HGamTruthPhotonsAuxDyn_eta;
   vector<float>   *HGamTruthPhotonsAuxDyn_py;
   vector<float>   *HGamTruthPhotonsAuxDyn_m;
   vector<float>   *HGamTruthPhotonsAuxDyn_pz;
   vector<float>   *HGamTruthPhotonsAuxDyn_e;
  */
  
   // List of branches
  /*
  TBranch        *b_EventInfo;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium;   //!
   TBranch        *b_EventInfoAuxDyn_runNumber;   //!
   TBranch        *b_EventInfoAuxDyn_eventTypeBitmask;   //!
   TBranch        *b_EventInfoAuxDyn_eventNumber;   //!
   TBranch        *b_EventInfoAuxDyn_lumiBlock;   //!
   TBranch        *b_EventInfoAuxDyn_mcChannelNumber;   //!
   TBranch        *b_EventInfoAuxDyn_averageInteractionsPerCrossing;   //!
   TBranch        *b_EventInfoAuxDyn_mcEventWeights;   //!
   TBranch        *b_EventInfoAuxDyn_passTrig_HLT_2g20_tight;   //!
   TBranch        *b_HGamAntiKt4EMTopoJets;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60_Eff;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_60;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_pt;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_eta;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70_Eff;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_phi;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_70;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_m;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77_Eff;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_77;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85_Eff;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_85;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_Jvt;   //!
   TBranch        *b_HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta;   //!
   TBranch        *b_HGamAntiKt4TruthJets;   //!
   TBranch        *b_HGamAntiKt4TruthJetsAuxDyn_pt;   //!
   TBranch        *b_HGamAntiKt4TruthJetsAuxDyn_eta;   //!
   TBranch        *b_HGamAntiKt4TruthJetsAuxDyn_phi;   //!
   TBranch        *b_HGamAntiKt4TruthJetsAuxDyn_m;   //!
   TBranch        *b_HGamElectrons;   //!
  */
   TBranch        *b_HGamElectronsAuxDyn_pt;   //!
   /*
   TBranch        *b_HGamElectronsAuxDyn_eta;   //!
   TBranch        *b_HGamElectronsAuxDyn_phi;   //!
   TBranch        *b_HGamElectronsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamElectronsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamElectronsAuxDyn_m;   //!
   TBranch        *b_HGamEventInfo;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pTt_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_jj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dy_j_j;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Dphi_yy_jj;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedBasic;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassed;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedEventClean;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isDalitz;   //!
   TBranch        *b_HGamEventInfoAuxDyn_m_yy_resolution;   //!
   TBranch        *b_HGamEventInfoAuxDyn_NLoosePhotons;   //!
   */
   TBranch        *b_HGamEventInfoAuxDyn_m_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_pT_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_cosTS_yy;   //!
   TBranch        *b_HGamEventInfoAuxDyn_Njets;   //!
   TBranch        *b_HGamEventInfoAuxDyn_cutFlow;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weightInitial;   //!
   TBranch        *b_HGamEventInfoAuxDyn_weight;   //!
   TBranch        *b_HGamEventInfoAuxDyn_TST_met;   //!
   /*
   TBranch        *b_HGamEventInfoAuxDyn_weightCategory;   //!
   TBranch        *b_HGamEventInfoAuxDyn_category;   //!
   TBranch        *b_HGamEventInfoAuxDyn_TST_sumet;   //!
   TBranch        *b_HGamEventInfoAuxDyn_numberOfPrimaryVertices;   //!
   TBranch        *b_HGamEventInfoAuxDyn_CST_met;   //!
   TBranch        *b_HGamEventInfoAuxDyn_selectedVertexZ;   //!
   TBranch        *b_HGamEventInfoAuxDyn_CST_sumet;   //!
   TBranch        *b_HGamEventInfoAuxDyn_hardestVertexZ;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedPreselection;   //!
   TBranch        *b_HGamEventInfoAuxDyn_truthVertexZ;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedPID;   //!
   TBranch        *b_HGamEventInfoAuxDyn_eventShapeDensity;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedIsolation;   //!
   TBranch        *b_HGamEventInfoAuxDyn_mu;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedRelPtCuts;   //!
   TBranch        *b_HGamEventInfoAuxDyn_isPassedMassCut;   //!
   */
   TBranch        *b_HGamEventInfoAuxDyn_crossSectionBRfilterEff;   //!
   /*
   TBranch        *b_HGamEventInfoAuxDyn_yAbs_yy;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopo;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_source;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_name;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet;   //!
   TBranch        *b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx;   //!
   TBranch        *b_HGamMET_Truth;   //!
   TBranch        *b_HGamMET_TruthAuxDyn_mpy;   //!
   TBranch        *b_HGamMET_TruthAuxDyn_source;   //!
   TBranch        *b_HGamMET_TruthAuxDyn_name;   //!
   TBranch        *b_HGamMET_TruthAuxDyn_sumet;   //!
   TBranch        *b_HGamMET_TruthAuxDyn_mpx;   //!
   TBranch        *b_HGamMuons;   //!
  */
   TBranch        *b_HGamMuonsAuxDyn_pt;   //!
   /*
   TBranch        *b_HGamMuonsAuxDyn_passIPCut;   //!
   TBranch        *b_HGamMuonsAuxDyn_eta;   //!
   TBranch        *b_HGamMuonsAuxDyn_phi;   //!
   TBranch        *b_HGamMuonsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamMuonsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamMuonsInJets;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_pt;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_passIPCut;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_eta;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_phi;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_ptvarcone20;   //!
   TBranch        *b_HGamMuonsInJetsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamPhotons;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoCone40;   //!
   TBranch        *b_HGamPhotonsAuxDyn_topoetcone20;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isEMTight;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isTight;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoCone40CaloOnly;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isConv;   //!
   TBranch        *b_HGamPhotonsAuxDyn_truthOrigin;   //!
   TBranch        *b_HGamPhotonsAuxDyn_truthType;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isTight_nofudge;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isEMTight_nofudge;   //!
   TBranch        *b_HGamPhotonsAuxDyn_ptcone20;   //!
   TBranch        *b_HGamPhotonsAuxDyn_eta_s2;   //!
   TBranch        *b_HGamPhotonsAuxDyn_pt;   //!
   TBranch        *b_HGamPhotonsAuxDyn_eta;   //!
   TBranch        *b_HGamPhotonsAuxDyn_phi;   //!
   TBranch        *b_HGamPhotonsAuxDyn_m;   //!
   TBranch        *b_HGamPhotonsAuxDyn_relEreso;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoCone20;   //!
   TBranch        *b_HGamPhotonsAuxDyn_isIsoCone20Higgs;   //!
   TBranch        *b_HGamTruthElectrons;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_pt;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_eta;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_m;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_px;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_py;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_pz;   //!
   TBranch        *b_HGamTruthElectronsAuxDyn_e;   //!
   TBranch        *b_HGamTruthEventInfo;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_m_h1;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_yAbs_yy;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_m_h2;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_pTt_yy;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_m_yy;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_pT_yy;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_cosTS_yy;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_Njets;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_m_jj;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_Dy_j_j;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_Dphi_yy_jj;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_isFiducial;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_isFiducialKinOnly;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_TruthNonInt_met;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_TruthInt_sumet;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_pT_h1;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_pT_h2;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_y_h1;   //!
   TBranch        *b_HGamTruthEventInfoAuxDyn_y_h2;   //!
   TBranch        *b_HGamTruthHiggsBosons;   //!
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_pt;   //!
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_eta;   //!
   */
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_m;   //!
   /*
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_px;   //!
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_py;   //!
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_pz;   //!
   TBranch        *b_HGamTruthHiggsBosonsAuxDyn_e;   //!
   TBranch        *b_HGamTruthMuons;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_pt;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_eta;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_m;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_px;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_py;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_pz;   //!
   TBranch        *b_HGamTruthMuonsAuxDyn_e;   //!
   TBranch        *b_HGamTruthPhotons;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_isIsolated;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_etcone20;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_etcone40;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_truthOrigin;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_truthType;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_pt;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_px;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_eta;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_py;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_m;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_pz;   //!
   TBranch        *b_HGamTruthPhotonsAuxDyn_e;   //!
  */

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("PowhegPy8_VBF125.MxAOD.p2421.h008.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("PowhegPy8_VBF125.MxAOD.p2421.h008.root");
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
   /*
   //EventInfo = 0;
   EventInfoAuxDyn_mcEventWeights = 0;
   // HGamAntiKt4EMTopoJets = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60_Eff = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_60 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_pt = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_eta = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70_Eff = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_phi = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_70 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_m = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77_Eff = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_77 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85_Eff = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_85 = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_Jvt = 0;
   HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta = 0;
   //HGamAntiKt4TruthJets = 0;
   HGamAntiKt4TruthJetsAuxDyn_pt = 0;
   HGamAntiKt4TruthJetsAuxDyn_eta = 0;
   HGamAntiKt4TruthJetsAuxDyn_phi = 0;
   HGamAntiKt4TruthJetsAuxDyn_m = 0;
   //HGamElectrons = 0;
   */
   HGamElectronsAuxDyn_pt = 0;
   /*
   HGamElectronsAuxDyn_eta = 0;
   HGamElectronsAuxDyn_phi = 0;
   HGamElectronsAuxDyn_ptvarcone20 = 0;
   HGamElectronsAuxDyn_topoetcone20 = 0;
   HGamElectronsAuxDyn_m = 0;
   //HGamEventInfo = 0;
   //HGamMET_Reference_AntiKt4EMTopo = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_source = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_name = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet = 0;
   HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx = 0;
   //HGamMET_Truth = 0;
   HGamMET_TruthAuxDyn_mpy = 0;
   HGamMET_TruthAuxDyn_source = 0;
   HGamMET_TruthAuxDyn_name = 0;
   HGamMET_TruthAuxDyn_sumet = 0;
   HGamMET_TruthAuxDyn_mpx = 0;
   //HGamMuons = 0;
   */
   HGamMuonsAuxDyn_pt = 0;
   /*
   HGamMuonsAuxDyn_passIPCut = 0;
   HGamMuonsAuxDyn_eta = 0;
   HGamMuonsAuxDyn_phi = 0;
   HGamMuonsAuxDyn_ptvarcone20 = 0;
   HGamMuonsAuxDyn_topoetcone20 = 0;
   //HGamMuonsInJets = 0;
   HGamMuonsInJetsAuxDyn_pt = 0;
   HGamMuonsInJetsAuxDyn_passIPCut = 0;
   HGamMuonsInJetsAuxDyn_eta = 0;
   HGamMuonsInJetsAuxDyn_phi = 0;
   HGamMuonsInJetsAuxDyn_ptvarcone20 = 0;
   HGamMuonsInJetsAuxDyn_topoetcone20 = 0;
   //HGamPhotons = 0;
   HGamPhotonsAuxDyn_isIsoCone40 = 0;
   HGamPhotonsAuxDyn_topoetcone20 = 0;
   HGamPhotonsAuxDyn_isEMTight = 0;
   HGamPhotonsAuxDyn_isTight = 0;
   HGamPhotonsAuxDyn_isIsoCone40CaloOnly = 0;
   HGamPhotonsAuxDyn_isConv = 0;
   HGamPhotonsAuxDyn_truthOrigin = 0;
   HGamPhotonsAuxDyn_truthType = 0;
   HGamPhotonsAuxDyn_isTight_nofudge = 0;
   HGamPhotonsAuxDyn_isEMTight_nofudge = 0;
   HGamPhotonsAuxDyn_ptcone20 = 0;
   HGamPhotonsAuxDyn_eta_s2 = 0;
   HGamPhotonsAuxDyn_pt = 0;
   HGamPhotonsAuxDyn_eta = 0;
   HGamPhotonsAuxDyn_phi = 0;
   HGamPhotonsAuxDyn_m = 0;
   HGamPhotonsAuxDyn_relEreso = 0;
   HGamPhotonsAuxDyn_isIsoCone20 = 0;
   HGamPhotonsAuxDyn_isIsoCone20Higgs = 0;
   //HGamTruthElectrons = 0;
   HGamTruthElectronsAuxDyn_pt = 0;
   HGamTruthElectronsAuxDyn_eta = 0;
   HGamTruthElectronsAuxDyn_m = 0;
   HGamTruthElectronsAuxDyn_px = 0;
   HGamTruthElectronsAuxDyn_py = 0;
   HGamTruthElectronsAuxDyn_pz = 0;
   HGamTruthElectronsAuxDyn_e = 0;
   //HGamTruthEventInfo = 0;
   //HGamTruthHiggsBosons = 0;
   HGamTruthHiggsBosonsAuxDyn_pt = 0;
   HGamTruthHiggsBosonsAuxDyn_eta = 0;
   */
   HGamTruthHiggsBosonsAuxDyn_m = 0;
   /*
   HGamTruthHiggsBosonsAuxDyn_px = 0;
   HGamTruthHiggsBosonsAuxDyn_py = 0;
   HGamTruthHiggsBosonsAuxDyn_pz = 0;
   HGamTruthHiggsBosonsAuxDyn_e = 0;
   //HGamTruthMuons = 0;
   HGamTruthMuonsAuxDyn_pt = 0;
   HGamTruthMuonsAuxDyn_eta = 0;
   HGamTruthMuonsAuxDyn_m = 0;
   HGamTruthMuonsAuxDyn_px = 0;
   HGamTruthMuonsAuxDyn_py = 0;
   HGamTruthMuonsAuxDyn_pz = 0;
   HGamTruthMuonsAuxDyn_e = 0;
   //HGamTruthPhotons = 0;
   HGamTruthPhotonsAuxDyn_isIsolated = 0;
   HGamTruthPhotonsAuxDyn_etcone20 = 0;
   HGamTruthPhotonsAuxDyn_etcone40 = 0;
   HGamTruthPhotonsAuxDyn_truthOrigin = 0;
   HGamTruthPhotonsAuxDyn_truthType = 0;
   HGamTruthPhotonsAuxDyn_pt = 0;
   HGamTruthPhotonsAuxDyn_px = 0;
   HGamTruthPhotonsAuxDyn_eta = 0;
   HGamTruthPhotonsAuxDyn_py = 0;
   HGamTruthPhotonsAuxDyn_m = 0;
   HGamTruthPhotonsAuxDyn_pz = 0;
   HGamTruthPhotonsAuxDyn_e = 0;
   */

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   
   /*
   //fChain->SetBranchAddress("EventInfo", &EventInfo, &b_EventInfo);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_g35_loose_g25_loose", &EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose, &b_EventInfoAuxDyn_passTrig_HLT_g35_loose_g25_loose);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_g35_medium_g25_medium", &EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium, &b_EventInfoAuxDyn_passTrig_HLT_g35_medium_g25_medium);
   fChain->SetBranchAddress("EventInfoAuxDyn.runNumber", &EventInfoAuxDyn_runNumber, &b_EventInfoAuxDyn_runNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.eventTypeBitmask", &EventInfoAuxDyn_eventTypeBitmask, &b_EventInfoAuxDyn_eventTypeBitmask);
   fChain->SetBranchAddress("EventInfoAuxDyn.eventNumber", &EventInfoAuxDyn_eventNumber, &b_EventInfoAuxDyn_eventNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.lumiBlock", &EventInfoAuxDyn_lumiBlock, &b_EventInfoAuxDyn_lumiBlock);
   fChain->SetBranchAddress("EventInfoAuxDyn.mcChannelNumber", &EventInfoAuxDyn_mcChannelNumber, &b_EventInfoAuxDyn_mcChannelNumber);
   fChain->SetBranchAddress("EventInfoAuxDyn.averageInteractionsPerCrossing", &EventInfoAuxDyn_averageInteractionsPerCrossing, &b_EventInfoAuxDyn_averageInteractionsPerCrossing);
   fChain->SetBranchAddress("EventInfoAuxDyn.mcEventWeights", &EventInfoAuxDyn_mcEventWeights, &b_EventInfoAuxDyn_mcEventWeights);
   fChain->SetBranchAddress("EventInfoAuxDyn.passTrig_HLT_2g20_tight", &EventInfoAuxDyn_passTrig_HLT_2g20_tight, &b_EventInfoAuxDyn_passTrig_HLT_2g20_tight);
   //fChain->SetBranchAddress("HGamAntiKt4EMTopoJets", &HGamAntiKt4EMTopoJets, &b_HGamAntiKt4EMTopoJets);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_60", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_60_Eff", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60_Eff, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_60_Eff);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.SF_MV2c20_60", &HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_60, &b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_60);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.pt", &HGamAntiKt4EMTopoJetsAuxDyn_pt, &b_HGamAntiKt4EMTopoJetsAuxDyn_pt);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_70", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.eta", &HGamAntiKt4EMTopoJetsAuxDyn_eta, &b_HGamAntiKt4EMTopoJetsAuxDyn_eta);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_70_Eff", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70_Eff, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_70_Eff);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.phi", &HGamAntiKt4EMTopoJetsAuxDyn_phi, &b_HGamAntiKt4EMTopoJetsAuxDyn_phi);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.SF_MV2c20_70", &HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_70, &b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_70);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.m", &HGamAntiKt4EMTopoJetsAuxDyn_m, &b_HGamAntiKt4EMTopoJetsAuxDyn_m);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_77", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_77_Eff", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77_Eff, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_77_Eff);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.SF_MV2c20_77", &HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_77, &b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_77);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_85", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.MV2c20_85_Eff", &HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85_Eff, &b_HGamAntiKt4EMTopoJetsAuxDyn_MV2c20_85_Eff);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.SF_MV2c20_85", &HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_85, &b_HGamAntiKt4EMTopoJetsAuxDyn_SF_MV2c20_85);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.Jvt", &HGamAntiKt4EMTopoJetsAuxDyn_Jvt, &b_HGamAntiKt4EMTopoJetsAuxDyn_Jvt);
   fChain->SetBranchAddress("HGamAntiKt4EMTopoJetsAuxDyn.DetectorEta", &HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta, &b_HGamAntiKt4EMTopoJetsAuxDyn_DetectorEta);
   //fChain->SetBranchAddress("HGamAntiKt4TruthJets", &HGamAntiKt4TruthJets, &b_HGamAntiKt4TruthJets);
   fChain->SetBranchAddress("HGamAntiKt4TruthJetsAuxDyn.pt", &HGamAntiKt4TruthJetsAuxDyn_pt, &b_HGamAntiKt4TruthJetsAuxDyn_pt);
   fChain->SetBranchAddress("HGamAntiKt4TruthJetsAuxDyn.eta", &HGamAntiKt4TruthJetsAuxDyn_eta, &b_HGamAntiKt4TruthJetsAuxDyn_eta);
   fChain->SetBranchAddress("HGamAntiKt4TruthJetsAuxDyn.phi", &HGamAntiKt4TruthJetsAuxDyn_phi, &b_HGamAntiKt4TruthJetsAuxDyn_phi);
   fChain->SetBranchAddress("HGamAntiKt4TruthJetsAuxDyn.m", &HGamAntiKt4TruthJetsAuxDyn_m, &b_HGamAntiKt4TruthJetsAuxDyn_m);
   //fChain->SetBranchAddress("HGamElectrons", &HGamElectrons, &b_HGamElectrons);
   */
   fChain->SetBranchAddress("HGamElectronsAuxDyn.pt", &HGamElectronsAuxDyn_pt, &b_HGamElectronsAuxDyn_pt);
   /*
   fChain->SetBranchAddress("HGamElectronsAuxDyn.eta", &HGamElectronsAuxDyn_eta, &b_HGamElectronsAuxDyn_eta);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.phi", &HGamElectronsAuxDyn_phi, &b_HGamElectronsAuxDyn_phi);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.ptvarcone20", &HGamElectronsAuxDyn_ptvarcone20, &b_HGamElectronsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.topoetcone20", &HGamElectronsAuxDyn_topoetcone20, &b_HGamElectronsAuxDyn_topoetcone20);
   fChain->SetBranchAddress("HGamElectronsAuxDyn.m", &HGamElectronsAuxDyn_m, &b_HGamElectronsAuxDyn_m);
   //fChain->SetBranchAddress("HGamEventInfo", &HGamEventInfo, &b_HGamEventInfo);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pTt_yy", &HGamEventInfoAuxDyn_pTt_yy, &b_HGamEventInfoAuxDyn_pTt_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_jj", &HGamEventInfoAuxDyn_m_jj, &b_HGamEventInfoAuxDyn_m_jj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dy_j_j", &HGamEventInfoAuxDyn_Dy_j_j, &b_HGamEventInfoAuxDyn_Dy_j_j);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Dphi_yy_jj", &HGamEventInfoAuxDyn_Dphi_yy_jj, &b_HGamEventInfoAuxDyn_Dphi_yy_jj);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedBasic", &HGamEventInfoAuxDyn_isPassedBasic, &b_HGamEventInfoAuxDyn_isPassedBasic);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassed", &HGamEventInfoAuxDyn_isPassed, &b_HGamEventInfoAuxDyn_isPassed);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedEventClean", &HGamEventInfoAuxDyn_isPassedEventClean, &b_HGamEventInfoAuxDyn_isPassedEventClean);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isDalitz", &HGamEventInfoAuxDyn_isDalitz, &b_HGamEventInfoAuxDyn_isDalitz);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy_resolution", &HGamEventInfoAuxDyn_m_yy_resolution, &b_HGamEventInfoAuxDyn_m_yy_resolution);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.NLoosePhotons", &HGamEventInfoAuxDyn_NLoosePhotons, &b_HGamEventInfoAuxDyn_NLoosePhotons);
   */
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.m_yy", &HGamEventInfoAuxDyn_m_yy, &b_HGamEventInfoAuxDyn_m_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.pT_yy", &HGamEventInfoAuxDyn_pT_yy, &b_HGamEventInfoAuxDyn_pT_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.cosTS_yy", &HGamEventInfoAuxDyn_cosTS_yy, &b_HGamEventInfoAuxDyn_cosTS_yy);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.Njets", &HGamEventInfoAuxDyn_Njets, &b_HGamEventInfoAuxDyn_Njets);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.cutFlow", &HGamEventInfoAuxDyn_cutFlow, &b_HGamEventInfoAuxDyn_cutFlow);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightInitial", &HGamEventInfoAuxDyn_weightInitial, &b_HGamEventInfoAuxDyn_weightInitial);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weight", &HGamEventInfoAuxDyn_weight, &b_HGamEventInfoAuxDyn_weight);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.TST_met", &HGamEventInfoAuxDyn_TST_met, &b_HGamEventInfoAuxDyn_TST_met);
   /*
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.weightCategory", &HGamEventInfoAuxDyn_weightCategory, &b_HGamEventInfoAuxDyn_weightCategory);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.category", &HGamEventInfoAuxDyn_category, &b_HGamEventInfoAuxDyn_category);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.TST_sumet", &HGamEventInfoAuxDyn_TST_sumet, &b_HGamEventInfoAuxDyn_TST_sumet);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.numberOfPrimaryVertices", &HGamEventInfoAuxDyn_numberOfPrimaryVertices, &b_HGamEventInfoAuxDyn_numberOfPrimaryVertices);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.CST_met", &HGamEventInfoAuxDyn_CST_met, &b_HGamEventInfoAuxDyn_CST_met);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.selectedVertexZ", &HGamEventInfoAuxDyn_selectedVertexZ, &b_HGamEventInfoAuxDyn_selectedVertexZ);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.CST_sumet", &HGamEventInfoAuxDyn_CST_sumet, &b_HGamEventInfoAuxDyn_CST_sumet);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.hardestVertexZ", &HGamEventInfoAuxDyn_hardestVertexZ, &b_HGamEventInfoAuxDyn_hardestVertexZ);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPreselection", &HGamEventInfoAuxDyn_isPassedPreselection, &b_HGamEventInfoAuxDyn_isPassedPreselection);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.truthVertexZ", &HGamEventInfoAuxDyn_truthVertexZ, &b_HGamEventInfoAuxDyn_truthVertexZ);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedPID", &HGamEventInfoAuxDyn_isPassedPID, &b_HGamEventInfoAuxDyn_isPassedPID);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.eventShapeDensity", &HGamEventInfoAuxDyn_eventShapeDensity, &b_HGamEventInfoAuxDyn_eventShapeDensity);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedIsolation", &HGamEventInfoAuxDyn_isPassedIsolation, &b_HGamEventInfoAuxDyn_isPassedIsolation);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.mu", &HGamEventInfoAuxDyn_mu, &b_HGamEventInfoAuxDyn_mu);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedRelPtCuts", &HGamEventInfoAuxDyn_isPassedRelPtCuts, &b_HGamEventInfoAuxDyn_isPassedRelPtCuts);
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.isPassedMassCut", &HGamEventInfoAuxDyn_isPassedMassCut, &b_HGamEventInfoAuxDyn_isPassedMassCut);
   */
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff", &HGamEventInfoAuxDyn_crossSectionBRfilterEff, &b_HGamEventInfoAuxDyn_crossSectionBRfilterEff);
   /*
   fChain->SetBranchAddress("HGamEventInfoAuxDyn.yAbs_yy", &HGamEventInfoAuxDyn_yAbs_yy, &b_HGamEventInfoAuxDyn_yAbs_yy);
   //fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopo", &HGamMET_Reference_AntiKt4EMTopo, &b_HGamMET_Reference_AntiKt4EMTopo);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.mpy", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpy);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.source", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_source, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_source);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.name", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_name, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_name);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.sumet", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_sumet);
   fChain->SetBranchAddress("HGamMET_Reference_AntiKt4EMTopoAuxDyn.mpx", &HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx, &b_HGamMET_Reference_AntiKt4EMTopoAuxDyn_mpx);
   //fChain->SetBranchAddress("HGamMET_Truth", &HGamMET_Truth, &b_HGamMET_Truth);
   fChain->SetBranchAddress("HGamMET_TruthAuxDyn.mpy", &HGamMET_TruthAuxDyn_mpy, &b_HGamMET_TruthAuxDyn_mpy);
   fChain->SetBranchAddress("HGamMET_TruthAuxDyn.source", &HGamMET_TruthAuxDyn_source, &b_HGamMET_TruthAuxDyn_source);
   fChain->SetBranchAddress("HGamMET_TruthAuxDyn.name", &HGamMET_TruthAuxDyn_name, &b_HGamMET_TruthAuxDyn_name);
   fChain->SetBranchAddress("HGamMET_TruthAuxDyn.sumet", &HGamMET_TruthAuxDyn_sumet, &b_HGamMET_TruthAuxDyn_sumet);
   fChain->SetBranchAddress("HGamMET_TruthAuxDyn.mpx", &HGamMET_TruthAuxDyn_mpx, &b_HGamMET_TruthAuxDyn_mpx);
   //fChain->SetBranchAddress("HGamMuons", &HGamMuons, &b_HGamMuons);
   */
   fChain->SetBranchAddress("HGamMuonsAuxDyn.pt", &HGamMuonsAuxDyn_pt, &b_HGamMuonsAuxDyn_pt);
   /*
   fChain->SetBranchAddress("HGamMuonsAuxDyn.passIPCut", &HGamMuonsAuxDyn_passIPCut, &b_HGamMuonsAuxDyn_passIPCut);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.eta", &HGamMuonsAuxDyn_eta, &b_HGamMuonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.phi", &HGamMuonsAuxDyn_phi, &b_HGamMuonsAuxDyn_phi);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.ptvarcone20", &HGamMuonsAuxDyn_ptvarcone20, &b_HGamMuonsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamMuonsAuxDyn.topoetcone20", &HGamMuonsAuxDyn_topoetcone20, &b_HGamMuonsAuxDyn_topoetcone20);
   //fChain->SetBranchAddress("HGamMuonsInJets", &HGamMuonsInJets, &b_HGamMuonsInJets);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.pt", &HGamMuonsInJetsAuxDyn_pt, &b_HGamMuonsInJetsAuxDyn_pt);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.passIPCut", &HGamMuonsInJetsAuxDyn_passIPCut, &b_HGamMuonsInJetsAuxDyn_passIPCut);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.eta", &HGamMuonsInJetsAuxDyn_eta, &b_HGamMuonsInJetsAuxDyn_eta);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.phi", &HGamMuonsInJetsAuxDyn_phi, &b_HGamMuonsInJetsAuxDyn_phi);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.ptvarcone20", &HGamMuonsInJetsAuxDyn_ptvarcone20, &b_HGamMuonsInJetsAuxDyn_ptvarcone20);
   fChain->SetBranchAddress("HGamMuonsInJetsAuxDyn.topoetcone20", &HGamMuonsInJetsAuxDyn_topoetcone20, &b_HGamMuonsInJetsAuxDyn_topoetcone20);
   //fChain->SetBranchAddress("HGamPhotons", &HGamPhotons, &b_HGamPhotons);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoCone40", &HGamPhotonsAuxDyn_isIsoCone40, &b_HGamPhotonsAuxDyn_isIsoCone40);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.topoetcone20", &HGamPhotonsAuxDyn_topoetcone20, &b_HGamPhotonsAuxDyn_topoetcone20);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isEMTight", &HGamPhotonsAuxDyn_isEMTight, &b_HGamPhotonsAuxDyn_isEMTight);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isTight", &HGamPhotonsAuxDyn_isTight, &b_HGamPhotonsAuxDyn_isTight);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoCone40CaloOnly", &HGamPhotonsAuxDyn_isIsoCone40CaloOnly, &b_HGamPhotonsAuxDyn_isIsoCone40CaloOnly);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isConv", &HGamPhotonsAuxDyn_isConv, &b_HGamPhotonsAuxDyn_isConv);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.truthOrigin", &HGamPhotonsAuxDyn_truthOrigin, &b_HGamPhotonsAuxDyn_truthOrigin);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.truthType", &HGamPhotonsAuxDyn_truthType, &b_HGamPhotonsAuxDyn_truthType);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isTight_nofudge", &HGamPhotonsAuxDyn_isTight_nofudge, &b_HGamPhotonsAuxDyn_isTight_nofudge);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isEMTight_nofudge", &HGamPhotonsAuxDyn_isEMTight_nofudge, &b_HGamPhotonsAuxDyn_isEMTight_nofudge);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.ptcone20", &HGamPhotonsAuxDyn_ptcone20, &b_HGamPhotonsAuxDyn_ptcone20);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.eta_s2", &HGamPhotonsAuxDyn_eta_s2, &b_HGamPhotonsAuxDyn_eta_s2);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.pt", &HGamPhotonsAuxDyn_pt, &b_HGamPhotonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.eta", &HGamPhotonsAuxDyn_eta, &b_HGamPhotonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.phi", &HGamPhotonsAuxDyn_phi, &b_HGamPhotonsAuxDyn_phi);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.m", &HGamPhotonsAuxDyn_m, &b_HGamPhotonsAuxDyn_m);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.relEreso", &HGamPhotonsAuxDyn_relEreso, &b_HGamPhotonsAuxDyn_relEreso);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoCone20", &HGamPhotonsAuxDyn_isIsoCone20, &b_HGamPhotonsAuxDyn_isIsoCone20);
   fChain->SetBranchAddress("HGamPhotonsAuxDyn.isIsoCone20Higgs", &HGamPhotonsAuxDyn_isIsoCone20Higgs, &b_HGamPhotonsAuxDyn_isIsoCone20Higgs);
   //fChain->SetBranchAddress("HGamTruthElectrons", &HGamTruthElectrons, &b_HGamTruthElectrons);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.pt", &HGamTruthElectronsAuxDyn_pt, &b_HGamTruthElectronsAuxDyn_pt);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.eta", &HGamTruthElectronsAuxDyn_eta, &b_HGamTruthElectronsAuxDyn_eta);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.m", &HGamTruthElectronsAuxDyn_m, &b_HGamTruthElectronsAuxDyn_m);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.px", &HGamTruthElectronsAuxDyn_px, &b_HGamTruthElectronsAuxDyn_px);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.py", &HGamTruthElectronsAuxDyn_py, &b_HGamTruthElectronsAuxDyn_py);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.pz", &HGamTruthElectronsAuxDyn_pz, &b_HGamTruthElectronsAuxDyn_pz);
   fChain->SetBranchAddress("HGamTruthElectronsAuxDyn.e", &HGamTruthElectronsAuxDyn_e, &b_HGamTruthElectronsAuxDyn_e);
   //fChain->SetBranchAddress("HGamTruthEventInfo", &HGamTruthEventInfo, &b_HGamTruthEventInfo);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.m_h1", &HGamTruthEventInfoAuxDyn_m_h1, &b_HGamTruthEventInfoAuxDyn_m_h1);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.yAbs_yy", &HGamTruthEventInfoAuxDyn_yAbs_yy, &b_HGamTruthEventInfoAuxDyn_yAbs_yy);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.m_h2", &HGamTruthEventInfoAuxDyn_m_h2, &b_HGamTruthEventInfoAuxDyn_m_h2);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.pTt_yy", &HGamTruthEventInfoAuxDyn_pTt_yy, &b_HGamTruthEventInfoAuxDyn_pTt_yy);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.m_yy", &HGamTruthEventInfoAuxDyn_m_yy, &b_HGamTruthEventInfoAuxDyn_m_yy);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.pT_yy", &HGamTruthEventInfoAuxDyn_pT_yy, &b_HGamTruthEventInfoAuxDyn_pT_yy);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.cosTS_yy", &HGamTruthEventInfoAuxDyn_cosTS_yy, &b_HGamTruthEventInfoAuxDyn_cosTS_yy);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.Njets", &HGamTruthEventInfoAuxDyn_Njets, &b_HGamTruthEventInfoAuxDyn_Njets);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.m_jj", &HGamTruthEventInfoAuxDyn_m_jj, &b_HGamTruthEventInfoAuxDyn_m_jj);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.Dy_j_j", &HGamTruthEventInfoAuxDyn_Dy_j_j, &b_HGamTruthEventInfoAuxDyn_Dy_j_j);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.Dphi_yy_jj", &HGamTruthEventInfoAuxDyn_Dphi_yy_jj, &b_HGamTruthEventInfoAuxDyn_Dphi_yy_jj);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.isFiducial", &HGamTruthEventInfoAuxDyn_isFiducial, &b_HGamTruthEventInfoAuxDyn_isFiducial);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.isFiducialKinOnly", &HGamTruthEventInfoAuxDyn_isFiducialKinOnly, &b_HGamTruthEventInfoAuxDyn_isFiducialKinOnly);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.TruthNonInt_met", &HGamTruthEventInfoAuxDyn_TruthNonInt_met, &b_HGamTruthEventInfoAuxDyn_TruthNonInt_met);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.TruthInt_sumet", &HGamTruthEventInfoAuxDyn_TruthInt_sumet, &b_HGamTruthEventInfoAuxDyn_TruthInt_sumet);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.pT_h1", &HGamTruthEventInfoAuxDyn_pT_h1, &b_HGamTruthEventInfoAuxDyn_pT_h1);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.pT_h2", &HGamTruthEventInfoAuxDyn_pT_h2, &b_HGamTruthEventInfoAuxDyn_pT_h2);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.y_h1", &HGamTruthEventInfoAuxDyn_y_h1, &b_HGamTruthEventInfoAuxDyn_y_h1);
   fChain->SetBranchAddress("HGamTruthEventInfoAuxDyn.y_h2", &HGamTruthEventInfoAuxDyn_y_h2, &b_HGamTruthEventInfoAuxDyn_y_h2);
   //fChain->SetBranchAddress("HGamTruthHiggsBosons", &HGamTruthHiggsBosons, &b_HGamTruthHiggsBosons);
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.pt", &HGamTruthHiggsBosonsAuxDyn_pt, &b_HGamTruthHiggsBosonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.eta", &HGamTruthHiggsBosonsAuxDyn_eta, &b_HGamTruthHiggsBosonsAuxDyn_eta);
   */
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.m", &HGamTruthHiggsBosonsAuxDyn_m, &b_HGamTruthHiggsBosonsAuxDyn_m);
   /*
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.px", &HGamTruthHiggsBosonsAuxDyn_px, &b_HGamTruthHiggsBosonsAuxDyn_px);
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.py", &HGamTruthHiggsBosonsAuxDyn_py, &b_HGamTruthHiggsBosonsAuxDyn_py);
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.pz", &HGamTruthHiggsBosonsAuxDyn_pz, &b_HGamTruthHiggsBosonsAuxDyn_pz);
   fChain->SetBranchAddress("HGamTruthHiggsBosonsAuxDyn.e", &HGamTruthHiggsBosonsAuxDyn_e, &b_HGamTruthHiggsBosonsAuxDyn_e);
   //fChain->SetBranchAddress("HGamTruthMuons", &HGamTruthMuons, &b_HGamTruthMuons);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.pt", &HGamTruthMuonsAuxDyn_pt, &b_HGamTruthMuonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.eta", &HGamTruthMuonsAuxDyn_eta, &b_HGamTruthMuonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.m", &HGamTruthMuonsAuxDyn_m, &b_HGamTruthMuonsAuxDyn_m);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.px", &HGamTruthMuonsAuxDyn_px, &b_HGamTruthMuonsAuxDyn_px);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.py", &HGamTruthMuonsAuxDyn_py, &b_HGamTruthMuonsAuxDyn_py);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.pz", &HGamTruthMuonsAuxDyn_pz, &b_HGamTruthMuonsAuxDyn_pz);
   fChain->SetBranchAddress("HGamTruthMuonsAuxDyn.e", &HGamTruthMuonsAuxDyn_e, &b_HGamTruthMuonsAuxDyn_e);
   //fChain->SetBranchAddress("HGamTruthPhotons", &HGamTruthPhotons, &b_HGamTruthPhotons);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.isIsolated", &HGamTruthPhotonsAuxDyn_isIsolated, &b_HGamTruthPhotonsAuxDyn_isIsolated);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.etcone20", &HGamTruthPhotonsAuxDyn_etcone20, &b_HGamTruthPhotonsAuxDyn_etcone20);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.etcone40", &HGamTruthPhotonsAuxDyn_etcone40, &b_HGamTruthPhotonsAuxDyn_etcone40);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.truthOrigin", &HGamTruthPhotonsAuxDyn_truthOrigin, &b_HGamTruthPhotonsAuxDyn_truthOrigin);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.truthType", &HGamTruthPhotonsAuxDyn_truthType, &b_HGamTruthPhotonsAuxDyn_truthType);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.pt", &HGamTruthPhotonsAuxDyn_pt, &b_HGamTruthPhotonsAuxDyn_pt);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.px", &HGamTruthPhotonsAuxDyn_px, &b_HGamTruthPhotonsAuxDyn_px);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.eta", &HGamTruthPhotonsAuxDyn_eta, &b_HGamTruthPhotonsAuxDyn_eta);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.py", &HGamTruthPhotonsAuxDyn_py, &b_HGamTruthPhotonsAuxDyn_py);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.m", &HGamTruthPhotonsAuxDyn_m, &b_HGamTruthPhotonsAuxDyn_m);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.pz", &HGamTruthPhotonsAuxDyn_pz, &b_HGamTruthPhotonsAuxDyn_pz);
   fChain->SetBranchAddress("HGamTruthPhotonsAuxDyn.e", &HGamTruthPhotonsAuxDyn_e, &b_HGamTruthPhotonsAuxDyn_e);
   */
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
