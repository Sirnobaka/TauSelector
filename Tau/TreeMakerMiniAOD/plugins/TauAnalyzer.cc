// -*- C++ -*-
//
// Package:    Tau/TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc Tau/TauAnalyzer/plugins/TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrei Oskin
//         Created:  Sat, 14 Sep 2019 01:15:10 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <stdio.h>

// user include files
//#include "GeneratorInterface/TauolaInterface/interface/TauSpinnerCMS.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"

static double dphi(double phi1, double phi2) {
	double result = TMath::Abs(phi1 - phi2);
	if (result < TMath::Pi()) return result;
	return 2 * TMath::Pi() - result;
}

using namespace std;
using namespace edm;
using namespace reco;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

      bool FindGenTau (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool GenTauDM   (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool DecaychannelMatch (vector<MySimpleParticle> &particles, int p1, int p2=0, int p3=0, int p4=0, int p5=0, int p6=0);
      void GenTauProducts (vector<MySimpleParticle> particles);
      bool FindRecoTau (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool AddVertex   (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool CheckMuon    (const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool CheckElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool AddMET       (const edm::Event& iEvent, const edm::EventSetup& iSetup);

      // ----------member data ---------------------------
};


edm::EDGetTokenT<double> TauSpinnerWTToken_;
edm::EDGetTokenT<double> TauSpinnerWTFlipToken_;
edm::EDGetTokenT<double> TauSpinnerWThminusToken_;
edm::EDGetTokenT<double> TauSpinnerWThplusToken_;
edm::EDGetTokenT<bool>   TauSpinnerWTisValidToken_;
edm::EDGetTokenT<std::vector<double>> PolarimetricVectorToken_;

edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
edm::EDGetTokenT<reco::PFTauCollection> TauCollectionToken_;
edm::EDGetTokenT<reco::PFTauDiscriminator> Token_absIso;
edm::EDGetTokenT<reco::TrackCollection> TrackToken_;
edm::EDGetTokenT<reco::VertexCollection> PVToken_;
edm::EDGetTokenT<reco::MuonCollection> MuonCollectionToken_;
edm::EDGetTokenT<reco::GsfElectronCollection> ElectronCollectionToken_;
edm::EDGetTokenT<reco::PFMETCollection> MetCollectionToken_;
edm::EDGetTokenT<reco::PFJetCollection> JetCollectionToken_;
edm::EDGetTokenT<reco::JetTracksAssociationCollection> m_associator;


TTree * tree;

int t_Run;
int t_Event;
int null;
double nullD;
bool monitoring;
double tauPtMin;
double piPtMin;
double tauEtaMax;
double tauDzMax;
double METcut;

double WT;
double WTFlip;
double WThminus;
double WThplus;
bool   WTisValid;
int Decaychannel;

double gentau_pt, gentau_e, gentau_eta, gentau_phi;
double genpiCh_pt, genpiCh_e, genpiCh_eta, genpiCh_phi;
double genpi0_pt, genpi0_e, genpi0_eta, genpi0_phi;
double gentaunu_pt, gentaunu_e, gentaunu_eta, gentaunu_phi;

double tau_pt, tau_e, tau_eta, tau_phi, tau_dm, tau_dz, tau_q, tau_m, tau_absIso;
double piChar_pt, piChar_e, piChar_eta, piChar_phi, piChar_m, piChar_q;
double piZero_pt, piZero_e, piZero_eta, piZero_phi, piZero_m;

double hx, hy, hz, ht;

double picharged_pt;
double picharged_e;
double x_MC;
std::vector<double> PolarimetricVector;

math::XYZPoint pv_position;
int nVtx;
double met;
double met_phi;
int nTauC;
int nPi0;
int nTau;
int tau_found;
double pipiMass;
double upsilon;
double m_t;
double dPhi;
int nTrks;
double tau_ptDif, tau_etaDif, tau_phiDif;
double piCh_ptDif, piCh_etaDif, piCh_phiDif;
double pi0_ptDif, pi0_etaDif, pi0_phiDif;

TH1D * WTHist;
TH1D * WTFlipHist;
TH1D * WThminusHist;
TH1D * WThplusHist;
TH1D * xHist;
TH1D * xHist_Weight;
TH1D * xHist_unWeight;
TH1D * xHist_reWeight;

TH1D * allTauPt;

TH1D * x_tauP_polMinus;
TH1D * x_tauM_polMinus;
TH1D * x_tauP_polPlus;
TH1D * x_tauM_polPlus;

// Poalrinetric vector components
TH1D * hx_Hist;
TH1D * hy_Hist;
TH1D * hz_Hist;
TH1D * ht_Hist;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig) {

  null            = iConfig.getParameter<int>("null");
  nullD           = iConfig.getParameter<double>("nullD");
  monitoring	  = iConfig.getParameter<bool>("monitoring");
  tauPtMin 	  = iConfig.getParameter<double>("tauPtMin");
  piPtMin 	  = iConfig.getParameter<double>("piPtMin");
  tauEtaMax 	  = iConfig.getParameter<double>("tauEtaMax");
  tauDzMax 	  = iConfig.getParameter<double>("tauDzMax");
  METcut          = iConfig.getParameter<double>("METcut");

  TauSpinnerWTToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("WTCollection"));
  TauSpinnerWTFlipToken_    = consumes<double>(iConfig.getParameter<edm::InputTag>("WTFlipCollection"));
  TauSpinnerWThminusToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("WThminusCollection"));
  TauSpinnerWThplusToken_   = consumes<double>(iConfig.getParameter<edm::InputTag>("WThplusCollection"));
  TauSpinnerWTisValidToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("WTisValidCollection"));
  PolarimetricVectorToken_  = consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("PolarimetricVector"));

  GenParticleToken_         = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  TauCollectionToken_       = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("tauCollection"));
  Token_absIso              = consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("absIsoDiscriminator"));
  PVToken_ 		    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  TrackToken_		    = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  MuonCollectionToken_      = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonCollection"));
  ElectronCollectionToken_  = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electronCollection"));
  MetCollectionToken_       = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metCollection"));
  JetCollectionToken_       = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetCollection"));
  m_associator              = consumes<reco::JetTracksAssociationCollection>(iConfig.getParameter<edm::InputTag>("jetTracks"));

   //now do what ever initialization is needed
   //usesResource("TFileService");

}


TauAnalyzer::~TauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  t_Run   = iEvent.id().run();
  t_Event = iEvent.id().event();

  edm::Handle<double> WTHandle;
  iEvent.getByToken(TauSpinnerWTToken_, WTHandle);
  edm::Handle<double> WTFlipHandle;
  iEvent.getByToken(TauSpinnerWTFlipToken_, WTFlipHandle);
  edm::Handle<double> WThminusHandle;
  iEvent.getByToken(TauSpinnerWThminusToken_, WThminusHandle);
  edm::Handle<double> WThplusHandle;
  iEvent.getByToken(TauSpinnerWThplusToken_, WThplusHandle);
  edm::Handle<bool>   WTisValidHandle;
  iEvent.getByToken(TauSpinnerWTisValidToken_, WTisValidHandle);
  edm::Handle<std::vector<double>> PolarimetricVectorHandle;
  iEvent.getByToken(PolarimetricVectorToken_, PolarimetricVectorHandle);

  WT        = *WTHandle;
  WTFlip    = *WTFlipHandle;
  WThminus  = *WThminusHandle;
  WThplus   = *WThplusHandle;
  WTisValid = *WTisValidHandle;
  PolarimetricVector = *PolarimetricVectorHandle;

  hx = PolarimetricVector[0];
  hy = PolarimetricVector[1];
  hz = PolarimetricVector[2];
  ht = PolarimetricVector[3];

  if (!AddVertex(iEvent, iSetup)) {
    std::cout << "Vertex" << std::endl;
    return;
  }

  if (CheckMuon(iEvent, iSetup)) {
    std::cout << "Muon" << std::endl;
    return;
  }

  if (CheckElectron(iEvent, iSetup)) {
    std::cout << "Electron" << std::endl;
    return;
  }

  if (!AddMET(iEvent, iSetup)) {
    std::cout << "MET" << std::endl;
    return;
  }

  if(!GenTauDM(iEvent, iSetup)) {
    std::cout << "Decay mode" << std::endl;
    return;
  }

  if(!FindRecoTau(iEvent, iSetup)) {
    std::cout << "Reconstruction" << std::endl;
    return;
  }

  if((WT >= 0) && (WT <= 2)) {
    WTHist->Fill(WT);
    WTFlipHist->Fill(WTFlip);
  }

  TLorentzVector pi0, pi1;
  pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
  pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

  pipiMass = (pi0 + pi1).M();
//  upsilon  = (pi1.E() - pi0.E()) / tau_pt;
  upsilon  = 2 * piChar_pt / tau_pt - 1;
  //RHT10    = tau_pt / jetPtSum10;
  //RHT15    = tau_pt / jetPtSum15;
  //RHT20    = tau_pt / jetPtSum20;
  dPhi     = dphi(tau_phi, met_phi);
  m_t      = sqrt(2 * tau_pt * met * (1 - cos(dPhi)));

  tau_ptDif = tau_pt - gentau_pt;
  tau_etaDif = tau_eta - gentau_eta;
  tau_phiDif = tau_phi - gentau_phi;

  piCh_ptDif = piChar_pt - genpiCh_pt;
  piCh_phiDif = piChar_phi - genpiCh_phi;
  piCh_etaDif = piChar_eta - genpiCh_eta;

  pi0_ptDif = piZero_pt - genpi0_pt;
  pi0_phiDif = piZero_phi - genpi0_phi;
  pi0_etaDif = piZero_eta - genpi0_eta;

  tree->Fill();

}


bool TauAnalyzer::FindGenTau(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bool tauFound = false;
  edm::Handle<reco::GenParticleCollection> GP;
  iEvent.getByToken(GenParticleToken_, GP);
  
  if(GP.isValid()) {
    std::cout << "GenParticles collection is valid" << std::endl;
    for(unsigned i = 0 ; i < GP->size() ; i++) {
      if (abs((*GP)[i].pdgId())==15) {
	int tau_pdgid = (*GP)[i].pdgId();

        std::cout << "Initial tau daughters parameters:" << std::endl;
	std::vector<MySimpleParticle> tau_daughters;
        for (unsigned k = 0; k < (*GP)[i].numberOfDaughters(); k++) {
          MySimpleParticle tp((*GP)[i].daughter(k)->p4().Px(), (*GP)[i].daughter(k)->p4().Py(), (*GP)[i].daughter(k)->p4().Pz(), (*GP)[i].daughter(k)->p4().E(), (*GP)[i].daughter(k)->pdgId());
          tau_daughters.push_back(tp);
        }

        vector<int>  pdgid;
        double *HH = new double[4];

        for (unsigned int l = 0; l < tau_daughters.size(); l++) {
          pdgid.push_back(tau_daughters[l].pdgid());
        }

        if( Decaychannel == 4) {
	  std::cout << "tau -> pi + tau_nu = " << endl;
	  double Tau_e = (*GP)[i].energy();
	  double PiCh_e = tau_daughters[1].e();
	  picharged_pt = sqrt(TMath::Power(tau_daughters[1].px(), 2) + TMath::Power(tau_daughters[1].py(), 2));
	  picharged_e = tau_daughters[1].e();
	  x_MC = PiCh_e/Tau_e;
	  std::cout << "Polarimetric vector of tau: " << std::endl;
	  for (unsigned l = 0; l < PolarimetricVector.size(); l++) {
	    std::cout << "HH[" << l << "] = " << PolarimetricVector[l] << std::endl;
          }
	  if ((x_MC >= 0) && (x_MC <= 1)) {
              std::cout << "WT = " << WT << std::endl;
              std::cout << "WTFlip = " << WTFlip << std::endl;
              xHist->Fill(x_MC);
              xHist_Weight->Fill(x_MC, WT);
              xHist_unWeight->Fill(x_MC, 1/WT);
              xHist_reWeight->Fill(x_MC, WTFlip);
	      hx_Hist->Fill(hx);
	      hy_Hist->Fill(hy);
	      hz_Hist->Fill(hz);
	      ht_Hist->Fill(ht);
            }
          tauFound = true;
        }	
      }
    }
  }
  return tauFound;
}

bool TauAnalyzer::GenTauDM (const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::GenParticleCollection> GP;
  iEvent.getByToken(GenParticleToken_, GP);

  int channel = null;
  bool ChannelOK = false;

  if(GP.isValid()) {
    if (monitoring) std::cout << "GenParticles collection is valid" << std::endl;
    for(unsigned i = 0 ; i < GP->size() ; i++) {
      if (abs((*GP)[i].pdgId())==15) {
        channel = -1;
	int tau_pdgid = (*GP)[i].pdgId();
        gentau_e   = (*GP)[i].p4().E();
        gentau_pt  = (*GP)[i].pt();
        gentau_eta = (*GP)[i].eta();
        gentau_phi = (*GP)[i].phi();

        //std::vector<SimpleParticle> tau_daughters_simple;
	std::vector<MySimpleParticle> tau_daughters;
        for (unsigned k = 0; k < (*GP)[i].numberOfDaughters(); k++) {
          MySimpleParticle tp((*GP)[i].daughter(k)->p4().Px(), (*GP)[i].daughter(k)->p4().Py(), (*GP)[i].daughter(k)->p4().Pz(), (*GP)[i].daughter(k)->p4().E(), (*GP)[i].daughter(k)->pdgId());
          tau_daughters.push_back(tp);
        }

        std::vector<int>  pdgid;
        for (unsigned int l = 0; l < tau_daughters.size(); l++) {
          pdgid.push_back(tau_daughters[l].pdgid());
        }
        if( pdgid.size() == 2 &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-211) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 211) ) ||
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-321) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 321) )
            )
          ) {
          channel = 3;
          if (abs(pdgid[1]) == 321) channel = 6;
        }
        else if ( pdgid.size() == 3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 111) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 130) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 130) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 310) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 310) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-321, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 321, 111) )
            )
          ) {
 
          channel = 4;
          if(abs(pdgid[1])==321 || abs(pdgid[2])==130 || abs(pdgid[2])==310) channel = 7;
	  GenTauProducts(tau_daughters);
          ChannelOK = true;
	  if (monitoring) {
	    std::cout << "gentau_pt  = " << gentau_pt << std::endl;
	    std::cout << "gentau_eta = " << gentau_eta << std::endl;
	    std::cout << "gentau_phi = " << gentau_phi << std::endl;
          }
        }
        else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111,-211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 211) )
            )
          ) {
          channel = 5;
        }
        else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211) )
            )
          ) {
          channel = 6;
        }
        else if( pdgid.size()==5 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 111) )
            )
          ) {
          channel = 8;
        }
        else if( pdgid.size()==5 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 211) )
            )
          ) {
          channel = 9;
        }
        else if( pdgid.size()==3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12) )
            )
          ) {
          channel = 1;
        }
        else if( pdgid.size()==4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12, 22) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12, 22) )
            )
          ) {
          channel = 1;
        }
        else if( pdgid.size()==3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14) )
            )
          ) {
          channel = 2;
        }
        else if( pdgid.size()==4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14, 22) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14, 22) )
            )
          ) {
          channel = 2;
        }
      }
      if (channel == -1 && abs((*GP)[i].pdgId())==15 && monitoring) {
        std::cout << "Unidentified decay channel of tau was found" << std::endl;
	std::cout << "List of daughters:" << std::endl;
        for (unsigned t = 0; t < (*GP)[i].numberOfDaughters(); t++) {
          std::cout << "TauDaughter[" << t << "] pdgId = " << (*GP)[i].daughter(t)->pdgId() << std::endl;
        }
      }
    }
  }
  Decaychannel = channel;
  return ChannelOK;
}

bool TauAnalyzer::FindRecoTau (const edm::Event& event, const edm::EventSetup& setup) {

  tau_found = 0;
  nTauC     = 0;

  edm::Handle<reco::PFTauCollection> taus;
  event.getByToken(TauCollectionToken_, taus);
  if (!taus.isValid() || taus->size() == 0) return false;

  edm::Handle<reco::PFTauDiscriminator> absIso;
  event.getByToken(Token_absIso, absIso);


  size_t index  = taus->size();
    for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
     auto& tau = (*taus)[i];
     cut(tau.pt() > tauPtMin); // tau transverse momentum
     cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity

     auto& chargedHadrons = tau.signalPFChargedHadrCands();
     cut(chargedHadrons.size() == 1); // single charged hadron (pi+)
     auto& pi0s = tau.signalPiZeroCandidates();
     cut(pi0s.size() > 0); // at least one pi0

     if (monitoring) {
       std::cout << "Pi0 from tau candidate:" << std::endl;
       for (unsigned k = 0; k < pi0s.size(); k++) {
         std::cout << k << " pt = " << pi0s[k].pt() << std::endl;
       }
     }

     auto& pi0 = pi0s.front();
     auto& pi1 = *chargedHadrons.front();

     cut(pi1.pt() > piPtMin); // pi+ transverse momentum

     ++nTauC;
     //allTauPt->Fill(tau.pt());

     cut((pv_position - tau.vertex()).R() < tauDzMax); // tau vertex displacement

     if (index < taus->size()) {
       double iso = absIso->value(i);
       double minIso = absIso->value(index);
       cut( iso < minIso || (iso == minIso && tau.pt() > (*taus)[index].pt()) )
     }

    index = i;
#undef cut
  };

  if (index == taus->size()) return false;

  auto& tau = (*taus)[index];
  auto& pi0 = tau.signalPiZeroCandidates().front();
  auto& pi1 = *tau.signalPFChargedHadrCands().front();
  tau_pt                     = tau.pt();
  tau_eta                    = tau.eta();
  tau_phi                    = tau.phi();
  tau_dm                     = tau.decayMode();
  tau_dz                     = (pv_position - tau.vertex()).R();
  tau_q                      = tau.charge();
  tau_m                      = tau.mass();
  tau_absIso                 = absIso->value(index);

  piZero_pt  = pi0.pt();
  piZero_eta = pi0.eta();
  piZero_phi = pi0.phi();
  piZero_m   = pi0.mass();

  piChar_pt  = pi1.pt();
  piChar_eta = pi1.eta();
  piChar_phi = pi1.phi();
  piChar_q   = pi1.charge();
  piChar_m   = pi1.mass();

  nPi0      = tau.signalPiZeroCandidates().size();
  nTau      = taus->size();
  tau_found = 1;

  if (monitoring) {
    std::cout << "tau_pt  = " << tau_pt << std::endl;
    std::cout << "tau_eta = " << tau_eta << std::endl;
    std::cout << "tau_phi = " << tau_phi << std::endl;
  }

  return true;

}

bool TauAnalyzer::CheckMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::MuonCollection> muons;
	iEvent.getByToken(MuonCollectionToken_, muons);
	if (!muons.isValid()) return false;
	for (auto& muon: *muons) if (muon.pt() > 15) return true;
	return false;
}

bool TauAnalyzer::CheckElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::GsfElectronCollection> electrons;
	iEvent.getByToken(ElectronCollectionToken_, electrons);
	if (!electrons.isValid()) return false;
	for (auto& electron: *electrons) if (electron.pt() > 15) return true;
	return false;
}

bool TauAnalyzer::AddMET(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::PFMETCollection> mets;
	iEvent.getByToken(MetCollectionToken_, mets);
	if (!mets.isValid() || !mets->size()) return false;

	auto& MET = mets->front();
	met = MET.pt();
	met_phi = MET.phi();
	if (met < METcut) return false;
	return true;
}

bool TauAnalyzer::AddVertex(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(PVToken_, vertices);
	if (!vertices.isValid()) return false;

	nVtx = vertices->size();
	if (nVtx == 0) return false;
	pv_position = vertices->front().position();
	return true;
}

bool TauAnalyzer::DecaychannelMatch(vector<MySimpleParticle> &particles, int p1, int p2, int p3, int p4, int p5, int p6) {
  // Copy pdgid-s of all particles
  vector<int> list;
  
  // Fill vector list with PGDId-s of tau daughters
  for (unsigned int i = 0; i < particles.size(); i++) list.push_back(particles[i].pdgid());
  
  // Create array out of pdgid-s
  int p[6] = { p1, p2, p3, p4, p5, p6 };

  // 1) Check if 'particles' contain all pdgid-s on the list 'p'
  for (int i = 0; i < 6; i++) {
    // if the pdgid is zero - finish (the list of searched particles ends)
    if(p[i] == 0) break;
    bool found = false;
    // Loop over tau daughters
    for (unsigned int j = 0; j < list.size(); j++) {
      // if pdgid is found - erese it from the list and search for the next one
      if(list[j] == p[i]) {
        found = true;
        list.erase(list.begin()+j);
        break;
      }
    }
    
    if(!found) return false;
  }
  
  // if there are more particles on the list - there is no match
  if(list.size()!=0) return false;

  
  // 2) Rearrange particles to match the order of pdgid-s listed in array 'p'

  vector<MySimpleParticle> newList;
  
  for(int i = 0; i < 6; i++) {
    // if the pdgid is zero - finish
    if(p[i] == 0) break;
    for(unsigned int j = 0; j < particles.size(); j++) {
      // if pdgid is found - copy it to new list and erese from the old one
      if(particles[j].pdgid() == p[i]) {
        newList.push_back(particles[j]);
        particles.erase(particles.begin()+j);
        break;
      }
    }
  }
  
  particles = newList;

  return true;
}

void TauAnalyzer::GenTauProducts (vector<MySimpleParticle> particles) {

  genpiCh_pt = 0; genpiCh_e = 0; genpiCh_phi = 0; genpiCh_eta = 0;
  genpi0_pt = 0;  genpi0_e  = 0; genpi0_phi  = 0; genpi0_eta = 0;
  gentaunu_pt = 0; gentaunu_e = 0; gentaunu_phi = 0; gentaunu_eta = 0;

  for (unsigned int i = 0; i < particles.size(); i++) {
    if (abs(particles[i].pdgid()) == 211 || abs(particles[i].pdgid()) == 321) {
      genpiCh_pt = particles[i].pt();
      genpiCh_e  = particles[i].e();
      genpiCh_phi  = particles[i].getAnglePhi();
      genpiCh_eta  = particles[i].getAngleEta();
    } else if (abs(particles[i].pdgid()) == 111) {
      genpi0_pt = particles[i].pt();
      genpi0_e  = particles[i].e();
      genpi0_phi  = particles[i].getAnglePhi();
      genpi0_eta  = particles[i].getAngleEta();
    } else if (abs(particles[i].pdgid()) == 16) {
      gentaunu_pt = particles[i].pt();
      gentaunu_e  = particles[i].e();
      gentaunu_phi  = particles[i].getAnglePhi();
      gentaunu_eta  = particles[i].getAngleEta();
    } else continue;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void 
TauAnalyzer::beginJob()
{
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("tree", "tree", 1);

  tree->Branch("Decaychannel", &Decaychannel, "Decaychannel/I");
  tree->Branch("WT",&WT,"WT/D");
  tree->Branch("WTFlip",&WTFlip,"WTFlip/D");
  tree->Branch("WThminus",&WThminus,"WThminus/D");
  tree->Branch("WThplus",&WThplus,"WThplus/D");
  tree->Branch("WTisValid",&WTisValid,"WTisValid/D");
  //tree->Branch("gentau_pt",&gentau_pt,"gentau_pt/D");
  //tree->Branch("gentau_e",&gentau_e,"gentau_e/D");
  //tree->Branch("picharged_pt",&picharged_pt,"picharged_pt/D");
  //tree->Branch("picharged_e",&picharged_e,"picharged_e/D");
  //tree->Branch("x_MC",&x_MC,"x_MC/D");
  tree->Branch("hx",&hx,"hx");
  tree->Branch("hy",&hy,"hy");
  tree->Branch("hz",&hz,"hz");
  tree->Branch("ht",&ht,"ht");

  // Reconstructed parameters
  tree->Branch("tau_pt",&tau_pt,"tau_pt/D");
  tree->Branch("tau_eta",&tau_eta,"tau_eta/D");
  tree->Branch("tau_phi",&tau_phi,"tau_phi/D");
  tree->Branch("tau_q",&tau_q,"tau_q/D");
  tree->Branch("tau_m",&tau_m,"tau_m/D");
  tree->Branch("tau_dm",&tau_dm,"tau_dm/D");
  tree->Branch("tau_dz",&tau_dz,"tau_dz/D");
  tree->Branch("tau_absIso",&tau_absIso,"tau_absIso/D");

  tree->Branch("piChar_pt", &piChar_pt, "piChar_pt/D");
  tree->Branch("piChar_eta", &piChar_eta, "piChar_eta/D");
  tree->Branch("piChar_phi", &piChar_phi, "piChar_phi/D");
  tree->Branch("piChar_q", &piChar_q, "piChar_q/D");
  tree->Branch("piChar_m", &piChar_m, "piChar_m/D");

  tree->Branch("piZero_pt", &piZero_pt, "piZero_pt/D");
  tree->Branch("piZero_eta", &piZero_eta, "piZero_eta/D");
  tree->Branch("piZero_phi", &piZero_phi, "piZero_phi/D");
  tree->Branch("piZero_m", &piZero_m, "piZero_m/D");

  tree->Branch("pipiMass", &pipiMass, "pipiMass/D");
  tree->Branch("upsilon", &upsilon, "upsilon/D");

  tree->Branch("met", &met, "met/D");
  tree->Branch("met_phi", &met_phi, "met_phi/D");
  tree->Branch("m_t", &m_t, "m_t/D");
  tree->Branch("dPhi", &dPhi, "dPhi/D");
  
  tree->Branch("nVtx",&nVtx,"nVtx/I");
  tree->Branch("nTrks",&nTrks,"nTrks/I");
  tree->Branch("nTau",&nTau,"nTau/I");
  tree->Branch("nTauC",&nTauC,"nTauC/I");

  tree->Branch("nPi0",&nPi0,"nPi0/I");

  // Generated parameters
  tree->Branch("gentau_pt", &gentau_pt, "gentau_pt/D");
  tree->Branch("gentau_e", &gentau_e, "gentau_e/D");
  tree->Branch("gentau_phi", &gentau_phi, "gentau_phi/D");
  tree->Branch("gentau_eta", &gentau_eta, "gentau_eta/D");

  tree->Branch("genpiCh_pt", &genpiCh_pt, "genpiCh_pt/D");
  tree->Branch("genpiCh_e", &genpiCh_e, "genpiCh_e/D");
  tree->Branch("genpiCh_phi", &genpiCh_phi, "genpiCh_phi/D");
  tree->Branch("genpiCh_eta", &genpiCh_eta, "genpiCh_eta/D");

  tree->Branch("genpi0_pt", &genpi0_pt, "genpi0_pt/D");
  tree->Branch("genpi0_e", &genpi0_e, "genpi0_e/D");
  tree->Branch("genpi0_phi", &genpi0_phi, "genpi0_phi/D");
  tree->Branch("genpi0_eta", &genpi0_eta, "genpi0_eta/D");

  tree->Branch("gentaunu_pt", &gentaunu_pt, "gentaunu_pt/D");
  tree->Branch("gentaunu_e", &gentaunu_e, "gentaunu_e/D");
  tree->Branch("gentaunu_phi", &gentaunu_phi, "gentaunu_phi/D");
  tree->Branch("gentaunu_eta", &gentaunu_eta, "gentaunu_eta/D");

  tree->Branch("tau_ptDif", &tau_ptDif, "tau_ptDif/D");
  tree->Branch("tau_phiDif", &tau_phiDif, "tau_phiDif/D");
  tree->Branch("tau_etaDif", &tau_etaDif, "tau_etaDif/D");

  tree->Branch("piCh_ptDif", &piCh_ptDif, "piCh_ptDif/D");
  tree->Branch("piCh_phiDif", &piCh_phiDif, "piCh_phiDif/D");
  tree->Branch("piCh_etaDif", &piCh_etaDif, "piCh_etaDif/D");

  tree->Branch("pi0_ptDif", &pi0_ptDif, "pi0_ptDif/D");
  tree->Branch("pi0_phiDif", &pi0_phiDif, "pi0_phiDif/D");
  tree->Branch("pi0_etaDif", &pi0_etaDif, "pi0_etaDif/D");

  WTHist        = FS->make<TH1D>("WTHist","WTHist",100,0,2);
  WTFlipHist    = FS->make<TH1D>("WTFlipHist","WTFlipHist",100,0,2);
  WThminusHist  = FS->make<TH1D>("WThminusHist","WThminusHist",100,0,2);
  WThplusHist   = FS->make<TH1D>("WThplusHist","WThplusHist",100,0,2);
  xHist           = FS->make<TH1D>("xHist","x",100,0,1);
  xHist_unWeight  = FS->make<TH1D>("xHist_unWeight","x/WT",100,0,1);
  xHist_Weight    = FS->make<TH1D>("xHist_Weight","x*WT",100,0,1);
  xHist_reWeight  = FS->make<TH1D>("xHist_reWeight","x*WTFlip",100,0,1);
  hx_Hist         = FS->make<TH1D>("hx_Hist","hx",100,-1,1);
  hy_Hist         = FS->make<TH1D>("hy_Hist","hy",100,-1,1);
  hz_Hist         = FS->make<TH1D>("hz_Hist","hz",100,-1,1);
  ht_Hist         = FS->make<TH1D>("ht_Hist","ht",100,-1,1);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void TauAnalyzer::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  //descriptions.addDefault(desc);
  descriptions.add("tauanalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
