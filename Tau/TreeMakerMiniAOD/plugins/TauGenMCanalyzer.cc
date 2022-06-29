// -*- C++ -*-
//
// Package:    Tau/TauGenMCanalyzer
// Class:      TauGenMCanalyzer
//
/**\class TauGenMCanalyzer TauGenMCanalyzer.cc Tau/TauGenMCanalyzer/plugins/TauGenMCanalyzer.cc

 Description: [one line class summary]

 Implementation:
         [Notes on implementation]
*/

// system include files
#include <memory>
#include <iostream>
#include <iomanip> // std::setw
#include <string>
#include <vector>
#include <chrono>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "PhysicsTools/Heppy/interface/TriggerBitChecker.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/PFTauDecayModeTools.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Pi Zero libs
#include "RecoTauTag/RecoTau/interface/RecoTauPiZeroPlugins.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/CombinatoricGenerator.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "CommonTools/CandUtils/src/AddFourMomenta.cc"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/AssociativeIterator.h"

// For Btag calibration
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

// PileUp reweighting
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// deltaPhi
//#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

// Tau ID scale factors calculation
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "TauPOG/TauIDSFs/src/TauIDSFTool.cc"

// Muons corrections
#include "RoccoR/RoccoR.h"
#include "RoccoR/RoccoR.cc"

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/VectorUtil.h"

#include "Tau/TreeMakerMiniAOD/plugins/ParticleMonitor.h"
#include "Tau/TreeMakerMiniAOD/plugins/BJetCandidate.h"
#include "Tau/TreeMakerMiniAOD/plugins/LeptonCandidate.h"
#include "Tau/TreeMakerMiniAOD/plugins/GenRecoMonitor.h"
#include "Tau/TreeMakerMiniAOD/plugins/PU_distributions.h"
#include "Tau/TreeMakerMiniAOD/plugins/TriggerMatching.h"
#include "Tau/TreeMakerMiniAOD/plugins/TauESCorr.cc"
//#include "Tau/TreeMakerMiniAOD/plugins/PiZeroReconstructor.h"
//#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"

// Uncomment this for debugging
//#define DEBUG

void GetLastDaughter(const reco::Candidate* &particle) {
    if (particle->numberOfDaughters() == 1) {
        particle = particle->daughter(0);
        GetLastDaughter(particle);
    } else return;
}

void GetLastMother(const reco::Candidate* &particle) {
    if (particle->mother()->pdgId() == particle->pdgId()) {
        particle = particle->mother();
        GetLastMother(particle);
    }
    else return;
}

void SortGenTaus(std::vector <reco::GenParticle> &items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].pt() < items[i].pt()) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

int CalculateID (pat::Tau TauParticle, std::vector<std::string> WorkingPoints) {
    int TauIDValue = 0;
    for (std::string WP: WorkingPoints) {
        if (TauParticle.tauID(WP) > 0) {
            TauIDValue++;
        }
    }
    return TauIDValue;
};

int CalculateElectronID (pat::Electron EleParticle, std::vector<std::string> WorkingPoints) {
    int EleIDValue = 0;
    for (std::string WP: WorkingPoints) {
        if (EleParticle.electronID(WP) > 0) {
            EleIDValue++;
        }
    }
    return EleIDValue;
};
//
// class declaration
//

class TauGenMCanalyzer : public edm::EDAnalyzer {
public:
    explicit TauGenMCanalyzer(const edm::ParameterSet&);
    ~TauGenMCanalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    void FindGenTau        (const edm::Event&);
    int  AddTaus           (const edm::Event&);
    int  AddLeptons        (const edm::Event&);
    bool AddMET            (const edm::Event&);
    bool AddJets           (const edm::Event&);
    bool AddVertex         (const edm::Event&);
    void AddTriggers       (const edm::Event&);
    void CountTracks       (const edm::Event&);
    void ClearVectors      ();
    void GetGenWeight      (const edm::Event&);
    void GetPuMCWeight     (const edm::Event&);

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    
    // ----------member data ---------------------------

    edm::EDGetTokenT<pat::TauCollection> TauCollectionToken_;
    edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
    edm::EDGetTokenT<pat::ElectronCollection> ElectronCollectionToken_;
    edm::EDGetTokenT<pat::JetCollection> PuppiJetCollectionToken_;
    edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
    edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
    edm::EDGetTokenT<pat::METCollection> PuppiMetCollectionToken_;
    edm::EDGetTokenT<reco::VertexCollection> PVToken_;
    // SV coolection
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SVToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
    edm::EDGetTokenT<reco::GenJetCollection> tok_GenAK4Jets_;
    // PF candidates collection
    edm::EDGetTokenT<pat::PackedCandidateCollection> PackedCandidateCollectionToken_;
    
    edm::EDGetTokenT<pat::IsolatedTrackCollection> TrackToken_;
    edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoToken_;
    edm::EDGetTokenT<GenRunInfoProduct> GenRunInfoToken_;
    edm::EDGetTokenT<LHEEventProduct> LHEInfoToken_;

    edm::LumiReWeighting LumiWeights_;
    edm::EDGetTokenT<edm::View<pat::Tau> > tauSrcToken_;

    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> tok_PuInfo;
    edm::EDGetTokenT<double> rhoTag;
    EffectiveAreas* effectiveAreas;

    // HLT
    edm::InputTag theTriggerResultsLabel;
    edm::EDGetTokenT<edm::TriggerResults> tok_trigRes;
    std::vector<std::string>  trigNamesTarget1;
    std::vector<std::string>  trigNamesTarget2;
    std::vector<std::string>  trigNamesTarget3;
    std::vector<std::string>  trigNames1;
    std::vector<std::string>  trigNames2;
    std::vector<std::string>  trigNames3;
    std::vector<std::string>  trigNames4;
    std::vector<std::string>  trigNames5;
    std::vector<std::string>  trigNames6;
    std::vector<std::string>  trigNames7;
    std::vector<std::string>  trigNamesSelected;
    // trigger prescales
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> tok_triggerObjects;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescales;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescalesL1min;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescalesL1max;

    //////////////////////////////////////////////////////////////////////////////////////////////
    
    TTree * tree;
    TH1F * allTauPt;
    
    UInt_t t_Run;
    UInt_t t_Event;

    // Generated particles parameters
    std::vector <const reco::GenParticle*> GenTauCandidates;
    std::vector <const reco::GenParticle*> GenLepCandidates;
    // Gen taus
    std::vector <double> genTauPt;
    std::vector <double> genTauEta;
    std::vector <double> genTauPhi;
    std::vector <double> genTauE;
    std::vector <int> genTauCharge;
    std::vector <int> genTauPDGId;
    std::vector <int> genTauStatus;
    std::vector <int> genTauIsLastCopy;
    std::vector <int> genTauDM;
    std::vector <int> genTauDecayK0;
    std::vector <int> genTauMother;
    std::vector <int> genTauTrueMother;
    std::vector <int> genTauFromW;
    std::vector <int> genTauFromWFromt;
    // PiChar and PiZero
    std::vector <double> genPiCharPt;
    std::vector <double> genPiCharE;
    std::vector <double> genPiCharEta;
    std::vector <double> genPiCharPhi;
    std::vector <int>    genPiCharQ;
    std::vector <double> genPi0Pt;
    std::vector <double> genPi0E;
    std::vector <double> genPi0Eta;
    std::vector <double> genPi0Phi;
    std::vector <double> genTauVisPt;
    std::vector <double> genTauVisEta;
    std::vector <double> genTauVisPhi;
    std::vector <double> genTauVisEnergy;

    // Gen leptons
    std::vector <double> genLeptonPt;
    std::vector <double> genLeptonEta;
    std::vector <double> genLeptonPhi;
    std::vector <double> genLeptonE;
    std::vector <int> genLeptonCharge;
    std::vector <int> genLeptonPDGId;
    std::vector <int> genLeptonStatus;
    std::vector <int> genLeptonIsLastCopy;
    std::vector <int> genLeptonMother;
    //std::vector <int> genLeptonTrueMother;
    std::vector <int> genLeptonIsTauDecayProduct;
    std::vector <int> genLeptonFromW;
    std::vector <int> genLeptonFromWFromt;

    // Gen quarks
    std::vector <double> genQuarkPt;
    std::vector <double> genQuarkEta;
    std::vector <double> genQuarkPhi;
    std::vector <double> genQuarkEnergy;
    std::vector <int>    genQuarkMother;
    std::vector <int>    genQuarkPDGId;
    std::vector <int>    genQuarkStatus;
    std::vector <int>    genQuarkFromT;
    //std::vector <int>    genQuarkCharge;
    std::vector <double> genTopPt;
    std::vector <double> genTopEta;
    std::vector <double> genTopPhi;
    std::vector <double> genTopEnergy;
    std::vector <int>    genTopPDGId;
    std::vector <int>    genTopStatus;
    std::vector <int>    genTopWdaughterPDGId;

    // Reco taus
    std::vector <double> TauPt;
    std::vector <double> TauEta;
    std::vector <double> TauPhi;
    std::vector <int>    TauDM;
    std::vector <int>    TauCharge;
    std::vector <double> TauMass;
    std::vector <double> TauDz;

    // New tau DM IDs
    std::vector <int>   TauMVADM2017_v1;
    std::vector <float> TauMVADM2017_v1_DM0raw;
    std::vector <float> TauMVADM2017_v1_DM1raw;
    std::vector <float> TauMVADM2017_v1_DM2raw;
    std::vector <float> TauMVADM2017_v1_DM10raw;
    std::vector <float> TauMVADM2017_v1_DM11raw;
    std::vector <float> TauMVADM2017_v1_DMOtherraw;

    std::vector <int> TaunGamma;

    std::vector <double> TauDeep2017v2p1ElectronRejection;
    std::vector <double> TauDeep2017v2p1MuonRejection;
    std::vector <double> TauDeep2017v2p1JetRejection;
    std::vector <double> TauAbsIso;
    std::vector <int>    TauFindingNewDMs;

    std::vector <int> TauDeepTau2017v2p1VSjetWP;
    std::vector <int> TauDeepTau2017v2p1VSmuWP;
    std::vector <int> TauDeepTau2017v2p1VSeWP;
    // Tau decay products
    std::vector <double> PiCharPt;
    std::vector <double> PiCharEta;
    std::vector <double> PiCharPhi;
    std::vector <int>    PiCharQ;
    std::vector <double> PiCharMass;

    std::vector <double> PiZeroPt;
    std::vector <double> PiZeroEta;
    std::vector <double> PiZeroPhi;
    std::vector <double> PiZeroMass;

    const reco::Vertex* Primary_vertex;
    math::XYZPoint pv_position;

    // Reco electrons
    std::vector <double> ElePt;
    std::vector <double> EleEta;
    std::vector <double> ElePhi;
    std::vector <double> EleEnergy;
    std::vector <int>    EleCharge;
    std::vector <int>    EleFlavor;
    std::vector <double> EleDz;
    std::vector <double> EleDxy;
    std::vector <double> EleClosestCtfDz;
    std::vector <double> EleClosestCtfDxy;
    std::vector <double> EleSuperClusterEta;
    std::vector <double> EleRelIso;
    std::vector <int>    EleCutBasedID;
    std::vector <int>    EleMVAisoID;
    std::vector <int>    EleMVAnoIsoID;
    std::vector <double> EleEcalTrkEnergyPostCorr;
    std::vector <double> EleEcalTrkEnergyErrPostCorr;

    // Reco muons
    std::vector <double> MuPt;
    std::vector <double> MuEta;
    std::vector <double> MuPhi;
    std::vector <double> MuEnergy;
    std::vector <int>    MuCharge;
    std::vector <int>    MuFlavor;
    std::vector <double> MuDz;
    std::vector <double> MuDxy;
    std::vector <double> MuRelIso;
    std::vector <int>    MuCutBasedId;
    std::vector <int>    MuMVAId;
    std::vector <int>    MuPFIso;
    std::vector <int>    MuTkIso;
    std::vector <int>    MuMiniIso;
    std::vector <int>    MuMultiIso;

    // MET
    double Met;
    //double Met_phi;
    double Met_eta;
    //double Met_energy;
    double Met_significance;
    double Met_mEtSig;

    int nVtx;
    int nTrks;

    std::vector <double> JetPt;
    std::vector <double> JetEta;
    std::vector <double> JetPhi;
    std::vector <double> JetMass;
    std::vector <double> JetEnergy;
    std::vector <double> JetBprob;
    std::vector <double> JetBBprob;
    std::vector <double> JetLepbprob;
    std::vector <double> JetBprobCSV;
    std::vector <double> JetBBprobCSV;
    std::vector <int>    JetHadronFlavour;
    std::vector <double> JetPVdR;
    std::vector <int>    JetCharge;

    //////////////////////////////////////////////////////

    int resultTriggerWeight;
    int triggerPrescaleHLT;
    int triggerPrescaleL1max;
    int triggerPrescaleL1min;

    UInt_t TriggerBit1; // 32 bit
    UInt_t TriggerBit2;
    UInt_t TriggerBit3;
    UInt_t TriggerBit4;
    UInt_t TriggerBit5;
    UInt_t TriggerBit6;
    UInt_t TriggerBit7;
    UInt_t TriggerBit_selected;

    double PU_weight;
    int    nPUv;
    float  Tnpv;

    double GenEventInfoWeight;
    double GenRunInfoCrossSection;
    double GenRunInfofilterEfficiency;
    double LHEEventInfoWeight;

    //////////////////////////////////////////////////////
    bool isMC;
    bool noRejection;
    bool noTauSelection;
    bool noLepSelection;
    double BtagThr;
    bool SelectBquarks;
    double tauPtMin;
    double piPtMin;
    double tauEtaMax;
    double tauDzMax;
    double BJetPtMin;
    double JetEtaMax;
    double MuElePtMin;
    double EtaMax;
    double null;
    bool UsePuppiJets;

    std::vector< float > MCPileUp;
    std::vector< float > DataPileUp;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauGenMCanalyzer::TauGenMCanalyzer(const edm::ParameterSet& iConfig) {
    //now do what ever initialization is needed

    isMC                        = iConfig.getParameter<bool>("isMC");
    noRejection                 = iConfig.getParameter<bool>("noRejection");
    noTauSelection              = iConfig.getParameter<bool>("noTauSelection");
    noLepSelection              = iConfig.getParameter<bool>("noLepSelection");
    SelectBquarks               = iConfig.getParameter<bool>("SelectBquarks");
    BtagThr                     = iConfig.getParameter<double>("BtagThr");
    tauPtMin                    = iConfig.getParameter<double>("tauPtMin");
    piPtMin                     = iConfig.getParameter<double>("piPtMin");
    tauEtaMax                   = iConfig.getParameter<double>("tauEtaMax");
    tauDzMax                    = iConfig.getParameter<double>("tauDzMax");
    BJetPtMin                   = iConfig.getParameter<double>("BJetPtMin"); // 30
    JetEtaMax                   = iConfig.getParameter<double>("JetEtaMax");
    MuElePtMin                  = iConfig.getParameter<double>("MuElePtMin"); // 20
    EtaMax                      = iConfig.getParameter<double>("EtaMax"); // 2.4
    null                        = iConfig.getParameter<double>("null");

    std::string tauCollection         = iConfig.getParameter<std::string>("tauCollection");
    std::string muonCollection        = iConfig.getParameter<std::string>("muonCollection");
    std::string electronCollection    = iConfig.getParameter<std::string>("electronCollection");
    std::string PuppijetCollection    = iConfig.getParameter<std::string>("PuppijetCollection");
    std::string jetCollection         = iConfig.getParameter<std::string>("jetCollection");
    std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
    std::string PuppimetCollection    = iConfig.getParameter<std::string>("PuppimetCollection");
    std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
    std::string SVCollection          = iConfig.getParameter<std::string>("SVCollection");
    std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
    std::string GenJetsCollection     = iConfig.getParameter<std::string>("GenJetsCollection");
    std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
    std::string PackedCandidateCollection = iConfig.getParameter<std::string>("PackedCandidateCollection");
    theTriggerResultsLabel            = edm::InputTag("TriggerResults","","HLT");
    trigNamesTarget1                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget1");
    trigNamesTarget2                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget2");
    trigNamesTarget3                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget3");
    trigNames1                        = iConfig.getParameter<std::vector<std::string>>("Triggers1");
    trigNames2                        = iConfig.getParameter<std::vector<std::string>>("Triggers2");
    trigNames3                        = iConfig.getParameter<std::vector<std::string>>("Triggers3");
    trigNames4                        = iConfig.getParameter<std::vector<std::string>>("Triggers4");
    trigNames5                        = iConfig.getParameter<std::vector<std::string>>("Triggers5");
    trigNames6                        = iConfig.getParameter<std::vector<std::string>>("Triggers6");
    trigNames7                        = iConfig.getParameter<std::vector<std::string>>("Triggers7");
    trigNamesSelected                 = iConfig.getParameter<std::vector<std::string>>("SelectedTriggers");
    
    TauCollectionToken_         = consumes<pat::TauCollection>(edm::InputTag(tauCollection));
    tauSrcToken_                = consumes<edm::View<pat::Tau>>(edm::InputTag(tauCollection));
    MuonCollectionToken_        = consumes<pat::MuonCollection>(edm::InputTag(muonCollection));
    ElectronCollectionToken_    = consumes<pat::ElectronCollection>(edm::InputTag(electronCollection));
    PuppiJetCollectionToken_    = consumes<pat::JetCollection>(edm::InputTag(PuppijetCollection));
    JetCollectionToken_         = consumes<pat::JetCollection>(edm::InputTag(jetCollection));
    MetCollectionToken_         = consumes<pat::METCollection>(edm::InputTag(metCollection));
    PuppiMetCollectionToken_    = consumes<pat::METCollection>(edm::InputTag(PuppimetCollection));
    PVToken_                    = consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
    SVToken_                    = consumes<reco::VertexCompositePtrCandidateCollection>(edm::InputTag(SVCollection));
    GenParticleToken_           = consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
    tok_GenAK4Jets_             = consumes<reco::GenJetCollection>(edm::InputTag(GenJetsCollection));
    TrackToken_                 = consumes<pat::IsolatedTrackCollection>(edm::InputTag(trackCollection));
    PackedCandidateCollectionToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag(PackedCandidateCollection));
    tok_trigRes                 = consumes<edm::TriggerResults>(theTriggerResultsLabel);
    tok_triggerPrescales        = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")); 
    tok_triggerPrescalesL1min   = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescalesL1min")); 
    tok_triggerPrescalesL1max   = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescalesL1max"));
    tok_triggerObjects          = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("Triggerobjects"));
    std::string GenEventInfo          = iConfig.getParameter<std::string>("GenEventInfo");
    std::string LHEEventInfo          = iConfig.getParameter<std::string>("LHEEventInfo");
    
    if (isMC) {
        GenEventInfoToken_      = consumes<GenEventInfoProduct>(edm::InputTag(GenEventInfo));
        GenRunInfoToken_        = consumes<GenRunInfoProduct, edm::InRun>(edm::InputTag(GenEventInfo));
        LHEInfoToken_           = consumes<LHEEventProduct>(edm::InputTag(LHEEventInfo));
    }

    tok_PuInfo = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupInfo"));
    rhoTag     = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoTag"));
    effectiveAreas = new EffectiveAreas((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath());
    const std::string PU_MC_File = (iConfig.getParameter<edm::FileInPath>("PU_MC_file")).fullPath();
    const std::string PU_data_File = (iConfig.getParameter<edm::FileInPath>("PU_data_file")).fullPath();

    for( int i=0; i<100; ++i) {
        MCPileUp.push_back(MC_pileip_f1[i]);
        DataPileUp.push_back(Data_pileip_f1[i]);
    }
    //LumiWeights_ = edm::LumiReWeighting("/afs/cern.ch/work/a/aoskin/Tau_Works/crab_TreeMaker/miniAOD_ttbar_UL2017_February2022_fullMC/cfg_files/pileup_data/PU_MC_UL2017.root",
    //"/afs/cern.ch/work/a/aoskin/Tau_Works/crab_TreeMaker/miniAOD_ttbar_UL2017_February2022_fullMC/cfg_files/pileup_data/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root", "pileup", "pileup");
    LumiWeights_ = edm::LumiReWeighting(MCPileUp, DataPileUp);
}


TauGenMCanalyzer::~TauGenMCanalyzer() {
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void TauGenMCanalyzer::analyze(const edm::Event& event, const edm::EventSetup&) {
    t_Run   = event.id().run();
    t_Event = event.id().event();

    bool VertexFound = AddVertex(event);
    if (!VertexFound) {
        ClearVectors();
        return;
    }
    int nTaus = AddTaus(event);
    if (nTaus < 1) {
        ClearVectors();
        return;
    }

    FindGenTau(event);

    int nLeptons = AddLeptons(event);
    if ((nLeptons < 1 && nTaus < 2)) {
        ClearVectors();
        return;
    }
    bool MetFound = AddMET(event);
    if (!MetFound) {
        ClearVectors();
        return;
    }
    CountTracks(event);
    if (!AddJets(event)) {
        ClearVectors();
        return;
    }
    AddTriggers(event);
    GetPuMCWeight(event);
    GetGenWeight(event);

    tree->Fill();
    // Clear vectros
    ClearVectors();
    
};

// ------------ method called once each job just before starting event loop  ------------
void TauGenMCanalyzer::beginJob() {
    // declaring the tree and its branches.
    edm::Service<TFileService> FS;
    tree = FS->make<TTree>("tree", "tree", 1);
    
    tree->Branch("t_Run",  &t_Run,  "t_Run/i");
    tree->Branch("t_Event",&t_Event,"t_Event/i");

    if (isMC) {

        // gen Taus
        tree->Branch("genTauPt", &genTauPt);
        tree->Branch("genTauEta", &genTauEta);
        tree->Branch("genTauPhi", &genTauPhi);
        tree->Branch("genTauE", &genTauE);
        tree->Branch("genTauCharge", &genTauCharge);
        tree->Branch("genTauPDGId", &genTauPDGId);
        tree->Branch("genTauStatus", &genTauStatus);
        tree->Branch("genTauIsLastCopy", &genTauIsLastCopy);
        tree->Branch("genTauDM", &genTauDM);
        tree->Branch("genTauDecayK0", &genTauDecayK0);
        tree->Branch("genTauMother", &genTauMother);
        tree->Branch("genTauFromW", &genTauFromW);
        tree->Branch("genTauFromWFromt", &genTauFromWFromt);
        //
        tree->Branch("genPiCharPt", &genPiCharPt);
        tree->Branch("genPiCharEta", &genPiCharEta);
        tree->Branch("genPiCharPhi", &genPiCharPhi);
        tree->Branch("genPiCharE", &genPiCharE);
        tree->Branch("genPiCharQ", &genPiCharQ);
        tree->Branch("genPi0Pt", &genPi0Pt);
        tree->Branch("genPi0Eta", &genPi0Eta);
        tree->Branch("genPi0Phi", &genPi0Phi);
        tree->Branch("genPi0E", &genPi0E);
        tree->Branch("genTauVisPt", &genTauVisPt);
        tree->Branch("genTauVisEta", &genTauVisEta);
        tree->Branch("genTauVisPhi", &genTauVisPhi);
        tree->Branch("genTauVisEnergy", &genTauVisEnergy);
        // Gen Leptons
        tree->Branch("genLeptonPt", &genLeptonPt);
        tree->Branch("genLeptonEta", &genLeptonEta);
        tree->Branch("genLeptonPhi", &genLeptonPhi);
        tree->Branch("genLeptonE", &genLeptonE);
        tree->Branch("genLeptonCharge", &genLeptonCharge);
        tree->Branch("genLeptonPDGId", &genLeptonPDGId);
        tree->Branch("genLeptonStatus", &genLeptonStatus);
        tree->Branch("genLeptonIsLastCopy", &genLeptonIsLastCopy);
        tree->Branch("genLeptonIsTauDecayProduct", &genLeptonIsTauDecayProduct);
        tree->Branch("genLeptonMother", &genLeptonMother);
        tree->Branch("genLeptonFromW", &genLeptonFromW);
        tree->Branch("genLeptonFromWFromt", &genLeptonFromWFromt);
        // Gen quarks
        tree->Branch("genQuarkPt", &genQuarkPt);
        tree->Branch("genQuarkEta", &genQuarkEta);
        tree->Branch("genQuarkPhi", &genQuarkPhi);
        tree->Branch("genQuarkEnergy", &genQuarkEnergy);
        tree->Branch("genQuarkMother", &genQuarkMother);
        tree->Branch("genQuarkPDGId", &genQuarkPDGId);
        tree->Branch("genQuarkStatus", &genQuarkStatus);
        //tree->Branch("genQuarkCharge", &genQuarkCharge);
        tree->Branch("genQuarkFromT", &genQuarkFromT);
        //
        tree->Branch("genTopPt", &genTopPt);
        tree->Branch("genTopEta", &genTopEta);
        tree->Branch("genTopPhi", &genTopPhi);
        tree->Branch("genTopEnergy", &genTopEnergy);
        tree->Branch("genTopPDGId", &genTopPDGId);
        tree->Branch("genTopStatus", &genTopStatus);
        tree->Branch("genTopWdaughterPDGId", &genTopWdaughterPDGId);
    }

    // Reco Taus
    tree->Branch("TauPt", &TauPt);
    tree->Branch("TauEta", &TauEta);
    tree->Branch("TauPhi", &TauPhi);
    tree->Branch("TauDM", &TauDM);
    tree->Branch("TauCharge", &TauCharge);
    tree->Branch("TauMass", &TauMass);
    tree->Branch("TauDz", &TauDz);
    tree->Branch("TauMVADM2017_v1", &TauMVADM2017_v1);
    tree->Branch("TauMVADM2017_v1_DM0raw", &TauMVADM2017_v1_DM0raw);
    tree->Branch("TauMVADM2017_v1_DM1raw", &TauMVADM2017_v1_DM1raw);
    tree->Branch("TauMVADM2017_v1_DM2raw", &TauMVADM2017_v1_DM2raw);
    tree->Branch("TauMVADM2017_v1_DM10raw", &TauMVADM2017_v1_DM10raw);
    tree->Branch("TauMVADM2017_v1_DM11raw", &TauMVADM2017_v1_DM11raw);
    tree->Branch("TauMVADM2017_v1_DMOtherraw", &TauMVADM2017_v1_DMOtherraw);
    tree->Branch("TaunGamma", &TaunGamma);
    tree->Branch("TauDeep2017v2p1ElectronRejection", &TauDeep2017v2p1ElectronRejection);
    tree->Branch("TauDeep2017v2p1MuonRejection", &TauDeep2017v2p1MuonRejection);
    tree->Branch("TauDeep2017v2p1JetRejection", &TauDeep2017v2p1JetRejection);
    tree->Branch("TauAbsIso", &TauAbsIso);
    tree->Branch("TauFindingNewDMs", &TauFindingNewDMs);
    tree->Branch("TauDeepTau2017v2p1VSjetWP", &TauDeepTau2017v2p1VSjetWP);
    tree->Branch("TauDeepTau2017v2p1VSmuWP", &TauDeepTau2017v2p1VSmuWP);
    tree->Branch("TauDeepTau2017v2p1VSeWP", &TauDeepTau2017v2p1VSeWP);
    tree->Branch("PiCharPt", &PiCharPt);
    tree->Branch("PiCharEta", &PiCharEta);
    tree->Branch("PiCharPhi", &PiCharPhi);
    tree->Branch("PiCharQ", &PiCharQ);
    tree->Branch("PiCharMass", &PiCharMass);
    tree->Branch("PiZeroPt", &PiZeroPt);
    tree->Branch("PiZeroEta", &PiZeroEta);
    tree->Branch("PiZeroPhi", &PiZeroPhi);
    tree->Branch("PiZeroMass", &PiZeroMass);
    // Leptons
    tree->Branch("ElePt", &ElePt);
    tree->Branch("EleEta", &EleEta);
    tree->Branch("ElePhi", &ElePhi);
    tree->Branch("EleEnergy", &EleEnergy);
    tree->Branch("EleCharge", &EleCharge);
    tree->Branch("EleFlavor", &EleFlavor);
    tree->Branch("EleDz", &EleDz);
    tree->Branch("EleDxy", &EleDxy);
    tree->Branch("EleClosestCtfDz", &EleClosestCtfDz);
    tree->Branch("EleClosestCtfDxy", &EleClosestCtfDxy);
    tree->Branch("EleSuperClusterEta", &EleSuperClusterEta);
    tree->Branch("EleRelIso", &EleRelIso);
    tree->Branch("EleCutBasedID", &EleCutBasedID);
    tree->Branch("EleMVAisoID", &EleMVAisoID);
    tree->Branch("EleMVAnoIsoID", &EleMVAnoIsoID);
    tree->Branch("EleEcalTrkEnergyPostCorr", &EleEcalTrkEnergyPostCorr);
    tree->Branch("EleEcalTrkEnergyErrPostCorr", &EleEcalTrkEnergyErrPostCorr);
    //
    tree->Branch("MuPt", &MuPt);
    tree->Branch("MuEta", &MuEta);
    tree->Branch("MuPhi", &MuPhi);
    tree->Branch("MuEnergy", &MuEnergy);
    tree->Branch("MuCharge", &MuCharge);
    tree->Branch("MuFlavor", &MuFlavor);
    tree->Branch("MuDz", &MuDz);
    tree->Branch("MuDxy", &MuDxy);
    tree->Branch("MuRelIso", &MuRelIso);
    tree->Branch("MuCutBasedId", &MuCutBasedId);
    tree->Branch("MuMVAId", &MuMVAId);
    tree->Branch("MuPFIso", &MuPFIso);
    tree->Branch("MuTkIso", &MuTkIso);
    tree->Branch("MuMiniIso", &MuMiniIso);
    tree->Branch("MuMultiIso", &MuMultiIso);
    // MET
    tree->Branch("Met", &Met, "Met/D");
    //tree->Branch("Met_phi", &Met_phi, "Met_phi/D");
    tree->Branch("Met_eta", &Met_eta, "Met_eta/D");
    tree->Branch("Met_significance", &Met_significance, "Met_significance/D");
    tree->Branch("Met_mEtSig", &Met_mEtSig, "Met_mEtSig/D");
    //tree->Branch("Met_energy", &Met_energy, "Met_energy/D");

    tree->Branch("nVtx",&nVtx,"nVtx/I");
    tree->Branch("nTrks",&nTrks,"nTrks/I");

    // Jets
    tree->Branch("JetPt", &JetPt);
    tree->Branch("JetEta", &JetEta);
    tree->Branch("JetPhi", &JetPhi);
    tree->Branch("JetMass", &JetMass);
    tree->Branch("JetEnergy", &JetEnergy);
    tree->Branch("JetBprob", &JetBprob);
    tree->Branch("JetBBprob", &JetBBprob);
    tree->Branch("JetLepbprob", &JetLepbprob);
    tree->Branch("JetBprobCSV", &JetBprobCSV);
    tree->Branch("JetBBprobCSV", &JetBBprobCSV);
    tree->Branch("JetHadronFlavour", &JetHadronFlavour);
    tree->Branch("JetPVdR", &JetPVdR);
    tree->Branch("JetCharge", &JetCharge);

    tree->Branch("TriggerBit1", &TriggerBit1, "TriggerBit1/i");
    tree->Branch("TriggerBit2", &TriggerBit2, "TriggerBit2/i");
    tree->Branch("TriggerBit3", &TriggerBit3, "TriggerBit3/i");
    tree->Branch("TriggerBit4", &TriggerBit4, "TriggerBit4/i");
    tree->Branch("TriggerBit5", &TriggerBit5, "TriggerBit5/i");
    tree->Branch("TriggerBit6", &TriggerBit6, "TriggerBit6/i");
    tree->Branch("TriggerBit7", &TriggerBit7, "TriggerBit7/i");
    tree->Branch("TriggerBit_selected", &TriggerBit_selected, "TriggerBit_selected/i");

    tree->Branch("PU_weight",&PU_weight,"PU_weight/D");
    tree->Branch("Tnpv",&Tnpv,"Tnpv/F");
    tree->Branch("nPUv",&nPUv,"nPUv/I");

    tree->Branch("GenEventInfoWeight",&GenEventInfoWeight,"GenEventInfoWeight/D");
    tree->Branch("GenRunInfoCrossSection",&GenRunInfoCrossSection,"GenRunInfoCrossSection/D");
    tree->Branch("GenRunInfofilterEfficiency",&GenRunInfofilterEfficiency,"GenRunInfofilterEfficiency/D");
    tree->Branch("LHEEventInfoWeight",&LHEEventInfoWeight,"LHEEventInfoWeight/D");

    allTauPt = FS->make<TH1F>("allTauPt","allTauPt",300,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TauGenMCanalyzer::endJob() {

};

void TauGenMCanalyzer::FindGenTau(const edm::Event& event) {

    if (!isMC) return;

    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByToken(GenParticleToken_, genParticles);
    if (!genParticles.isValid()) {
        return;
    }

    const int pdg_tau      = 15;
    const int pdg_pi0      = 111;
    const int pdg_pi1      = 211;
    const int pdg_rho_plus = 213;
    const int pdg_W        = 24;
    const int pdg_nu_tau   = 16;
    const int pdg_electron = 11;
    const int pdg_nu_ele   = 12;
    const int pdg_mu       = 13;
    const int pdg_nu_mu    = 14;
    const int pdg_dquark   = 1;
    const int pdg_uquark   = 2;
    const int pdg_squark   = 3;
    const int pdg_cquark   = 4;
    const int pdg_bquark   = 5;
    const int pdg_tquark   = 6;
    const int pdg_gluon    = 21;
    //
    const int pdg_K0 = 311;
    const int pdg_K0S = 310;
    const int pdg_K0L = 130;
    const int pdg_KChar = 321;

    int nTop = 0;

    #define cut(condition) if (!(condition)) continue;

    for (auto& particle: *genParticles) {
        if (abs(particle.pdgId()) == pdg_tau) {
            cut(particle.pt() > 8.);
            // If tau radiate photon
            if (particle.numberOfDaughters() >= 2 && (abs(particle.daughter(0)->pdgId()) == pdg_tau || abs(particle.daughter(1)->pdgId()) == pdg_tau)) continue;
            //// investigate the source of tau
            const reco::Candidate* pi0 = nullptr;
            const reco::Candidate* pi1 = nullptr;
            math::XYZTLorentzVector TauVisP4;
            //
            int GenTauParticleFromW = 0;
            int GenTauParticleFromWFromt = 0;
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                if (abs(p->pdgId()) == pdg_W) {
                    GenTauParticleFromW = 1;
                    for (auto p1 = p->mother(); p1; p1 = p1->mother()) {
                        if (abs(p1->pdgId()) == pdg_tquark) {
                            GenTauParticleFromWFromt = 1;
                            break;
                        }
                    }
                    break;
                };
            };
            bool genTauMotherFound = false;
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                if (abs(p->pdgId()) == abs((&particle)->pdgId())) {
                    continue;
                } else {
                    genTauTrueMother.push_back(p->pdgId());
                    genTauMotherFound = true;
                    break;
                } 
            }
            if (!genTauMotherFound) {
                genTauTrueMother.push_back(-10);
            }
            genTauFromW.push_back(GenTauParticleFromW);
            genTauFromWFromt.push_back(GenTauParticleFromWFromt);
            genTauPt.push_back(particle.pt());
            genTauEta.push_back(particle.eta());
            genTauPhi.push_back(particle.phi());
            genTauE.push_back(particle.energy());
            genTauCharge.push_back(particle.charge());
            genTauPDGId.push_back(particle.pdgId());
            genTauStatus.push_back(particle.status());
            genTauIsLastCopy.push_back(particle.isLastCopy());
            int genDecayMode = GenTauDecayMode(particle);
            genTauDM.push_back(genDecayMode);
            //genTauMother.push_back(particle.mother()->pdgId());
            int nK0 = 0;
            if (GenTauParticleFromW > 0) {
                GenTauCandidates.push_back(&particle);
            }
            for (unsigned i = 0; i < particle.numberOfDaughters(); ++i) {
                const reco::Candidate* daughter = particle.daughter(i);
                int id = abs(daughter->pdgId());
                if (id != pdg_nu_tau) {
                    TauVisP4 += daughter->p4();
                }
                if (id == pdg_pi0 && !pi0)
                    pi0 = daughter;
                else if (id == pdg_pi1 && !pi1)
                    pi1 = daughter;
                if (id == pdg_K0 || id == pdg_K0S || id == pdg_K0L) {
                    nK0++;
                }
                //else if (id == pdg_nu_tau)
                //    nu_tau = daughter;
            };
            genTauDecayK0.push_back(nK0);
            genTauVisPt.push_back(TauVisP4.pt());
            genTauVisEta.push_back(TauVisP4.eta());
            genTauVisPhi.push_back(TauVisP4.phi());
            genTauVisEnergy.push_back(TauVisP4.energy());
            if (pi1) {
                genPiCharPt.push_back(pi1->pt());
                genPiCharEta.push_back(pi1->eta());
                genPiCharPhi.push_back(pi1->phi());
                genPiCharE.push_back(pi1->energy());
                genPiCharQ.push_back(pi1->charge());
            } else {
                genPiCharPt.push_back(null);
                genPiCharEta.push_back(null);
                genPiCharPhi.push_back(null);
                genPiCharE.push_back(null);
                genPiCharQ.push_back(null);
            }
            if (pi0) {
                genPi0Pt.push_back(pi0->pt());
                genPi0Eta.push_back(pi0->eta());
                genPi0Phi.push_back(pi0->phi());
                genPi0E.push_back(pi0->energy());
            } else {
                genPi0Pt.push_back(null);
                genPi0Eta.push_back(null);
                genPi0Phi.push_back(null);
                genPi0E.push_back(null);
            }
            ////
        } else if (abs(particle.pdgId()) == pdg_mu || abs(particle.pdgId()) == pdg_electron) {
            cut(particle.pt() > 8.);
            cut(particle.statusFlags().isLastCopy());
            //// investigate the source of tau
            int GenLepParticleFromW = 0;
            int GenLepParticleFromWFromt = 0;
            bool genLeptonMotherFound = false;
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                //genTauMother = p->pdgId();
                if (abs(p->pdgId()) == pdg_W) {
                    GenLepParticleFromW = 1;
                    for (auto p1 = p->mother(); p1; p1 = p1->mother()) {
                        if (abs(p1->pdgId()) == pdg_tquark) {
                            GenLepParticleFromWFromt = 1;
                            break;
                        }
                    }
                    break;
                };
            };
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                if (abs(p->pdgId()) == abs((&particle)->pdgId())) {
                    continue;
                } else {
                    //genLeptonTrueMother.push_back(p->pdgId());
                    genLeptonMother.push_back(p->pdgId());
                    genLeptonMotherFound = true;
                    break;
                } 
            }
            if (!genLeptonMotherFound) {
                genLeptonMother.push_back(particle.mother()->pdgId());
            }
            genLeptonFromW.push_back(GenLepParticleFromW);
            genLeptonFromWFromt.push_back(GenLepParticleFromWFromt);
            genLeptonPt.push_back(particle.pt());
            genLeptonEta.push_back(particle.eta());
            genLeptonPhi.push_back(particle.phi());
            genLeptonE.push_back(particle.energy());
            genLeptonCharge.push_back(particle.charge());
            genLeptonPDGId.push_back(particle.pdgId());
            genLeptonStatus.push_back(particle.status());
            genLeptonIsLastCopy.push_back(particle.isLastCopy());
            //genLeptonMother.push_back(particle.mother()->pdgId());
            genLeptonIsTauDecayProduct.push_back(particle.statusFlags().isTauDecayProduct());

            if (GenLepParticleFromW > 0) {
                GenLepCandidates.push_back(&particle);
            }
            ////
        } else if (abs(particle.pdgId()) == pdg_bquark || abs(particle.pdgId()) == pdg_dquark || abs(particle.pdgId()) == pdg_uquark ||
                   abs(particle.pdgId()) == pdg_squark || abs(particle.pdgId()) == pdg_cquark || abs(particle.pdgId()) == pdg_gluon) {
            cut(particle.pt() > 20.);
            cut(particle.statusFlags().isLastCopy());
            if (SelectBquarks) {
                cut(abs(particle.pdgId()) == pdg_bquark);
            }
            //
            genQuarkPt.push_back(particle.pt());
            genQuarkEta.push_back(particle.eta());
            genQuarkPhi.push_back(particle.phi());
            genQuarkEnergy.push_back(particle.energy());
            genQuarkPDGId.push_back(particle.pdgId());
            genQuarkStatus.push_back(particle.status());
            //genQuarkCharge.push_back(particle.charge());
            int GenQuarkParticleFromt = 0;
            bool genQMotherFound = false;
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                if (abs(p->pdgId()) == abs((&particle)->pdgId())) {
                    continue;
                } else {
                    genQuarkMother.push_back(p->pdgId());
                    genQMotherFound = true;
                    break;
                } 
            }
            if (!genQMotherFound) {
                genQuarkMother.push_back(particle.mother()->pdgId());
            }
            for (auto p = (&particle)->mother(); p; p = p->mother()) {
                if (abs(p->pdgId()) == pdg_tquark) {
                    GenQuarkParticleFromt = 1;
                    break;
                };
            };
            genQuarkFromT.push_back(GenQuarkParticleFromt);
            //
        // searching for t-quarks for MC correction
        } else if (abs(particle.pdgId()) == pdg_tquark) {
            cut(particle.statusFlags().isLastCopy());
            //
            nTop++;
            //std::cout << "Top " << nTop << " (" << particle.pdgId() << ") - list of daughters:" << std::endl;
            for (size_t i = 0; i < particle.numberOfDaughters(); ++i) {
                const reco::Candidate* daughter = particle.daughter(i);
                //std::cout << daughter->pdgId() << "(" << i << ") --> ";
                if (abs(daughter->pdgId()) == pdg_W) {
                    for (size_t j = 0; j < daughter->numberOfDaughters(); ++j) {
                        //std::cout << daughter->daughter(j)->pdgId() << "(" << j << "), ";
                        for (auto p = daughter->daughter(0); p; p = p->daughter(0)) {
                            //std::cout << " -> " << p->pdgId();
                            if (abs(p->pdgId()) == pdg_tau || abs(p->pdgId()) == pdg_mu || abs(p->pdgId()) == pdg_electron ||
                                abs(p->pdgId()) == pdg_nu_tau || abs(p->pdgId()) == pdg_nu_mu || abs(p->pdgId()) == pdg_nu_ele) {
                                genTopWdaughterPDGId.push_back(p->pdgId());
                                for (size_t l = 1; l != p->mother()->numberOfDaughters(); ++l) {
                                    //std::cout << " , " << p->mother()->daughter(l)->pdgId();
                                }
                                break;
                            } else if (abs(p->pdgId()) != pdg_W) {
                                genTopWdaughterPDGId.push_back(p->pdgId());
                                for (size_t l = 1; l != p->mother()->numberOfDaughters(); ++l) {
                                    //std::cout << " , " << p->mother()->daughter(l)->pdgId();
                                }
                                break;
                            }
                        }
                    }
                    //std::cout << std::endl;
                } else {
                    //std::cout << std::endl;
                }
            }
            //std::cout << std::endl;
            //
            genTopPt.push_back(particle.pt());
            genTopEta.push_back(particle.eta());
            genTopPhi.push_back(particle.phi());
            genTopEnergy.push_back(particle.energy());
            genTopPDGId.push_back(particle.pdgId());
            genTopStatus.push_back(particle.status());
        }
    }
};

//-----------------------------------------------------------------------------------------

int TauGenMCanalyzer::AddTaus(const edm::Event& event) {

    edm::Handle<pat::TauCollection> taus;
    event.getByToken(TauCollectionToken_, taus);

    std::vector<std::string> JetWPs{"byVVLooseDeepTau2017v2p1VSjet", "byVLooseDeepTau2017v2p1VSjet", "byLooseDeepTau2017v2p1VSjet",
                                    "byMediumDeepTau2017v2p1VSjet", "byTightDeepTau2017v2p1VSjet", "byVTightDeepTau2017v2p1VSjet",
                                    "byVVTightDeepTau2017v2p1VSjet"};
    std::vector<std::string> MuWPs{"byLooseDeepTau2017v2p1VSmu", "byMediumDeepTau2017v2p1VSmu", "byTightDeepTau2017v2p1VSmu"};
    std::vector<std::string> EleWPs{"byVVLooseDeepTau2017v2p1VSe", "byVLooseDeepTau2017v2p1VSe", "byLooseDeepTau2017v2p1VSe",
                                    "byMediumDeepTau2017v2p1VSe", "byTightDeepTau2017v2p1VSe", "byVTightDeepTau2017v2p1VSe", "byVVTightDeepTau2017v2p1VSe"};

    if (taus.isValid() && taus->size() != 0) {
        for (auto& tau: *taus) {
            // Selection criteria
            int TauVsJetDeepID = CalculateID(tau, JetWPs);
            int TauVsMuDeepID = CalculateID(tau, MuWPs);
            int TauVsEleDeepID = CalculateID(tau, EleWPs);
            // Medium
            if (TauVsJetDeepID < 4) continue;
            // Medium
            if (TauVsMuDeepID < 2) continue;
            // Medium
            if (TauVsEleDeepID < 4) continue;
            if (tau.pt() < tauPtMin) continue;
            if (abs(tau.eta()) > tauEtaMax) continue;
            if (!tau.leadChargedHadrCand()) continue;
            if ((pv_position - tau.vertex()).R() > tauDzMax) continue;

            TauPt.push_back(tau.pt());
            TauEta.push_back(tau.eta());
            TauPhi.push_back(tau.phi());
            TauDM.push_back(tau.decayMode());
            TauCharge.push_back(tau.charge());
            TauMass.push_back(tau.mass());
            TauDz.push_back((pv_position - tau.vertex()).R());
            TauMVADM2017_v1.push_back(tau.tauID("MVADM2017v1"));
            TauMVADM2017_v1_DM0raw.push_back(tau.tauID("MVADM2017v1DM0raw"));
            TauMVADM2017_v1_DM1raw.push_back(tau.tauID("MVADM2017v1DM1raw"));
            TauMVADM2017_v1_DM2raw.push_back(tau.tauID("MVADM2017v1DM2raw"));
            TauMVADM2017_v1_DM10raw.push_back(tau.tauID("MVADM2017v1DM10raw"));
            TauMVADM2017_v1_DM11raw.push_back(tau.tauID("MVADM2017v1DM11raw"));
            TauMVADM2017_v1_DMOtherraw.push_back(tau.tauID("MVADM2017v1DMotherraw"));
            TaunGamma.push_back(tau.signalGammaCands().size());

            TauDeep2017v2p1ElectronRejection.push_back(tau.tauID("byDeepTau2017v2p1VSeraw"));
            TauDeep2017v2p1MuonRejection.push_back(tau.tauID("byDeepTau2017v2p1VSmuraw"));
            TauDeep2017v2p1JetRejection.push_back(tau.tauID("byDeepTau2017v2p1VSjetraw"));
            TauAbsIso.push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
            TauFindingNewDMs.push_back(tau.tauID("decayModeFindingNewDMs"));

            TauDeepTau2017v2p1VSjetWP.push_back(TauVsJetDeepID);
            TauDeepTau2017v2p1VSmuWP.push_back(TauVsMuDeepID);
            TauDeepTau2017v2p1VSeWP.push_back(TauVsEleDeepID);

            PiCharPt.push_back(tau.leadChargedHadrCand()->pt());
            PiCharEta.push_back(tau.leadChargedHadrCand()->eta());
            PiCharPhi.push_back(tau.leadChargedHadrCand()->phi());
            PiCharQ.push_back(tau.leadChargedHadrCand()->charge());
            PiCharMass.push_back(tau.leadChargedHadrCand()->mass());

            math::XYZTLorentzVector tau_p4 = tau.p4();
            math::XYZTLorentzVector piChar_p4 = tau.leadChargedHadrCand()->p4();
            math::XYZTLorentzVector piZero_p4 = tau_p4 - piChar_p4;

            PiZeroPt.push_back(piZero_p4.pt());
            PiZeroEta.push_back(piZero_p4.eta());
            PiZeroPhi.push_back(piZero_p4.phi());
            PiZeroMass.push_back(piZero_p4.M());
        }
    }
    return TauPt.size();
};

int TauGenMCanalyzer::AddLeptons(const edm::Event& event) {

    edm::Handle<double> pRho;
    event.getByToken(rhoTag, pRho);

    edm::Handle<pat::ElectronCollection> electrons;
    event.getByToken(ElectronCollectionToken_, electrons);

    edm::Handle<pat::MuonCollection> muons;
    event.getByToken(MuonCollectionToken_, muons);

    std::vector <std::string> EleCutBasedIDWPs{"cutBasedElectronID-Fall17-94X-V2-veto", "cutBasedElectronID-Fall17-94X-V2-loose", "cutBasedElectronID-Fall17-94X-V2-medium", "cutBasedElectronID-Fall17-94X-V2-tight"};
    std::vector <std::string> EleMVAisoWPs{"mvaEleID-Fall17-iso-V2-wpLoose", "mvaEleID-Fall17-iso-V2-wp90", "mvaEleID-Fall17-iso-V2-wp80"};
    std::vector <std::string> EleMVAnoIsoWPs{"mvaEleID-Fall17-noIso-V2-wpLoose", "mvaEleID-Fall17-noIso-V2-wp90", "mvaEleID-Fall17-noIso-V2-wp80"};

    if (electrons.isValid() && !(*electrons).empty()) {
        for (auto& electron : *electrons) {
            if (electron.pt() < MuElePtMin) continue;
            if (abs(electron.eta()) > 2.5) continue;
            reco::GsfTrackRef trackRef;
            reco::TrackRef ClosestCtftrackRef;
            ClosestCtftrackRef = electron.closestCtfTrackRef();
            trackRef = electron.gsfTrack();
            if (abs(electron.eta()) < 1.3 && trackRef.isNonnull()) {
                if (std::abs(trackRef->dxy(pv_position)) > 0.05) continue;
                if (std::abs(trackRef->dz(pv_position)) > 0.1) continue;
            } else if (trackRef.isNonnull()) {
                if (std::abs(trackRef->dxy(pv_position)) > 0.1) continue;
                if (std::abs(trackRef->dz(pv_position)) > 0.2) continue;
            }
            if (electron.electronID("mvaEleID-Fall17-iso-V2-wp90") < 1) continue;

            ElePt.push_back(electron.pt());
            EleEta.push_back(electron.eta());
            ElePhi.push_back(electron.phi());
            EleCharge.push_back(electron.charge());
            EleEnergy.push_back(electron.energy());
            EleFlavor.push_back(electron.pdgId());
            //reco::TrackRef trackRef;
            //trackRef = electron.track();
            if (trackRef.isNonnull()) {
                EleDz.push_back(abs(trackRef->dz(pv_position)));
                EleDxy.push_back(abs(trackRef->dxy(pv_position)));
            } else {
                EleDz.push_back(-1);
                EleDxy.push_back(-1);
            }
            if (ClosestCtftrackRef.isNonnull()) {
                EleClosestCtfDz.push_back(ClosestCtftrackRef->dz(pv_position));
                EleClosestCtfDxy.push_back(ClosestCtftrackRef->dxy(pv_position));
            } else {
                EleClosestCtfDz.push_back(-1);
                EleClosestCtfDxy.push_back(-1);
            }
            EleSuperClusterEta.push_back(electron.superCluster()->eta());
            // Rel Iso
            float rho   = pRho.isValid() ? (*pRho) : 0;
            float chad  = electron.pfIsolationVariables().sumChargedHadronPt;
            float nhad  = electron.pfIsolationVariables().sumNeutralHadronEt;
            float pho   = electron.pfIsolationVariables().sumPhotonEt;
            float eArea = effectiveAreas->getEffectiveArea(fabs(electron.superCluster()->eta()));
            double relIso = (chad + std::max(0.0f, nhad + pho - rho*eArea)) / electron.pt();
            EleRelIso.push_back(relIso);
            // IDs
            int cutBasedID = CalculateElectronID(electron, EleCutBasedIDWPs);
            int MVAisoID   = CalculateElectronID(electron, EleMVAisoWPs);
            int MVAnoIsoID = CalculateElectronID(electron, EleMVAnoIsoWPs);
            EleCutBasedID.push_back(cutBasedID);
            EleMVAisoID.push_back(MVAisoID);
            EleMVAnoIsoID.push_back(MVAnoIsoID);
            EleEcalTrkEnergyPostCorr.push_back(electron.userFloat("ecalTrkEnergyPostCorr"));
            EleEcalTrkEnergyErrPostCorr.push_back(electron.userFloat("ecalTrkEnergyErrPostCorr"));

        }
    }

    // Muons
    if (muons.isValid() && !(*muons).empty()) {
        for (auto& muon : *muons) {
            if (!noLepSelection) {
                if (muon.pt() < MuElePtMin) continue;
                if (abs(muon.eta()) > 2.5) continue;
                if (muon.passed(reco::Muon::MvaLoose) < 1) continue;
                if (!(muon.isPFMuon())) continue;
                if (!(muon.isGlobalMuon())) continue;
                if (muon.numberOfMatchedStations() <= 1) continue;
                reco::TrackRef trackRef = muon.innerTrack();
                if (trackRef.isNonnull()) {
                    if (std::abs(trackRef->dxy(pv_position)) > 0.05) continue;
                    if (std::abs(trackRef->dz(pv_position)) > 0.1) continue;
                    //cut(trackRef->found() > 5);  // more than 5 hits in inner tracker
                }
            }

            MuPt.push_back(muon.pt());
            MuEta.push_back(muon.eta());
            MuPhi.push_back(muon.phi());
            MuCharge.push_back(muon.charge());
            MuEnergy.push_back(muon.energy());
            MuFlavor.push_back(muon.pdgId());
            reco::TrackRef trackRef = muon.innerTrack();
            if (trackRef.isNonnull()) {
                MuDz.push_back(abs(trackRef->dz(pv_position)));
                MuDxy.push_back(abs(trackRef->dxy(pv_position)));
            } else {
                MuDz.push_back(-1);
                MuDxy.push_back(-1);
            }
            reco::MuonPFIsolation iso = muon.pfIsolationR04();
            double relIso = (iso.sumChargedHadronPt + std::max(0.,iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5*iso.sumPUPt)) / muon.pt();
            MuRelIso.push_back(relIso);
            //
            int muonCutBasedId = muon.passed(reco::Muon::CutBasedIdLoose) + muon.passed(reco::Muon::CutBasedIdMedium) + muon.passed(reco::Muon::CutBasedIdTight);
            int muonMVAId = muon.passed(reco::Muon::MvaLoose) + muon.passed(reco::Muon::MvaMedium) + muon.passed(reco::Muon::MvaTight);
            int muonPFIso = muon.passed(reco::Muon::PFIsoVeryLoose) + muon.passed(reco::Muon::PFIsoLoose) + muon.passed(reco::Muon::PFIsoMedium) + muon.passed(reco::Muon::PFIsoTight) + muon.passed(reco::Muon::PFIsoVeryTight) + muon.passed(reco::Muon::PFIsoVeryVeryTight);
            int muonTkIso = muon.passed(reco::Muon::TkIsoLoose) + muon.passed(reco::Muon::TkIsoTight);
            int muonMiniIso = muon.passed(reco::Muon::MiniIsoLoose) + muon.passed(reco::Muon::MiniIsoMedium) + muon.passed(reco::Muon::MiniIsoTight) + muon.passed(reco::Muon::MiniIsoVeryTight);
            int muonMultiIso = muon.passed(reco::Muon::MultiIsoLoose) + muon.passed(reco::Muon::MultiIsoMedium);
            MuCutBasedId.push_back(muonCutBasedId);
            MuMVAId.push_back(muonMVAId);
            MuPFIso.push_back(muonPFIso);
            MuTkIso.push_back(muonTkIso);
            MuMiniIso.push_back(muonMiniIso);
            MuMultiIso.push_back(muonMultiIso);
        }
    }
    return (ElePt.size() + MuPt.size());

};

bool TauGenMCanalyzer::AddVertex(const edm::Event& event) {
    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(PVToken_, vertices);
    //std::cout << "Adding vertex" << std::endl;
    if (!vertices.isValid()) return false;
    //std::cout << "Vertices are valid" << std::endl;

    nVtx = vertices->size();
    //std::cout << nVtx << " vetices" << std::endl;
    if (nVtx == 0) return false;
    Primary_vertex = &vertices->front();
    pv_position = vertices->front().position();
    return true;
};

bool TauGenMCanalyzer::AddJets(const edm::Event& event) {

    edm::Handle<pat::JetCollection> Jets;
    if (UsePuppiJets) {
        event.getByToken(PuppiJetCollectionToken_, Jets);
    } else {
        event.getByToken(JetCollectionToken_, Jets);
    }
    if (!Jets.isValid()) return false;
    for (auto& jet: *Jets) {
        if (jet.pt() < 30) continue;
        if (abs(jet.eta()) > JetEtaMax) continue;
        //double bTagValue = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") + jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
        //if (bTagValue < BtagThr) continue;
        if ((jet.neutralHadronEnergy() / jet.energy() < 0.90) &&
        (jet.neutralEmEnergyFraction() < 0.90) &&
        (jet.muonEnergyFraction() < 0.80) &&
        (jet.chargedHadronEnergyFraction() > 0.) &&
        (jet.chargedMultiplicity() + jet.neutralMultiplicity() > 0.) &&
        (jet.chargedEmEnergyFraction() < 0.80) &&
        (jet.nConstituents() > 1)) {
            //
            JetPt.push_back(jet.pt());
            JetEta.push_back(jet.eta());
            JetPhi.push_back(jet.phi());
            JetMass.push_back(jet.mass());
            JetEnergy.push_back(jet.energy());
            JetBprob.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probb"));
            JetBBprob.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:probbb"));
            JetLepbprob.push_back(jet.bDiscriminator("pfDeepFlavourJetTags:problepb"));
            JetBprobCSV.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probb"));
            JetBBprobCSV.push_back(jet.bDiscriminator("pfDeepCSVJetTags:probbb"));
            JetHadronFlavour.push_back(jet.hadronFlavour());
            JetPVdR.push_back((pv_position - jet.vertex()).R());
            JetCharge.push_back(jet.charge());
        }
    }

    return true;

};

void TauGenMCanalyzer::AddTriggers (const edm::Event& iEvent) {

    resultTriggerWeight = null;
    triggerPrescaleHLT = null;
    triggerPrescaleL1max = null;
    triggerPrescaleL1min = null;

    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(tok_triggerObjects, triggerObjects);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(tok_triggerPrescales, triggerPrescales);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;
    iEvent.getByToken(tok_triggerPrescalesL1min, triggerPrescalesL1min);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
    iEvent.getByToken(tok_triggerPrescalesL1max, triggerPrescalesL1max);
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(tok_trigRes, triggerResults);

    std::vector<std::string> trigNameVec;
    std::vector<bool> trigPassVec;
    std::vector<int> trigPsVec;
    std::vector<int> trigL1minPsVec;
    std::vector<int> trigL1maxPsVec;
    std::vector<int> trigPrescaleVec;

    TriggerBit1 = 0;
    TriggerBit2 = 0;
    TriggerBit3 = 0;
    TriggerBit4 = 0;
    TriggerBit5 = 0;
    TriggerBit6 = 0;
    TriggerBit7 = 0;
    TriggerBit_selected = 0;

    /////////////////////////////TriggerResults////////////////////////////////////
    if (triggerResults.isValid()) {
        const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
        const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
        for ( unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++ ) {
            int hlt    = triggerResults->accept(iHLT);
            //
            const std::string& trigName = triggerNames.triggerName(iHLT);
            int ps        = triggerPrescales->getPrescaleForIndex(iHLT);
            int psL1min   = triggerPrescalesL1min->getPrescaleForIndex(iHLT);
            int psL1max   = triggerPrescalesL1max->getPrescaleForIndex(iHLT);
            bool pass     = triggerResults->accept(iHLT);
            if ( hlt > 0 ) {
                // Write prescale weights to corresponding vectors
                trigPsVec.push_back(ps);
                trigL1minPsVec.push_back(psL1min);
                trigL1maxPsVec.push_back(psL1max);
                trigNameVec.push_back(trigName);
                trigPassVec.push_back(pass);
                trigPrescaleVec.push_back(ps * psL1min);
                for ( unsigned int i=0; i<trigNames1.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames1[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit1 = TriggerBit1 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames2.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames2[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit2 = TriggerBit2 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames3.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames3[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit3 = TriggerBit3 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames4.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames4[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit4 = TriggerBit4 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames5.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames5[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit5 = TriggerBit5 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames6.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames6[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit6 = TriggerBit6 + TMath::Power(2, i_signed);
                    }
                }
                for ( unsigned int i=0; i<trigNames7.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames7[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit7 = TriggerBit7 + TMath::Power(2, i_signed);
                    }
                }
                // Count selected triggers only
                for ( unsigned int i=0; i<trigNamesSelected.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNamesSelected[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit_selected = TriggerBit_selected + TMath::Power(2, i_signed);
                    }
                }
            }
        }
    }

    int requiredTrigPrescale = 9999999;
    int requiredTrigIndex = -1;
    std::string requiredTrigName;
    // Loop over prescale vectors
    for(unsigned int l = 0; l < trigPrescaleVec.size(); l++) {
        if (!trigPassVec[l]) continue;
        if (trigPrescaleVec[l] < requiredTrigPrescale) {
            requiredTrigPrescale = trigPrescaleVec[l];
            requiredTrigIndex = l;
            requiredTrigName = trigNameVec[l];
        }
    }

    if (requiredTrigIndex > -1) {
        resultTriggerWeight = requiredTrigPrescale;
        triggerPrescaleHLT = trigPsVec[requiredTrigIndex];
        triggerPrescaleL1max = trigL1maxPsVec[requiredTrigIndex];
        triggerPrescaleL1min = trigL1minPsVec[requiredTrigIndex];
    }

};

void TauGenMCanalyzer::GetPuMCWeight (const edm::Event& iEvent) {

    Tnpv = -1;
    nPUv = -1;
    PU_weight = 1;

    if (!isMC) return;

    edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
    iEvent.getByToken(tok_PuInfo, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
        int BX = PVI->getBunchCrossing();
        if(BX == 0) { 
            Tnpv = PVI->getTrueNumInteractions();
            nPUv = PVI->getPU_NumInteractions();
            continue;
        }
    }
    PU_weight = LumiWeights_.weight( Tnpv );
};

void TauGenMCanalyzer::GetGenWeight (const edm::Event& event) {

    GenEventInfoWeight = 0;
    LHEEventInfoWeight = 1;

    if (!isMC) return;

    edm::Handle<GenEventInfoProduct> genEvtInfo; 
    event.getByToken(GenEventInfoToken_, genEvtInfo);

    edm::Handle<LHEEventProduct> LHEEvtInfo; 
    event.getByToken(LHEInfoToken_, LHEEvtInfo);

    if (genEvtInfo.isValid()) GenEventInfoWeight = genEvtInfo->weight();
    if (LHEEvtInfo.isValid()) LHEEventInfoWeight = LHEEvtInfo->originalXWGTUP();
};

void TauGenMCanalyzer::CountTracks(const edm::Event& event) {
    nTrks = 0;
    edm::Handle<pat::IsolatedTrackCollection> tracks;
    event.getByToken(TrackToken_, tracks);
    if (!tracks.isValid()) return;
    nTrks = tracks->size();
};

void TauGenMCanalyzer::ClearVectors () {
    GenTauCandidates.clear();
    GenLepCandidates.clear();
    genTauFromW.clear();
    genTauFromWFromt.clear();
    genLeptonFromW.clear();
    genLeptonFromWFromt.clear();
    //
    genTauPt.clear();
    genTauEta.clear();
    genTauPhi.clear();
    genTauE.clear();
    genTauCharge.clear();
    genTauPDGId.clear();
    genTauStatus.clear();
    genTauIsLastCopy.clear();
    genTauDM.clear();
    genTauDecayK0.clear();
    genTauMother.clear();
    genTauTrueMother.clear();
    //
    genPiCharPt.clear();
    genPiCharE.clear();
    genPiCharEta.clear();
    genPiCharPhi.clear();
    genPiCharQ.clear();
    genPi0Pt.clear();
    genPi0E.clear();
    genPi0Eta.clear();
    genPi0Phi.clear();
    genTauVisPt.clear();
    genTauVisEta.clear();
    genTauVisPhi.clear();
    genTauVisEnergy.clear();
    //
    genLeptonPt.clear();
    genLeptonEta.clear();
    genLeptonPhi.clear();
    genLeptonE.clear();
    genLeptonCharge.clear();
    genLeptonPDGId.clear();
    genLeptonStatus.clear();
    genLeptonIsLastCopy.clear();
    genLeptonMother.clear();
    //genLeptonTrueMother.clear();
    genLeptonIsTauDecayProduct.clear();
    //
    genQuarkPt.clear();
    genQuarkEta.clear();
    genQuarkPhi.clear();
    genQuarkEnergy.clear();
    genQuarkMother.clear();
    genQuarkPDGId.clear();
    genQuarkStatus.clear();
    //genQuarkCharge.clear();
    genQuarkFromT.clear();

    genTopPt.clear();
    genTopEta.clear();
    genTopPhi.clear();
    genTopEnergy.clear();
    genTopPDGId.clear();
    genTopStatus.clear();
    genTopWdaughterPDGId.clear();

    TauPt.clear();
    TauEta.clear();
    TauPhi.clear();
    TauDM.clear();
    TauCharge.clear();
    TauMass.clear();
    TauDz.clear();

    TauMVADM2017_v1.clear();
    TauMVADM2017_v1_DM0raw.clear();
    TauMVADM2017_v1_DM1raw.clear();
    TauMVADM2017_v1_DM2raw.clear();
    TauMVADM2017_v1_DM10raw.clear();
    TauMVADM2017_v1_DM11raw.clear();
    TauMVADM2017_v1_DMOtherraw.clear();
    TaunGamma.clear();

    TauDeep2017v2p1ElectronRejection.clear();
    TauDeep2017v2p1MuonRejection.clear();
    TauDeep2017v2p1JetRejection.clear();
    TauAbsIso.clear();
    TauFindingNewDMs.clear();

    TauDeepTau2017v2p1VSjetWP.clear();
    TauDeepTau2017v2p1VSmuWP.clear();
    TauDeepTau2017v2p1VSeWP.clear();

    PiCharPt.clear();
    PiCharEta.clear();
    PiCharPhi.clear();
    PiCharQ.clear();
    PiCharMass.clear();

    PiZeroPt.clear();
    PiZeroEta.clear();
    PiZeroPhi.clear();
    PiZeroMass.clear();

    ElePt.clear();
    EleEta.clear();
    ElePhi.clear();
    EleEnergy.clear();
    EleCharge.clear();
    EleFlavor.clear();
    EleDz.clear();
    EleDxy.clear();
    EleClosestCtfDz.clear();
    EleClosestCtfDxy.clear();
    EleSuperClusterEta.clear();
    EleRelIso.clear();
    EleCutBasedID.clear();
    EleMVAisoID.clear();
    EleMVAnoIsoID.clear();
    EleEcalTrkEnergyPostCorr.clear();
    EleEcalTrkEnergyErrPostCorr.clear();

    // Reco muons
    MuPt.clear();
    MuEta.clear();
    MuPhi.clear();
    MuEnergy.clear();
    MuCharge.clear();
    MuFlavor.clear();
    MuDz.clear();
    MuDxy.clear();
    MuRelIso.clear();
    MuCutBasedId.clear();
    MuMVAId.clear();
    MuPFIso.clear();
    MuTkIso.clear();
    MuMiniIso.clear();
    MuMultiIso.clear();
    // Jets
    JetPt.clear();
    JetEta.clear();
    JetPhi.clear();
    JetMass.clear();
    JetEnergy.clear();
    JetBprob.clear();
    JetBBprob.clear();
    JetLepbprob.clear();
    JetBprobCSV.clear();
    JetBBprobCSV.clear();
    JetHadronFlavour.clear();
    JetPVdR.clear();
    JetCharge.clear();
};

bool TauGenMCanalyzer::AddMET(const edm::Event& event) {

    edm::Handle<pat::METCollection> mets;
    event.getByToken(MetCollectionToken_, mets);

    if (!mets.isValid() || !mets->size()) return false;

    auto& MET = mets->front();
    Met              = MET.pt();
    //Met_phi          = MET.phi();
    Met_eta          = MET.eta();
    //Met_energy       = MET.energy();
    Met_significance = MET.significance();
    Met_mEtSig       = MET.mEtSig();

    return true;
}

// ------------ method called when starting to processes a run  ------------

void TauGenMCanalyzer::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
    GenRunInfoCrossSection = 0;
    GenRunInfofilterEfficiency = 0;

    if (!isMC) return;
        
    edm::Handle<GenRunInfoProduct> genRunInfoProduct;
    Run.getByToken(GenRunInfoToken_, genRunInfoProduct);
    GenRunInfoCrossSection = genRunInfoProduct->crossSection();
    GenRunInfofilterEfficiency = genRunInfoProduct->filterEfficiency();
};

// ------------ method called when ending the processing of a run  ------------

void TauGenMCanalyzer::endRun(edm::Run const& Run, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauGenMCanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.add("TauGenMCanalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauGenMCanalyzer);