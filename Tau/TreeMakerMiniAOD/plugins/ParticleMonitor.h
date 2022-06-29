#include <memory>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include <cstdio>
#include <iostream>
#include <math.h>


// Function for monitoring tau particle parameters
//
void TauMonitor (const pat::Tau &tau, bool DeepTau, const math::XYZPoint &pv_position) {

    reco::CandidatePtrVector VectorPi0 = tau.signalNeutrHadrCands();
    reco::CandidatePtrVector VectorPiCh = tau.signalChargedHadrCands();
    reco::CandidatePtrVector VectorGamma = tau.signalGammaCands();
    reco::CandidatePtrVector VectorSignalCands = tau.signalCands();

    std::cout << "Pt                        = " << tau.pt() << std::endl
    << "eta                       = " << tau.eta() << std::endl
    << "phi                       = " << tau.phi() << std::endl
    << "mass                      = " << tau.mass() << std::endl
    << "Decay mode                    = " << tau.decayMode() << std::endl
    << "decayModeFindingNewDMs        = " << tau.tauID("decayModeFindingNewDMs") << std::endl
    << "decayModeFinding              = " << tau.tauID("decayModeFinding") << std::endl
    << "Size of SignalCands vector    = " << VectorSignalCands.size() << std::endl
    << "Size of PiCh vector           = " << VectorPiCh.size() << std::endl
    << "Size of Pi0 vector            = " << VectorPi0.size() << std::endl
    << "Size of signalNeutrHadrCands  = " << tau.signalNeutrHadrCands().size() << std::endl
    << "Size of signalGammaCands      = " << tau.signalGammaCands().size() << std::endl
    //byCombinedIsolationDeltaBetaCorrRaw3Hits
    << "byCombinedIsolationDeltaBetaCorrRaw3Hits = " << tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << std::endl
    << "chargedIsoPtSum                          = " << tau.tauID("chargedIsoPtSum") << std::endl;
    if (DeepTau) {
        std::cout << "byDeepTau2017v2p1VSeraw                    = " << tau.tauID("byDeepTau2017v2p1VSeraw") << std::endl
        << "byDeepTau2017v2p1VSmuraw                   = " << tau.tauID("byDeepTau2017v2p1VSmuraw") << std::endl
        << "byDeepTau2017v2p1VSjetraw                  = " << tau.tauID("byDeepTau2017v2p1VSjetraw") << std::endl
        << "byLooseDeepTau2017v2p1VSe                  = " << tau.tauID("byLooseDeepTau2017v2p1VSe") << std::endl
        << "byLooseDeepTau2017v2p1VSmu                 = " << tau.tauID("byLooseDeepTau2017v2p1VSmu") << std::endl
        << "byLooseDeepTau2017v2p1VSjet                = " << tau.tauID("byLooseDeepTau2017v2p1VSjet") << std::endl
        // temp
        << "byVLooseDeepTau2017v2p1VSmu      = " << tau.tauID("byVLooseDeepTau2017v2p1VSmu") << std::endl
        << "byVVVLooseDeepTau2017v2p1VSe     = " << tau.tauID("byVVVLooseDeepTau2017v2p1VSe") << std::endl
        << "byVVVLooseDeepTau2017v2p1VSjet   = " << tau.tauID("byVVVLooseDeepTau2017v2p1VSjet") << std::endl;
    } else {
        std::cout << "byIsolationMVArun2v1DBnewDMwLTraw        = " << tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw") << std::endl
        << "byIsolationMVArun2v1PWnewDMwLTraw        = " << tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw") << std::endl
        << "byIsolationMVArun2v1DBdR03oldDMwLTraw    = " << tau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw") << std::endl
        << "byIsolationMVArun2v1DBoldDMwLTraw        = " << tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw") << std::endl
        << "byIsolationMVArun2v1PWdR03oldDMwLTraw    = " << tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw") << std::endl
        << "byIsolationMVArun2v1PWoldDMwLTraw        = " << tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw") << std::endl   
        << "againstElectronMVA6Raw                   = " << tau.tauID("againstElectronMVA6Raw") << std::endl
        // temp
        << "byVLooseIsolationMVArun2v1DBoldDMwLT = " << tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") << std::endl
        << "againstElectronVLooseMVA6             = " << tau.tauID("againstElectronVLooseMVA6") << std::endl
        << "againstMuonLoose3                     = " << tau.tauID("againstMuonLoose3") << std::endl;
    }
    //std::cout << "SV pos   = (" << tau.secondaryVertexPos().x() << ", " << tau.secondaryVertexPos().y() << ", " << tau.secondaryVertexPos().z() << ")" << std::endl;
    //std::cout << "dz(SV, PV) = " << (tau.secondaryVertexPos() - pv_position).R() << std::endl;
    if (tau.hasSecondaryVertex()) {
        std::cout << "has Secondary vertex" << std::endl;
    } else {
        std::cout << "doesn't have Secondary vertex" << std::endl;
    }
    /*
    if (tau.secondaryVertex()->isValid()) {
        std::cout << "Tau secondaryVertex is valid" << std::endl;
        //std::cout << "point = (" << tau.secondaryVertex()->position().x()
        //<< ", " << tau.secondaryVertex()->position().y()
        //<< ", " << tau.secondaryVertex()->position().z() << std:: endl;
    } else {
        std::cout << "Tau secondaryVertex is not valid" << std::endl;
    }
    */
    // Hadrons from tau decay
    /*
    const reco::CandidatePtr leadPiCh = tau.leadChargedHadrCand();
    std::cout << "Lead Pi Charged parameters" << std::ednl <<
    << "pt   = " << leadPiCh->pt() << std::endl
    << "eta  = " << leadPiCh->eta() << std::endl
    << "phi  = " << leadPiCh->phi() << std::endl
    << "mass = " << leadPiCh->mass() << std::endl;
    */

    std::cout << "All candidates:" << std::endl;
    if (VectorSignalCands.size() > 0) {
        for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
            std::cout << "Signal Candidate number " << l << " (pdgId " << VectorSignalCands[l]->pdgId() << ")" << std::endl
            << "Pt    = " << VectorSignalCands[l]->pt() << std::endl
            << "eta   = " << VectorSignalCands[l]->eta() << std::endl
            << "phi   = " << VectorSignalCands[l]->phi() << std::endl
            << "mass  = " << VectorSignalCands[l]->mass() << std::endl;
        }
    }
    std::cout << "Photon candidates:" << std::endl;
    if (tau.signalGammaCands().size() > 0) {
        for (unsigned l = 0; l < tau.signalGammaCands().size(); l++) {
            std::cout << "Gamma Candidate number " << l << " (pdgId " << ((tau.signalGammaCands())[l])->pdgId() << ")" << std::endl
            << "Pt    = " << ((tau.signalGammaCands())[l])->pt() << std::endl
            << "eta   = " << ((tau.signalGammaCands())[l])->eta() << std::endl
            << "phi   = " << ((tau.signalGammaCands())[l])->phi() << std::endl
            << "mass  = " << ((tau.signalGammaCands())[l])->mass() << std::endl
            << "N mothers = " << ((tau.signalGammaCands())[l])->numberOfMothers() << std::endl;
            if (((tau.signalGammaCands())[l])->mother() != nullptr) {
                std::cout << "mother pdgID = " << ((tau.signalGammaCands())[l])->mother()->pdgId() << std::endl;
            }
        }
    }
    if (VectorPiCh.size() > 0) {
        std::cout << "Signal charged candidates:" << std::endl;
        for (unsigned n = 0; n < VectorPiCh.size(); n++) {
            std::cout << "Charged Candidate " << n << ": pdgID = " << VectorPiCh[n]->pdgId() << ", Pt = " << VectorPiCh[n]->pt() << std::endl;
        }
    }
    if (VectorPi0.size() > 0) {
        std::cout << "Signal neutral candidates:" << std::endl;
        for (unsigned n = 0; n < VectorPi0.size(); n++) {
            std::cout << "Charged Candidate " << n << ": pdgID = " << VectorPi0[n]->pdgId() << ", Pt = " << VectorPi0[n]->pt() << std::endl;
        }
    } else {
        std::cout << "No neutral candidates" << std::endl;
    }
}

/*
void GenTauMonitor (const reco::GenParticle &gentau) {

    std::cout << "gentau pt    = " << gentau.pt() << std::endl
    << "gentau eta   = " << gentau.eta() << std::endl
    << "gentau phi   = " << gentau.phi() << std::endl;
}
*/

void GenParticleMonitor (const reco::GenParticle &genparticle, const math::XYZPoint &pv_position) {

    std::cout << "GenParticle pdgId = " << genparticle.pdgId() << std::endl
    << "pt       = " << genparticle.pt() << std::endl
    << "eta      = " << genparticle.eta() << std::endl
    << "phi      = " << genparticle.phi() << std::endl
    << "dZ       = " << (pv_position - genparticle.vertex()).R() << std::endl;
};

void METMonitor (const pat::MET &MET) {

    std::cout << "pt           = " << MET.pt() << std::endl
    << "eta          = " << MET.eta() << std::endl
    << "phi          = " << MET.phi() << std::endl
    << "energy       = " << MET.energy() << std::endl
    << "significance = " << MET.significance() << std::endl;
};

void JetMonitor (const pat::Jet &Jet) {

    std::cout << "pt           = " << Jet.pt() << std::endl
    << "eta          = " << Jet.eta() << std::endl
    << "phi          = " << Jet.phi() << std::endl
    << "energy       = " << Jet.energy() << std::endl
    << "bprob        = " << Jet.bDiscriminator("pfDeepCSVJetTags:probb") + Jet.bDiscriminator("pfDeepCSVJetTags:probb") << std::endl;
};