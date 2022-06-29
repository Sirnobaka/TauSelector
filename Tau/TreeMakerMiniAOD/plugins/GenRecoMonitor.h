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
#include "DataFormats/JetReco/interface/GenJet.h"

#include "Tau/TreeMakerMiniAOD/plugins/MySimpleParticle.h"

#include <cstdio>
#include <iostream>
#include <math.h>

static inline double sqr(double x) {
    return x * x;
}

static double dphi(double phi1, double phi2) {
    double result = TMath::Abs(phi1 - phi2);
    if (result < TMath::Pi()) return result;
    return 2 * TMath::Pi() - result;
}

class GenRecoMonitor {
public:

    const reco::GenParticle *genParticle = nullptr;
    const reco::GenJet *genJet = nullptr;
    int    recoPdgId;
    int    recoHadronFlavour;
    double recoPt;
    double recoEta;
    double recoPhi;
    double recoE;
    double deltaR;
    double null = -10;
    //const math::XYZPoint pv_position;
    //int    recoCharge;
    /*
    double genPt;
    double genEta;
    double genPhi;
    int    genPdgId;
    */

    GenRecoMonitor();
    GenRecoMonitor(const reco::GenParticle &genparticle, int recoparticlepdgId, double recoparticlePt, double recoparticleEta, double recoparticlePhi)
        : genParticle(&genparticle), recoPdgId(recoparticlepdgId), recoPt(recoparticlePt), recoEta(recoparticleEta), recoPhi(recoparticlePhi)
    {
        deltaR = TMath::Sqrt( sqr(dphi(genParticle->phi(), recoPhi)) + sqr(genParticle->eta() - recoEta) );
        //std::cout << "reco::GenParticle" << std::endl;
    };
    GenRecoMonitor(const reco::GenParticle &genparticle)
        : genParticle(&genparticle), recoPdgId(0), recoPt(null), recoEta(null), recoPhi(null), deltaR(0)
    {
        //std::cout << "reco::GenParticle" << std::endl;
    };
    // reco::Candidate constructor
    GenRecoMonitor(const reco::Candidate &genparticle, int recoparticlepdgId, double recoparticlePt, double recoparticleEta, double recoparticlePhi)
        : recoPdgId(recoparticlepdgId), recoPt(recoparticlePt), recoEta(recoparticleEta), recoPhi(recoparticlePhi)
    {
        //std::cout << "reco::Candidate" << std::endl;
        genParticle = static_cast<const reco::GenParticle*>(&genparticle);
        deltaR = TMath::Sqrt( sqr(dphi(genParticle->phi(), recoPhi)) + sqr(genParticle->eta() - recoEta) );
        //std::cout << "casted" << std::endl;
    };
    GenRecoMonitor(const reco::Candidate &genparticle)
        : recoPdgId(0), recoPt(null), recoEta(null), recoPhi(null), deltaR(0)
    {
        //std::cout << "reco::Candidate" << std::endl;
        genParticle = static_cast<const reco::GenParticle*>(&genparticle);
        //std::cout << "casted" << std::endl;
    };
    // for genjets
    GenRecoMonitor(const reco::GenJet &genjet, int recojetHadronFlavour, double recojetPt, double recojetEta, double recojetPhi, double recojetE)
        : genJet(&genjet), recoHadronFlavour(recojetHadronFlavour), recoPt(recojetPt), recoEta(recojetEta), recoPhi(recojetPhi), recoE(recojetE)
    {
        deltaR = TMath::Sqrt( sqr(dphi(genJet->phi(), recoPhi)) + sqr(genJet->eta() - recoEta) );
    };
    void PrintParameters();
    void PrintDaughters();
    void PrintMothers();
    void PrintComp(bool mothers, bool daughters);
    void PrintGen(bool mothers, bool daughters);
    void PrintCompJets();
    virtual ~GenRecoMonitor();

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of GenRecoMonitor object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
GenRecoMonitor::GenRecoMonitor() {

};

GenRecoMonitor::~GenRecoMonitor() {
// Destructor
};

void GenRecoMonitor::PrintDaughters() {
    std::cout << "Daughters pdgIds: ";
    for (unsigned i = 0; i < genParticle->numberOfDaughters(); ++i) {
        const reco::Candidate* daughter = genParticle->daughter(i);
        std::cout << daughter->pdgId() << "(" << i << "), " ;
    }
    std::cout << std::endl;
};

void GenRecoMonitor::PrintMothers() {
    std::cout << "Mothers of " << genParticle->pdgId() << " <-- ";
    int nMothers = 0;
    for (auto mother = genParticle->mother(); mother; mother = mother->mother()) {
        std::cout << "|" << mother->pdgId() << "(" << nMothers << ") <-- ";
        nMothers++;
    }
    std::cout << std::endl;
};

void GenRecoMonitor::PrintComp(bool mothers, bool daughters) {
    std::cout << "pdgId Gen (Reco) = " << genParticle->pdgId() << " (" << recoPdgId << ")" << std::endl;
    std::cout << "status           = " << genParticle->status() << std::endl;
    std::cout << "Pt Gen (Reco)    = " << genParticle->pt() << " (" << recoPt << ")" << std::endl;
    std::cout << "delta R  = " << deltaR << std::endl;
    std::cout << "isDirectPromptTauDecayProduct = " << genParticle->statusFlags().isDirectPromptTauDecayProduct() << std::endl;
    std::cout << "isPrompt = " << genParticle->statusFlags().isPrompt() << std::endl;
    std::cout << "isDecayedLeptonHadron = " << genParticle->statusFlags().isDecayedLeptonHadron() << std::endl;
    std::cout << "isTauDecayProduct = " << genParticle->statusFlags().isTauDecayProduct() << std::endl;
    std::cout << "isPromptTauDecayProduct = " << genParticle->statusFlags().isPromptTauDecayProduct() << std::endl;
    std::cout << "isDirectTauDecayProduct = " << genParticle->statusFlags().isDirectTauDecayProduct() << std::endl;
    std::cout << "isDirectPromptTauDecayProduct = " << genParticle->statusFlags().isDirectPromptTauDecayProduct() << std::endl;
    std::cout << "isDirectHadronDecayProduct = " << genParticle->statusFlags().isDirectHadronDecayProduct() << std::endl;
    if (mothers) PrintMothers();
    if (daughters) PrintDaughters();
};

void GenRecoMonitor::PrintCompJets() {
    std::cout << "pdgId Gen (Reco)  = " << genJet->pdgId() << " (" << recoHadronFlavour << ")" << std::endl;
    //std::cout << "status           = " << genJet->status() << std::endl;
    std::cout << "Pt Gen (Reco)     = " << genJet->pt() << " (" << recoPt << ")" << std::endl;
    std::cout << "Energy Gen (Reco) = " << genJet->energy() << " (" << recoE << ")" << std::endl;
    std::cout << "delta R           = " << deltaR << std::endl;
};

void GenRecoMonitor::PrintGen(bool mothers, bool daughters) {
    std::cout << "pdgId Gen = " << genParticle->pdgId() << std::endl;
    std::cout << "status    = " << genParticle->status() << std::endl;
    if (mothers) PrintMothers();
    if (daughters) PrintDaughters();
};

bool DecaychannelMatch(std::vector<MySimpleParticle> &particles, int p1, int p2 = 0, int p3 = 0, int p4 = 0, int p5 = 0, int p6 = 0) {
    // Copy pdgid-s of all particles
    std::vector<int> list;
    
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
        for (unsigned int j = 0;j < list.size(); j++) {
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

    std::vector<MySimpleParticle> newList;
    
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
};

// Decay mode of Generated tau lepton

int GenTauDecayMode (const reco::GenParticle &genParticle) {

    int TauDM = -20;

    int tau_pdgid = genParticle.pdgId();

    std::vector<MySimpleParticle> tau_daughters;
    for (unsigned k = 0; k < genParticle.numberOfDaughters(); k++) {
        MySimpleParticle tp(genParticle.daughter(k)->p4().Px(), genParticle.daughter(k)->p4().Py(), genParticle.daughter(k)->p4().Pz(), genParticle.daughter(k)->p4().E(), genParticle.daughter(k)->pdgId());
        tau_daughters.push_back(tp);
    }
    if (tau_daughters.size() >= 2 && (abs(tau_daughters[0].pdgid()) == 15 || abs(tau_daughters[1].pdgid()) == 15)) return -100;

    std::vector<int> pdgid;
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
        //channel = 3;
        //gentau_nPi0 = 0;
        TauDM = 0;
        //if (abs(pdgid[1]) == 321) channel = 6;
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
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 321, 111) ) ||
                // new lines with neutral kaons
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-311) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 311) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-321, 311) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 321,-311) )
            )
        ) {
        //channel = 4;
        //gentau_nPi0 = 1;
        TauDM = 1;
        //if(abs(pdgid[1])==321 || abs(pdgid[2])==130 || abs(pdgid[2])==310) channel = 7;
    }
    else if( pdgid.size() == 4 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 211) ) ||
                // new lines with  kaons
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters, 16, 111,-211,-311) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 211, 311) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 311,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111,-311, 321) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 130, 310,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-130,-310, 211) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 130, 130,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-130,-130, 211) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111,-311,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 311, 211) )

            )
        ) {
        //channel = 5;
        //gentau_nPi0 = 2;
        TauDM = 2;
    }
    else if( pdgid.size() == 5 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 211) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 321) )
            )
        ) {
        TauDM = 3;
    }
    else if( pdgid.size() == 6 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111, 111,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 111, 211) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111, 111,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 111, 321) )
            )
        ) {
        TauDM = 4;
    }
    else if( pdgid.size() == 4 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211) ) ||
                // new DMs (with Kaons)
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 211,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211,-211, 321) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 321, 321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211,-321,-321) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 321,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 321,-321) )            )
        ) {
        //channel = 6;
        //gentau_nPi0 = 0;
        TauDM = 10;
    }
    else if( pdgid.size()==5 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211, 111) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 111) ) ||
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211,-311) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 311) )
            )
        ) {
        //channel = 8;
        //gentau_nPi0 = 1;
        TauDM = 11;
    }
    else if( pdgid.size()==5 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-211) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 211) ) ||
                // new kaons
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-321) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 321) )
            )
        ) {
        //channel = 9;
        //gentau_nPi0 = 3;
        TauDM = 3;
    } 
    else if( pdgid.size()==6 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211, 111, 111) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 111, 111) )
            )
        ) {
        // 3 Pi charged, 2 Pi zero
        TauDM = 12;
    }
    else if( pdgid.size()==3 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12) )
            )
        ) {
        //channel = 1;
        TauDM = -1;
    }
    else if( pdgid.size()==4 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12, 22) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12, 22) )
            )
        ) {
        //channel = 1;
        TauDM = -1;
    }
    else if( pdgid.size()==3 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14) )
            )
        ) {
        //channel = 2;
        TauDM = -2;
    }
    else if( pdgid.size()==4 &&
            (
                ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14, 22) ) ||
                ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14, 22) )
            )
        ) {
        //channel = 2;
        TauDM = -2;
    }

    return TauDM;
}
