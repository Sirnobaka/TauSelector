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

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/PFTauDecayModeTools.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Common/interface/AssociativeIterator.h"

// For Btag calibration
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

// PileUp reweighting
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//edm::Handle<pat::TriggerObjectStandAloneCollection> &triggerObjects

std::string decimal_to_binary_string(UInt_t n) {
    std::string result = "";
    do {
        result = (std::to_string(n % 2)) + result;
        n = n / 2;
    } while (n > 0);
    return result;
}

UInt_t TriggerMatchingFunc (edm::Handle<pat::TriggerObjectStandAloneCollection> &triggerObjects,
                            edm::Handle<edm::TriggerResults> &triggerResults, 
                            const edm::TriggerNames &triggerNames,
                            double Pt, double Eta, double Phi, 
                            std::vector<std::string> &triggerNamesGroup,
                            bool monitoring) 
{

    int null = -10;

    UInt_t Particle_TriggerMatch = 0;
    //int Particle_TriggerObject = null;
    std::vector<int> Particle_TriggerObjects;
    bool Particle_TriggerObject_found = false;

    //// for commented monitoring
    ////if (monitoring) std::cout << "Number of trigger objects = " << triggerObjects->size() << std::endl;

    //double Particle_TriggerObject_dRmin = 100;
    //double Particle_triggerObject_dPt = 100;

    for (size_t i = 0; i < triggerObjects->size(); ++i) {
        pat::TriggerObjectStandAlone patTriggerObjectStandAloneUnpacked(triggerObjects->at(i));
        auto& TriggerObject_i = (*triggerObjects)[i];
        patTriggerObjectStandAloneUnpacked.unpackPathNames(triggerNames);
        double ParticledR = sqrt(deltaR2(Eta, Phi, TriggerObject_i.eta(), TriggerObject_i.phi()));
        double ParticledPt = abs(Pt - TriggerObject_i.pt());
        //
        if (TriggerObject_i.pt() > 15) {
            for(unsigned name_i = 0; name_i < patTriggerObjectStandAloneUnpacked.pathNames().size(); name_i++) {
                // Loop over set of trigger names
                for (unsigned int trigger_i=0; trigger_i<triggerNamesGroup.size(); ++trigger_i) {
                    if (patTriggerObjectStandAloneUnpacked.pathNames()[name_i].find(triggerNamesGroup[trigger_i].c_str())!= std::string::npos ) {
                        ////if (monitoring && ParticledR < 0.4) {
                        ////    std::cout << triggerNamesGroup[trigger_i] << " (dR = " << ParticledR << ")" << std::endl;
                        ////}
                        if (ParticledR < 0.1) {
                            //if (monitoring) std::cout << "dR  = " << ParticledR << std::endl;
                            ////if (monitoring) std::cout << "dPt = " << ParticledPt << std::endl;
                            Particle_TriggerObjects.push_back(i);
                            //Particle_TriggerObject_dRmin = ParticledR;
                            //Particle_triggerObject_dPt = ParticledPt;
                            //Particle_TriggerObject = i;
                            Particle_TriggerObject_found = true;
                        }
                    } else continue;
                }
            }    
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Add already picked trigger object. There may be multiple trigger names.
    /*
    if (monitoring) {
        std::cout << "Particle trigger object number = " << Particle_TriggerObject;
        if (Particle_TriggerObject_found) std::cout << " (Selected)";
        std::cout << std::endl
        << "(Eta, Phi) = " << (*triggerObjects)[Particle_TriggerObject].eta() << ", " << (*triggerObjects)[Particle_TriggerObject].phi() << std::endl
        << "Pt         = " << (*triggerObjects)[Particle_TriggerObject].pt() << std::endl;
    }
    */
    // If at least one trigger object has dR < 0.1 with particle
    if (Particle_TriggerObject_found) { 
        // iterate over all these triggerr objects
        for (auto triggerObject_iter = Particle_TriggerObjects.begin(); triggerObject_iter != Particle_TriggerObjects.end(); triggerObject_iter++) {
            // get particular trigger object
            pat::TriggerObjectStandAlone patTriggerObjectStandAloneUnpacked_Particle(triggerObjects->at(*triggerObject_iter));
            patTriggerObjectStandAloneUnpacked_Particle.unpackPathNames(triggerNames);
            // iterate over every name from this trigger object
            for(unsigned name_i = 0; name_i < patTriggerObjectStandAloneUnpacked_Particle.pathNames().size(); name_i++) {
                ////if (monitoring) std::cout << patTriggerObjectStandAloneUnpacked_Particle.pathNames()[name_i] << ", ";
                // match it with names from my trigger group
                for (unsigned int trigger_i=0; trigger_i<triggerNamesGroup.size(); ++trigger_i) {
                    int i_signed = trigger_i;
                    UInt_t thisTrigger = TMath::Power(2, i_signed);
                    if (patTriggerObjectStandAloneUnpacked_Particle.pathNames()[name_i].find(triggerNamesGroup[trigger_i].c_str())!= std::string::npos ) {
                        //std::cout << triggerNamesGroup[trigger_i] << std::endl;
                        UInt_t LogicAndResult = thisTrigger & Particle_TriggerMatch;
                        // Calculate Long Integer of trigger matching
                        if (monitoring) {
                            ////std::cout << std::endl << "this Trigger  = " << decimal_to_binary_string(TMath::Power(2, i_signed)) << std::endl;
                            //" (" << TMath::Power(2, i_signed) << ")" << std::endl;
                            //std::cout << "current Match = " << decimal_to_binary_string(Particle_TriggerMatch) <<
                            //" (" << Particle_TriggerMatch << ")" << std::endl;
                            //std::cout << "result of & = " << decimal_to_binary_string(LogicAndResult) <<
                            //" (" << LogicAndResult << ")" << std::endl;
                        }
                        // Check if this trigger name has been already added
                        if (LogicAndResult == 0) {
                            //if (monitoring) std::cout << "Trigger added" << std::endl;
                            Particle_TriggerMatch = Particle_TriggerMatch + TMath::Power(2, i_signed);
                        } else {
                            continue;
                        }
                    }
                }
            }
            ////if (monitoring) std::cout << std::endl;
        }
}

    ////if (monitoring) std::cout << "Final Trigger matching = " << decimal_to_binary_string(Particle_TriggerMatch) << std::endl;

    return Particle_TriggerMatch;

}