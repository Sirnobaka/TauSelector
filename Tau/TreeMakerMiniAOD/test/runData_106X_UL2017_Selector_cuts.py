import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TauGenMCanalyzer")

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

# Single Electron 2016 17 Jul 2018 data
# For mc 2017 94X
#process.GlobalTag.globaltag='94X_mc2017_realistic_v14'
# For mc 2017 UL 106X (Summer20)
process.GlobalTag.globaltag='106X_mc2017_realistic_v9'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
    #'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
    #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    #'/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/F296B7FF-629F-E811-86FE-F4E9D4AF0AF0.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/90000/E40B3455-5939-E811-94C5-0CC47A7C34C8.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/E2A2C81B-0638-E811-89CE-008CFAC93EA8.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2017D/Tau/MINIAOD/31Mar2018-v1/00000/12F092A0-3F37-E811-AD27-7845C4F92C96.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2016C/SingleMuon/MINIAOD/17Jul2018-v1/20000/065D2BF3-9198-E811-844D-90E2BAC9B7A8.root'
    #'/store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    #'/store/mc/RunIISummer20UL17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v1/00000/1BDD5344-60CA-9642-9C68-5F7ED5235DA7.root'
    '/store/data/Run2017F/SingleMuon/MINIAOD/UL2017_MiniAODv2-v1/270000/083FA4C6-5737-A74B-BF56-0FF8F9DBD00F.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
)

#process.treemaker = cms.EDAnalyzer('TauGenMCanalyzer')

# MVA MET

#process.load("RecoJets.JetProducers.ak4PFJets_cfi")
#process.task.add(process.ak4PFJets)
#process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
#process.ak4PFJets.doAreaFastjet = cms.bool(True)

#from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

#process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
#process.task.add(process.pfMVAMEtTask)
#process.MVAMET = process.pfMVAMEtTask
#process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
#process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
#process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Id
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
    updatedTauName = updatedTauName,
    toKeep = [ "MVADM_2017_v1" ]) #other tauIDs can be added in
tauIdEmbedder.runTauID()
# Path and EndPath definitions
# Was not commented
#process.p = cms.Path(
#    process.rerunMvaIsolationSequence *
#    getattr(process,updatedTauName)
#)

# EGamma POG corrections (not sure in era is OK for MC samples we use now)
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
sys.path.append('/afs/cern.ch/user/a/aoskin/Tau_packeges/CMSSW_10_2_22/src/EgammaUser/EgammaPostRecoTools/python/')
from EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       applyEnergyCorrections=False,
                       applyVIDOnCorrectedEgamma=False,
                       isMiniAOD=True,
                       era="2017-UL",
                       #eleIDModules=_defaultEleIDModules,
                       #phoIDModules=_defaultPhoIDModules,
                       runVID=True, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       runEnergyCorrections=True,
                       applyEPCombBug=False,
                       autoAdjustParams=True,
                       computeHeepTrkPtIso=True)   
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

# Update jet collection
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
   btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg'
      ],
   postfix='NewDFTraining'
)

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)

# TreeMakder for miniAOD

process.TauGenMCanalyzer = cms.EDAnalyzer("TauGenMCanalyzer",
    isMC = cms.bool(False),
    noRejection = cms.bool(False),
    noTauSelection = cms.bool(False),
    noLepSelection = cms.bool(False),
    SelectBquarks = cms.bool(False),
    BtagThr = cms.double(-10), # loose DeepFlavour
    #useHLT = cms.bool(False),
    #useTargetHLT = cms.bool(False),
    tauPtMin = cms.double(20),
    BJetPtMin = cms.double(30),
    MuElePtMin = cms.double(20),
    EtaMax = cms.double(2.4),
    piPtMin = cms.double(0),
    tauEtaMax = cms.double(2.4),
    tauDzMax = cms.double(0.2),
    METcut = cms.double(0),
    null = cms.double(-10),
    # Parameters added at January 2021
    JetEtaMax = cms.double(2.4),
    TriggerMatchingOn = cms.bool(False),
    #NrequiredJets = cms.int32(2),
    #NrequiredBJets = cms.int32(1),
    #NrequiredLeptons = cms.int32(1),
    #requiredLeptonPDGID = cms.int32(13),
    UseTau = cms.bool(True),
    UsePuppiJets = cms.bool(False),
    ####
    #tauCollection = cms.string("slimmedTaus"),
    tauCollection = cms.string("slimmedTausNewID"),
    muonCollection = cms.string("slimmedMuons"),
    electronCollection = cms.string("slimmedElectrons"),
    jetCollection = cms.string("selectedUpdatedPatJetsNewDFTraining"), # slimmedJets
    PuppijetCollection = cms.string("slimmedJetsPuppi"),
    metCollection = cms.string("slimmedMETs"),
    # Puppi MET
    PuppimetCollection = cms.string("slimmedMETsPuppi"),
    vertexCollection = cms.string("offlineSlimmedPrimaryVertices"),
    SVCollection     = cms.string("slimmedSecondaryVertices"),
    #genParticleCollection = cms.string("genParticles"),
    #genParticleCollection = cms.string("packedGenParticles"),
    genParticleCollection = cms.string("prunedGenParticles"),
    GenJetsCollection = cms.string("slimmedGenJets"),
    trackCollection = cms.string("isolatedTracks"), # generalTracks for AOD
    PackedCandidateCollection = cms.string("packedPFCandidates"),
    GenEventInfo = cms.string("generator"),
    LHEEventInfo = cms.string("externalLHEProducer"),
    rhoTag = cms.InputTag("fixedGridRhoAll"),
    effAreasConfigFile = cms.FileInPath("Tau/TreeMakerMiniAOD/data/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
    PU_data_file = cms.FileInPath("Tau/TreeMakerMiniAOD/data/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root"),
    PU_MC_file = cms.FileInPath("Tau/TreeMakerMiniAOD/data/PU_MC_UL2017.root"),
    # for 2017 data
    Triggerobjects = cms.InputTag('slimmedPatTrigger', '', 'PAT'),
    prescales      = cms.InputTag('patTrigger', '', 'PAT'),
    prescalesL1min = cms.InputTag('patTrigger', 'l1min', 'PAT'),
    prescalesL1max = cms.InputTag('patTrigger', 'l1max', 'PAT'),
    # pileup info
    PileupInfo = cms.InputTag('slimmedAddPileupInfo'),
    # for 2016 data
    #Triggerobjects = cms.InputTag('slimmedPatTrigger', '', 'DQM'),
    #prescales      = cms.InputTag('patTrigger', '', 'DQM'),
    #prescalesL1min = cms.InputTag('patTrigger', 'l1min', 'DQM'),
    #prescalesL1max = cms.InputTag('patTrigger', 'l1max', 'DQM'),
    Triggers = cms.vstring(""),
    TriggerTarget1 = cms.vstring(
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",
    ),
    TriggerTarget2 = cms.vstring(
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
    ),
    TriggerTarget3 = cms.vstring(
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
    ),
    # SingleMuon 2017 part 1
    Triggers1 = cms.vstring(
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_TightID_SingleL1_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_TightID_SingleL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_TightID_SingleL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_v",
        "HLT_IsoMu24_v",
        "HLT_IsoMu27_v",
        "HLT_IsoMu30_v",
        "HLT_L1SingleMu18_v",
        "HLT_L1SingleMu25_v",
    ),
    # SingleMuon 2017 part 2
    Triggers2 = cms.vstring(
        "HLT_L1_DoubleJet30_Mass_Min400_Mu10_v",
        "HLT_L2Mu10_v",
        "HLT_L2Mu50_v",
        "HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v",
        "HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v",
        "HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v",
        "HLT_Mu15_IsoVVVL_PFHT450_v",
        "HLT_Mu15_IsoVVVL_PFHT600_v",
        "HLT_Mu20_v",
        "HLT_Mu27_v",
        "HLT_Mu3_PFJet40_v",
        "HLT_Mu50_IsoVVVL_PFHT450_v",
        "HLT_Mu50_v",
        "HLT_Mu55_v",
        "HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v",
        "HLT_OldMu100_v",
        "HLT_TkMu100_v",
    ),
    # SingleElectron 2017 part 1
    Triggers3 = cms.vstring(
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele135_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele145_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v",
        "HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v",
        "HLT_Ele15_IsoVVVL_PFHT450_v",
        "HLT_Ele15_IsoVVVL_PFHT600_v",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",
        "HLT_Ele200_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele250_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele27_WPTight_Gsf_v",
    ),
    # SingleElectron 2017 part 2
    Triggers4 = cms.vstring(
        "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",
        "HLT_Ele300_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v",
        "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
        "HLT_Ele32_WPTight_Gsf_v",
        "HLT_Ele35_WPTight_Gsf_L1EGMT_v",
        "HLT_Ele35_WPTight_Gsf_v",
        "HLT_Ele38_WPTight_Gsf_v",
        "HLT_Ele40_WPTight_Gsf_v",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v",
        "HLT_Ele50_IsoVVVL_PFHT450_v",
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v",
    ),
    # Tau
    Triggers5 = cms.vstring(
        "HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
        "HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
        "HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
        "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
        "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
        "HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
        "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v",
        "HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1_v",
        "HLT_MediumChargedIsoPFTau100HighPtRelaxedIso_Trk50_eta2p1_1pr_v",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v",
        "HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v",
        "HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
        "HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
        "HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
    ),
    # MET 2017 part 1
    Triggers6 = cms.vstring(
        "HLT_CaloMET100_HBHECleaned_v",
        "HLT_CaloMET100_NotCleaned_v",
        "HLT_CaloMET110_NotCleaned_v",
        "HLT_CaloMET250_HBHECleaned_v",
        "HLT_CaloMET250_NotCleaned_v",
        "HLT_CaloMET300_HBHECleaned_v",
        "HLT_CaloMET350_HBHECleaned_v",
        "HLT_CaloMET70_HBHECleaned_v",
        "HLT_CaloMET80_HBHECleaned_v",
        "HLT_CaloMET80_NotCleaned_v",
        "HLT_CaloMET90_HBHECleaned_v",
        "HLT_CaloMET90_NotCleaned_v",
        "HLT_CaloMHT90_v",
        "HLT_DiJet110_35_Mjj650_PFMET110_v",
        "HLT_DiJet110_35_Mjj650_PFMET120_v",
        "HLT_DiJet110_35_Mjj650_PFMET130_v",
        "HLT_L1ETMHadSeeds_v",
        "HLT_MET105_IsoTrk50_v",
        "HLT_MET120_IsoTrk50_v",
        "HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
        "HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_v",
        "HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_v",
        "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v",
        "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",
        "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v",
        "HLT_PFMET110_PFMHT110_IDTight_v",
        "HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v",
    ),
    # MET 2017 part 2
    Triggers7 = cms.vstring(
        "HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v",
        "HLT_PFMET130_PFMHT130_IDTight_v",
        "HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v",
        "HLT_PFMET140_PFMHT140_IDTight_v",
        "HLT_PFMET200_HBHECleaned_v",
        "HLT_PFMET200_HBHE_BeamHaloCleaned_v",
        "HLT_PFMET200_NotCleaned_v",
        "HLT_PFMET250_HBHECleaned_v",
        "HLT_PFMET300_HBHECleaned_v",
        "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
        "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",
        "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",
        "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v",
        "HLT_PFMETTypeOne110_PFMHT110_IDTight_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_HFCleaned_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_v",
        "HLT_PFMETTypeOne130_PFMHT130_IDTight_v",
        "HLT_PFMETTypeOne140_PFMHT140_IDTight_v",
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET110_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET120_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET130_v",
    ), 
    SelectedTriggers = cms.vstring(
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v",
        "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v",
        "HLT_IsoMu24_v",
        "HLT_IsoMu27_v",
        "HLT_IsoMu30_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele27_WPTight_Gsf_v",
        "HLT_Ele32_WPTight_Gsf_v",
        "HLT_Ele35_WPTight_Gsf_v",
        "HLT_Ele38_WPTight_Gsf_v",
        "HLT_Ele40_WPTight_Gsf_v",
        #"HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
        #"HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
        #"HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
        #"HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
        #"HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
        #"HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
        #"HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
        #"HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
        "HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1_v",
        "HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1_v",
        #"HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
        #"HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
        #"HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
    ),
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Data_SingleMuon_106X_UL2017_cuts.root')
)

process.p = cms.Path(process.jecSequence * process.rerunMvaIsolationSequence * getattr(process,updatedTauName) * process.egammaPostRecoSeq * process.TauGenMCanalyzer)
#process.p = cms.Path(process.rerunMvaIsolationSequence * getattr(process,updatedTauName) * process.egammaPostRecoSeq * process.TauGenMCanalyzer)
#process.p = cms.Path(process.rerunMvaIsolationSequence * getattr(process,updatedTauName) * process.TauGenMCanalyzer)
process.p.associate(process.patAlgosToolsTask)
