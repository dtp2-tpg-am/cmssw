import FWCore.ParameterSet.Config as cms

process = cms.Process("L1DTTrigPhase2Prod")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.DTGeometryESModule.applyAlignment = False

process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"

process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")

process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.debug = False
process.dtTriggerPhase2PrimitiveDigis.dumpMPsToFile = True
process.dtTriggerPhase2PrimitiveDigis.chi2Th = cms.untracked.double(0.16)

#scenario
process.dtTriggerPhase2PrimitiveDigis.scenario = 0 #0 is mc, 1 is data, 2 is slice test
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis"
process.CalibratedDigis.scenario = 0

process.dtTriggerPhase2BayesPrimitiveDigis = process.dtTriggerPhase2PrimitiveDigis.clone()
process.dtTriggerPhase2BayesPrimitiveDigis.grouping_code = 2 ## initial grouping
process.dtTriggerPhase2BayesPrimitiveDigis.PseudoBayesPattern.minNLayerHits = 4
process.dtTriggerPhase2BayesPrimitiveDigis.PseudoBayesPattern.minSingleSLHitsMax = 2 
process.dtTriggerPhase2BayesPrimitiveDigis.PseudoBayesPattern.minSingleSLHitsMin = 2 
process.dtTriggerPhase2BayesPrimitiveDigis.PseudoBayesPattern.minUncorrelatedHits = 3
process.dtTriggerPhase2BayesPrimitiveDigis.dumpMPsToFile = True

process.dtTriggerPhase2StdPrimitiveDigis   = process.dtTriggerPhase2PrimitiveDigis.clone()
process.dtTriggerPhase2StdPrimitiveDigis.grouping_code = 0 ## initial grouping
process.dtTriggerPhase2StdPrimitiveDigis.dumpMPsToFile = True

#process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:/eos/cms/store/user/folguera/P2L1TUpgrade/digis_segments_Run2016BSingleMuonRAW-RECO_camilo.root'))
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:/eos/cms/store/user/folguera/P2L1TUpgrade/Mu_FlatPt2to100-pythia8-gun_file.root'
)
    # 'file:/eos/cms/store/user/folguera/P2L1TUpgrade/digis_segments_Run2016BSingleMuonRAW-RECO_camilo.root'),skipEvents=cms.untracked.uint32(1)),
                        )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))


####################### SliceTest specials ##############################



process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring(
                                   'drop *',
                                   'keep *_CalibratedDigis_*_*',
                                   'keep *_dtTriggerPhase2BayesPrimitiveDigis_*_*',
                                   'keep *_dtTriggerPhase2StdPrimitiveDigis_*_*',
                                   'keep *_genParticles_*_*',
                               ),
                               fileName = cms.untracked.string('DTTriggerPhase2Primitives.root')
)

process.p = cms.Path(process.CalibratedDigis*process.dtTriggerPhase2BayesPrimitiveDigis*process.dtTriggerPhase2StdPrimitiveDigis)
process.this_is_the_end = cms.EndPath(process.out)






