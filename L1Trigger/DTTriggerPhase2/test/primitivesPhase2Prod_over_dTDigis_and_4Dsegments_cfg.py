import FWCore.ParameterSet.Config as cms

process = cms.Process("L1DTTrigPhase2Prod")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.DTGeometryESModule.applyAlignment = False

process.load("Phase2L1Trigger.DTTrigger.dtTriggerPhase2PrimitiveDigis_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
#process.GlobalTag.globaltag = "90X_dataRun2_Express_v2"
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"

#Calibrate Digis
process.load("Phase2L1Trigger.CalibratedDigis.CalibratedDigis_cfi")
#process.CalibratedDigis.flat_calib = 325 #turn to 0 to use the DB  , 325 for JM and Jorge benchmark

#DTTriggerPhase2
process.load("Phase2L1Trigger.DTTrigger.dtTriggerPhase2PrimitiveDigis_cfi")
#process.dtTriggerPhase2PrimitiveDigis.trigger_with_sl = 3  #4 means SL 1 and 3
#for the moment the part working in phase2 format is the slice test
process.dtTriggerPhase2PrimitiveDigis.p2_df = 0
#process.dtTriggerPhase2PrimitiveDigis.filter_primos = True
#for debugging
process.dtTriggerPhase2PrimitiveDigis.pinta = True
process.dtTriggerPhase2PrimitiveDigis.min_phinhits_match_segment = 4
#process.dtTriggerPhase2PrimitiveDigis.debug = True

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
        #'file:/eos/user/c/carrillo/digis_segments_Run2016BSingleMuonRAW-RECO.root'
        'file:/tmp/carrillo/digis_segments_Run2016BSingleMuonRAW-RECO.root'
        )
                            )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *'),
                               fileName = cms.untracked.string('/tmp/carrillo/DTTriggerPhase2Primitives.root')
)

process.p = cms.Path(process.CalibratedDigis*process.dtTriggerPhase2PrimitiveDigis)
#process.this_is_the_end = cms.EndPath(process.out)






