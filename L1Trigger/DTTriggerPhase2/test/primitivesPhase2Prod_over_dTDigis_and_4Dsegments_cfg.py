import FWCore.ParameterSet.Config as cms

process = cms.Process("L1DTTrigPhase2Prod")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.DTGeometryESModule.applyAlignment = False

process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
#process.GlobalTag.globaltag = "90X_dataRun2_Express_v2"
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"

#Calibrate Digis
process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
#process.CalibratedDigis.flat_calib = 325 #turn to 0 to use the DB  , 325 for JM and Jorge benchmark

#DTTriggerPhase2
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.scenario = 0 #0 is mc, 1 is data, 2 is slice test

process.dtTriggerPhase2BayesPrimitiveDigis = process.dtTriggerPhase2PrimitiveDigis.clone()
process.process.dtTriggerPhase2BayesPrimitiveDigis.grouping_code = 2 ## initial grouping

process.dtTriggerPhase2StdPrimitiveDigis   = process.dtTriggerPhase2PrimitiveDigis.clone()
process.dtTriggerPhase2StdPrimitiveDigis.grouping_code = 0 ## initial grouping


process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
        #'file:/eos/user/c/carrillo/digis_segments_Run2016BSingleMuonRAW-RECO.root'
    '/store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/PU200_106X_upgrade2023_realistic_v3-v2/70000/F42F882F-B3A8-4346-870D-3E62C930D076.root'
        )
                            )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring(
                                   'keep *'
                               ),
                               fileName = cms.untracked.string('DTTriggerPhase2Primitives.root')
)

process.p = cms.Path(process.CalibratedDigis*process.dtTriggerPhase2PrimitiveDigis)
process.this_is_the_end = cms.EndPath(process.out)






