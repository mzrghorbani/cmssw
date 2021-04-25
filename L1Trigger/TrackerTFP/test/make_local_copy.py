import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the source
process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring (
                                 "/store/relval/CMSSW_11_1_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v2_2026D49PU200_ext1-v1/20000/D3629C85-EA34-C147-AC4D-939C41DEC68A.root"
                             )
                             )

# tell the process to only run over 100 events (-1 would mean run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32 (1)
)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string (
                                   "D3629C85-EA34-C147-AC4D-939C41DEC68A.root"
                               )
                               )

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)