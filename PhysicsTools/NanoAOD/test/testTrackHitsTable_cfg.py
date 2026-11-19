# TrackHitsTableProducer on the tracks of TrackHitsTestInputProducer, with the Run 3 geometry and conditions.
#   cmsRun testTrackHitsTable_cfg.py [aod]   ("aod": on tracks whose track extras are not in the event)
import sys
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2025_cff import Run3_2025
from PhysicsTools.NanoAOD.trackHitsTableProducer_cfi import trackHitsTableProducer

aod = "aod" in sys.argv[1:]
process = cms.Process("TEST", Run3_2025)
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2025_realistic", "")
process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi")
process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi")

process.source = cms.Source("EmptySource")
process.maxEvents.input = 1
process.input = cms.EDProducer("TrackHitsTestInputProducer")
tracks = "input:aod" if aod else "input"
process.trackTable = cms.EDProducer("SimpleTrackFlatTableProducer", src=cms.InputTag(tracks), cut=cms.string(""),
    name=cms.string("Track"), doc=cms.string(""), singleton=cms.bool(False), extension=cms.bool(False),
    variables=cms.PSet(pt=cms.PSet(expr=cms.string("pt"), type=cms.string("float"), doc=cms.string("pt"),
                                   precision=cms.int32(-1), lazyEval=cms.untracked.bool(False))))
process.hitsTable = trackHitsTableProducer.clone(src=tracks, name="Track", hitsName="TrackHits")
process.cutTrackTable = process.trackTable.clone(cut="pt > 15", maxLen=cms.uint32(1), name="CutTrack")
process.cutHitsTable = process.hitsTable.clone(cut="pt > 15", maxLen=cms.uint32(1), name="CutTrack",
                                               hitsName="CutTrackHits")
process.output = cms.OutputModule("NanoAODOutputModule",
    fileName=cms.untracked.string("testTrackHitsTableAOD.root" if aod else "testTrackHitsTable.root"),
    outputCommands=cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*_*_TEST"))
process.end = cms.EndPath(process.output, cms.Task(process.input, process.trackTable, process.hitsTable,
                                                   process.cutTrackTable, process.cutHitsTable))
