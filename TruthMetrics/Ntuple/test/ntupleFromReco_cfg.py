# Write the trackster-building truth ntuple from an existing GEN-SIM-RECO (step3) file, without
# re-running RECO. The input must hold the truth graph and hit index (DIGI with enableTruth) and the
# CLUE3DHigh tracksters (RECO with ticl_dev). The CLUE3DHigh mask is not kept in FEVTDEBUGHLT, so it is
# recomputed here from the stored layer clusters and initial mask. The raw CLUEstering assignment needs
# the (transient) SoA layer clusters: rerunAssignment=1 rebuilds them from the stored rechits and reruns
# CLUEstering (for the emulator's raw closure); otherwise it is skipped and closure uses the stored
# tracksters only. clue3d=1 also runs the production (non-ticl_dev) CLUE3D on the same layer clusters and mask
# and stores its tracksters as ts_ticlTrackstersCLUE3DProd (a second baseline for the tuning).
#   cmsRun ntupleFromReco_cfg.py inputFiles=file:step3.root out=truthNtuple.root maxEvents=-1 [clue3d=1]
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9
from Configuration.ProcessModifiers.enableTruth_cff import enableTruth
from Configuration.ProcessModifiers.ticl_dev import ticl_dev

options = VarParsing("analysis")
options.register("out", "truthNtuple.root", VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "output ntuple (TFileService)")
options.register("rerunAssignment", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "rebuild the SoA layer clusters and rerun the CLUEstering assignment")
options.register("clue3d", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "also run the production CLUE3D (ticlTrackstersCLUE3DProd)")
options.parseArguments()

# without ticl_dev the CLUE3DHigh step is the production CLUE3D (the stored CLUEstering tracksters are read)
process = cms.Process("NTUPLE", Phase2C26I13M9, enableTruth, *([] if options.clue3d else [ticl_dev]))
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load("Configuration.Geometry.GeometryExtendedRun4D128Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Accelerators_cff")
process.load("RecoHGCal.TICL.TICLGeom_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T35", "")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(1))

from RecoHGCal.TICL.CLUE3DHighStep_cff import filteredLayerClustersCLUE3DHigh
process.filteredLayerClustersCLUE3DHigh = filteredLayerClustersCLUE3DHigh.clone()

from TruthMetrics.Ntuple.customiseTruthNtuple import customiseTruthNtuple
process.schedule = cms.Schedule()
process = customiseTruthNtuple(process, options.out)
process.tracksterTruthNtuplizer.tracksters = [cms.InputTag("ticlTrackstersCLUE3DHigh", "", "RECO")]
process.truthNtuple_step.associate(cms.Task(process.filteredLayerClustersCLUE3DHigh))
if options.rerunAssignment:
    import RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff as _lc
    from RecoHGCal.TICL.CLUE3DHighStep_cff import ticlTrackstersCLUEsteringAssignment
    task = cms.Task()
    for det in ("EE", "HSi", "HSci"):
        for stage in ("hgcalSoARecHits", "hgcalCLUEstering", "hgcalSoALayerClusters"):
            setattr(process, stage + det, getattr(_lc, stage + det).clone())
            task.add(getattr(process, stage + det))
    process.ticlTrackstersCLUEsteringAssignment = ticlTrackstersCLUEsteringAssignment.clone()
    task.add(process.ticlTrackstersCLUEsteringAssignment)
    process.truthNtuple_step.associate(task)
else:
    process.tracksterTruthNtuplizer.clueAssignment = ""
if options.clue3d:
    from RecoHGCal.TICL.CLUE3DHighStep_cff import ticlTrackstersCLUE3DHigh as _clue3d
    from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer
    from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal
    if _clue3d.patternRecognitionBy.value() != "CLUE3D":
        raise RuntimeError("clue3d=1: CLUE3DHigh is not CLUE3D (ticl_dev applied?)")
    process.ticlLayerTileProducer = ticlLayerTileProducer.clone()
    process.ticlSeedingGlobal = ticlSeedingGlobal.clone()
    process.ticlTrackstersCLUE3DProd = _clue3d.clone()
    process.tracksterTruthNtuplizer.tracksters.append(cms.InputTag("ticlTrackstersCLUE3DProd"))
    process.truthNtuple_step.associate(cms.Task(process.ticlLayerTileProducer, process.ticlSeedingGlobal,
                                                process.ticlTrackstersCLUE3DProd))
