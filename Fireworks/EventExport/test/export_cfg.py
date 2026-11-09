"""Standalone event-display export from RECO files.

cmsRun export_cfg.py inputFiles=file:step3.root outputFile=evd_nano.root presets=calo,hgcal,reco \
    geometry=ExtendedRun4D128 era=Phase2C26I13M9 globalTag=auto:phase2_realistic_T35
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

opts = VarParsing("analysis")
opts.register("presets", "calo,hgcal,reco", VarParsing.multiplicity.singleton, VarParsing.varType.string, "comma-separated presets")
opts.register("geometry", "ExtendedRun4D128", VarParsing.multiplicity.singleton, VarParsing.varType.string, "Configuration.Geometry.Geometry<X>Reco_cff")
opts.register("era", "Phase2C26I13M9", VarParsing.multiplicity.singleton, VarParsing.varType.string, "era of the input")
opts.register("globalTag", "auto:phase2_realistic_T35", VarParsing.multiplicity.singleton, VarParsing.varType.string, "global tag")
opts.register("withSim", True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "include sim objects and associations")
opts.register("threads", 4, VarParsing.multiplicity.singleton, VarParsing.varType.int, "number of threads")
opts.setDefault("outputFile", "evd_nano.root")
opts.parseArguments()

process = cms.Process("EVDEXPORT", getattr(eras, opts.era))
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.load(f"Configuration.Geometry.Geometry{opts.geometry}Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, "")

process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(opts.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(opts.maxEvents))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(opts.threads), wantSummary=cms.untracked.bool(True))

from Fireworks.EventExport.customise import addEventExport
addEventExport(process, presets=opts.presets.split(","), fileName=opts.outputFile, withSim=opts.withSim)
