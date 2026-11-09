"""Add the event-display export to a process: standalone (test/export_cfg.py) or inside a reco job
(cmsDriver.py ... --customise Fireworks/EventExport/customise.addEventExport)."""
import FWCore.ParameterSet.Config as cms


def addEventExport(process, presets=("calo", "hgcal", "reco"), fileName="evd_nano.root", withSim=True):
    process.load("Fireworks.EventExport.eventExport_cff")  # attaches the table producers to the process
    from Fireworks.EventExport.eventExport_cff import presets as known
    unknown = set(presets) - set(known)
    if unknown:
        raise ValueError(f"unknown event-export presets {sorted(unknown)}; known: {sorted(known)}")

    tasks = [known[p] for p in presets]
    if "hgcal" in presets:
        # TICL objects, tracks and associations from DPGAnalysis/HGCalNanoAOD. Its track tables propagate
        # tracks to the HGCal surfaces, which needs the propagator EventSetup modules.
        process.load("DPGAnalysis.HGCalNanoAOD.HGCalNanoAOD_cff")
        process.load("TrackingTools.MaterialEffects.MaterialPropagator_cfi")
        process.load("TrackingTools.MaterialEffects.OppositeMaterialPropagator_cfi")
        seq = process.hgcalNanoValidationSequence if withSim else process.hgcalNanoSequence
        process.eventExportHGCalPath = cms.Path(seq)

    process.eventExportTask = cms.Task(*tasks)
    process.eventExportOutput = cms.OutputModule("NanoAODOutputModule",
        fileName=cms.untracked.string(fileName),
        compressionAlgorithm=cms.untracked.string("ZLIB"),
        compressionLevel=cms.untracked.int32(4),
        outputCommands=cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table*_*_*"),
    )
    process.eventExportEndPath = cms.EndPath(process.eventExportOutput, process.eventExportTask)
    if hasattr(process, "schedule") and process.schedule is not None:
        if "hgcal" in presets:
            process.schedule.append(process.eventExportHGCalPath)
        process.schedule.append(process.eventExportEndPath)
    return process
