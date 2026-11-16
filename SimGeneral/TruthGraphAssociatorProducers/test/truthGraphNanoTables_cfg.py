# Constructs every module of truthGraphNanoTables_cff (no event), so each configuration is checked
# against the module's fillDescriptions.
import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.source = cms.Source("EmptySource")
process.maxEvents.input = 0
process.load("SimGeneral.TruthGraphAssociatorProducers.truthGraphNanoTables_cff")
process.truthGraphCaloBoundaryTask = cms.Task(process.truthGraphCaloBoundaryTable)
process.p = cms.Path(cms.Sequence(), process.truthGraphNanoTablesTask, process.truthGraphCaloBoundaryTask)
