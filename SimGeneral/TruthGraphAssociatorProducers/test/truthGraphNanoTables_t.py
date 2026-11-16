#!/usr/bin/env python3
"""truthGraphNanoTables_cff: the module labels, the table names, the inputs of the association
tables and the content of the tasks."""

import unittest

import FWCore.ParameterSet.Config as cms

from SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociationLabels_cff import recoLabels

_domains = ("tracksters", "tracks", "pfClustersEcal", "pfClustersHcal", "vertices", "secondaryVertices")


def _process():
    process = cms.Process("TEST")
    process.load("SimGeneral.TruthGraphAssociatorProducers.truthGraphNanoTables_cff")
    return process


class TestTruthGraphNanoTables(unittest.TestCase):
    def testAssociationTables(self):
        process = _process()
        from SimGeneral.TruthGraphAssociatorProducers.truthGraphNanoTables_cff import truthGraphAssociationTables
        labels = [label for domain in _domains for label in recoLabels(domain, "offline")]
        self.assertEqual(len(truthGraphAssociationTables), 2 * len(labels))
        for name, module in truthGraphAssociationTables.items():
            self.assertTrue(name.endswith("Table"))
            self.assertNotIn("_", name)
            self.assertNotIn(":", name)
            self.assertTrue(hasattr(process, name))
            vertex = module.type_() == "TruthBranchVertexAssociationFlatTableProducer"
            self.assertEqual(module.targets.getProductInstanceLabel(), "selectedVertices" if vertex else "selected")
            suffix = "RecoToTruthAdaptiveNominal" if module.direction.value() == "recoToTruth" else "TruthToReco"
            self.assertTrue(module.src.getProductInstanceLabel().endswith(suffix))
            self.assertTrue(module.name.value()[0].isupper())
        # the cfi templates are imported under private names, so process.load adds no stray module
        for template in ("truthGraphFlatTable", "truthBranchAssociationFlatTable", "truthBranchVertexAssociationFlatTable"):
            self.assertFalse(hasattr(process, template))

    def testTrackPayload(self):
        process = _process()
        for name in ("truthGraphGeneralTracksToTruthTable", "truthGraphTruthToGeneralTracksTable"):
            module = getattr(process, name)
            for payload in (module.recoToTruthPayload, module.truthToRecoPayload):
                self.assertEqual(payload.name.value(), "nSharedHits")
        self.assertEqual(process.truthGraphParticleFlowClusterECALToTruthTable.recoToTruthPayload.name.value(), "")

    def testSelection(self):
        process = _process()
        self.assertEqual(process.truthGraphTable.keepLevelMask.value(), 0)
        self.assertEqual(process.truthGraphCaloBoundaryTable.keepLevelMask.value(), 1 << 3)
        particleMaps = [label for domain in _domains[:4] for label in recoLabels(domain, "offline")]
        self.assertEqual(len(process.truthGraphTable.keepMatchedBy), len(particleMaps))
        self.assertEqual(len(process.truthGraphTable.keepMatchedByTruthToReco), len(particleMaps))
        vertexMaps = [label for domain in _domains[4:] for label in recoLabels(domain, "offline")]
        self.assertEqual(len(process.truthGraphTable.keepVerticesMatchedBy), len(vertexMaps))
        self.assertEqual(len(process.truthGraphTable.keepVerticesMatchedByTruthToReco), len(vertexMaps))

    def testTasks(self):
        process = _process()
        task = process.truthGraphNanoTablesTask
        self.assertTrue(task.contains(process.truthGraphTable))
        self.assertTrue(task.contains(process.truthBranchTargets))
        self.assertTrue(task.contains(process.allTrackToTruthBranchAssociators))
        self.assertFalse(task.contains(process.truthGraphCaloBoundaryTable))
        # The graph tables alone: no associator, so they run on any input.
        tables = process.truthGraphTablesTask
        self.assertTrue(tables.contains(process.truthGraphTable))
        self.assertFalse(tables.contains(process.truthBranchTargets))
        self.assertFalse(tables.contains(process.allTrackToTruthBranchAssociators))
        # The opt-in replaces truthGraphTable in place, so the other tables keep reading it.
        process.truthGraphTable = process.truthGraphCaloBoundaryTable.clone()
        self.assertTrue(task.contains(process.truthGraphTable))
        self.assertEqual(process.truthGraphTable.keepLevelMask.value(), 1 << 3)


if __name__ == "__main__":
    unittest.main()
