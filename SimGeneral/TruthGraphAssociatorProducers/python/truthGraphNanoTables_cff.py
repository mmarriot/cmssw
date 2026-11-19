# NanoAOD tables of the MC-truth graph (truth::Graph, PhysicsTools/TruthInfo): its particles,
# vertices and parent -> child links, and the reco <-> truth associations of this package in
# both directions.
#
# truthGraphTable (TruthGraphFlatTableProducer) writes
#   TruthParticles         the signal interaction, every particle a reco object is matched to
#                          (reco -> truth and truth -> reco maps), and their ancestors;
#   TruthVertices          the production vertex and the endpoint of each, and every truth vertex
#                          the vertex maps match (instance "vertices");
#   TruthParticleChildren  parent -> child links between them (instance "children");
# and the graph ids of the particle and vertex rows ("selected", "selectedVertices"), which the
# association tables read so their link indices are rows of these tables.
#
# truthGraphCaloBoundaryTable is the opt-in for a hit-level truth: the same selection plus every
# particle that entered the calorimeter (truth::LevelFlag CaloBoundary, pile-up included), which
# carry 99.8% of the simulated energy in the HGCAL rechits. Use it in place of truthGraphTable,
# so the other tables keep reading truthGraphTable:
#   process.truthGraphTable = process.truthGraphCaloBoundaryTable.clone()
#
# Association tables (TruthBranchNanoTableProducers.cc, the layout of truthNanoTables_cff): for
# the reco collections recoLabels names in the domains tracks, tracksters, pfClustersEcal,
# pfClustersHcal, vertices and secondaryVertices (not candidates and layerClusters),
# <Reco>ToTruth with a row per reco object and TruthTo<Reco> with a row per TruthParticles or
# TruthVertices row, each with n<A>Links / o<A>Links and a <A>Links table of index, graphId,
# payload and score. Every entry of the maps has its particle or vertex exported, so no link is
# left out.
#
# Tasks:
#   truthGraphTablesTask      the graph tables alone, no associator: the signal interaction and
#                             its ancestors (plus the keepLevelMask levels). Runs on any input; a
#                             file without the graph gives empty tables with all their columns.
#                             truthGraphTable still reads the maps it is configured with: if
#                             another Task of the job runs the associators, it keeps the matched
#                             particles and vertices too.
#   truthGraphNanoTablesTask  the graph tables, the association tables and the associators they
#                             read (truthGraphNanoAssociatorsTask, unscheduled).
# Inputs: the graph and the hit index the DIGI step built (truthLogicalGraphProducer,
# truthLogicalGraphHitIndexProducer: Phase-2 eras, which include the enableTruth modifier, keep
# them in FEVTDEBUGHLT and RECOSIM, not in AODSIM), and for the associators the RECO-level
# collections they read. The associators fail without these inputs; the table producers write
# empty tables. Every module label ends with Table, so "keep nanoaodFlatTable_*Table_*_*" keeps
# the tables.
#
# Cost at 200 pileup (TTbar, D104, 8 threads): the associators 3.5 s/event (tracks 2.1 s,
# tracksters 1.2 s); the tables 0.06 s/event, 2.2 MB/event for the graph tables (49k of the
# graph's 159k particles) and 3.1 MB/event for the association tables (ZSTD 5).
# truthGraphTablesTask alone: 3 ms/event and 0.14 MB/event (3.1k particles).
#
# Simulated energy per particle and hit: TruthGraphHitFlatTableProducer
# (PhysicsTools/TruthInfo truthGraphHitFlatTable_cfi) on a hit table with a uint32 DetId column.

import FWCore.ParameterSet.Config as cms

# The cfi templates under private names, so process.load does not add them as modules.
from PhysicsTools.TruthInfo.truthGraphFlatTable_cfi import truthGraphFlatTable as _truthGraphFlatTable
from PhysicsTools.TruthInfo.truthBranchAssociationFlatTable_cfi import (
    truthBranchAssociationFlatTable as _truthBranchAssociationFlatTable,
)
from PhysicsTools.TruthInfo.truthBranchVertexAssociationFlatTable_cfi import (
    truthBranchVertexAssociationFlatTable as _truthBranchVertexAssociationFlatTable,
)
from SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociationLabels_cff import recoLabels
from SimGeneral.TruthGraphAssociatorProducers.truthGraphAssociators_cff import (
    truthBranchTargets,
    allTrackToTruthBranchAssociators,
    allVertexToTruthBranchAssociators,
    allSecondaryVertexToTruthBranchAssociators,
    truthBranchTracksterAssociators,
    truthBranchPFClusterEcalAssociators,
    truthBranchPFClusterHcalAssociators,
)

# Mantissa bits kept in the payload and score columns of the association tables (a float has
# 23): 12 bits are a relative precision of 2.4e-4, and every hit count up to 8192 is exact. The
# graph tables keep the precisions of their producer's defaults.
_payloadPrecision = 12

# Reco <-> truth particle maps: (associator module, domain of recoLabels). The reco -> truth
# direction uses the adaptive nominal working point (each reco object's best truth particle), the
# truth -> reco direction every reco object sharing energy, or hits for tracks, with a particle.
_particleDomains = (
    ("truthBranchTracksterAssociators", "tracksters"),
    ("allTrackToTruthBranchAssociators", "tracks"),
    ("truthBranchPFClusterEcalAssociators", "pfClustersEcal"),
    ("truthBranchPFClusterHcalAssociators", "pfClustersHcal"),
)
# Reco <-> truth vertex maps: the payload is the share of the vertex's track pT^2, and a reco
# vertex links to every truth vertex its tracks are matched to.
_vertexDomains = (
    ("allVertexToTruthBranchAssociators", "vertices"),
    ("allSecondaryVertexToTruthBranchAssociators", "secondaryVertices"),
)
_recoToTruth = "RecoToTruthAdaptiveNominal"
# The track maps count hits (BranchHitAssociator Metric::SharedHits), in both directions.
_trackPayload = dict(
    recoToTruthPayload=dict(
        name="nSharedHits",
        doc="Number of the track hits on cells of the truth branch",
        scoreDoc="1 - the fraction of the track hits on cells of the branch (lower is better)",
    ),
    truthToRecoPayload=dict(
        name="nSharedHits",
        doc="Number of the track hits on cells of the truth branch",
        scoreDoc="1 - the fraction of the branch tracker cells the track covers (lower is better)",
    ),
)


def _productKey(label):
    """The key the associators name their products with: label and instance joined by '_'."""
    return label.replace(":", "_")


def _tableKey(label):
    """The same in camel case, for module labels and table names (no '_' or ':' in either)."""
    parts = label.split(":")
    return parts[0] + "".join(p[0].upper() + p[1:] for p in parts[1:] if p)


def _upper(key):
    return key[0].upper() + key[1:]


def _maps(domains, direction):
    return [
        cms.InputTag(module, _productKey(label) + direction)
        for module, domain in domains
        for label in recoLabels(domain, "offline")
    ]


truthGraphTable = _truthGraphFlatTable.clone(
    keepMatchedBy=_maps(_particleDomains, _recoToTruth),
    keepMatchedByTruthToReco=_maps(_particleDomains, "TruthToReco"),
    keepVerticesMatchedBy=_maps(_vertexDomains, _recoToTruth),
    keepVerticesMatchedByTruthToReco=_maps(_vertexDomains, "TruthToReco"),
    doc="MC-truth graph particles: the signal interaction, the particles reco objects are matched to, and their "
    "ancestors",
)
truthGraphCaloBoundaryTable = truthGraphTable.clone(
    keepLevelMask=1 << 3,  # truth::LevelFlag::CaloBoundary
    doc="MC-truth graph particles: the signal interaction, every particle that entered the calorimeter, the "
    "particles reco objects are matched to, and their ancestors",
)

truthGraphAssociationTables = {}
for _module, _domain in _particleDomains + _vertexDomains:
    _isVertex = (_module, _domain) in _vertexDomains
    _template = _truthBranchVertexAssociationFlatTable if _isVertex else _truthBranchAssociationFlatTable
    _truth = "TruthVertices" if _isVertex else "TruthParticles"
    _targets = cms.InputTag("truthGraphTable", "selectedVertices" if _isVertex else "selected")
    _payload = _trackPayload if _domain == "tracks" else {}
    for _label in recoLabels(_domain, "offline"):
        _key = _tableKey(_label)
        truthGraphAssociationTables[f"truthGraph{_upper(_key)}ToTruthTable"] = _template.clone(
            src=cms.InputTag(_module, _productKey(_label) + _recoToTruth),
            targets=_targets,
            direction="recoToTruth",
            name=f"{_upper(_key)}ToTruth",
            doc=f"{_label} -> {_truth} ({_module}, {_recoToTruth[len('RecoToTruth'):]}): "
            + ("every truth vertex its tracks are matched to, best first" if _isVertex else "the best match"),
            precision=_payloadPrecision,
            **_payload,
        )
        truthGraphAssociationTables[f"truthGraphTruthTo{_upper(_key)}Table"] = _template.clone(
            src=cms.InputTag(_module, _productKey(_label) + "TruthToReco"),
            targets=_targets,
            direction="truthToReco",
            name=f"TruthTo{_upper(_key)}",
            doc=f"{_truth} -> {_label} ({_module}): every reco object sharing with the truth object, best first",
            precision=_payloadPrecision,
            **_payload,
        )
globals().update(truthGraphAssociationTables)

truthGraphNanoAssociatorsTask = cms.Task(
    truthBranchTargets,
    allTrackToTruthBranchAssociators,
    allVertexToTruthBranchAssociators,
    allSecondaryVertexToTruthBranchAssociators,
    truthBranchTracksterAssociators,
    truthBranchPFClusterEcalAssociators,
    truthBranchPFClusterHcalAssociators,
)
truthGraphTablesTask = cms.Task(truthGraphTable)
truthGraphNanoTablesTask = cms.Task(
    truthGraphTablesTask,
    *truthGraphAssociationTables.values(),
    truthGraphNanoAssociatorsTask,
)
