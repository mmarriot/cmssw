# TruthInfo

> **Status: open for contributions.** The data model and the levels are covered by
> unit tests; a change to either should come with one. The package targets
> **Phase-2 (Run 4) only** (no Phase-1/Run-2 support).
>
> **Original author and maintainer:** Felice Pantaleo (CERN),
> <felice.pantaleo@cern.ch>.

A prototype **MC-truth graph** for CMS: a single, navigable, physics-oriented
abstraction of the generator + simulation truth history of an event, with
calorimeter and tracker hit indices layered on top. It replaces the need to
cross-navigate the many low-level truth collections (HepMC, GenParticles,
SimTracks/SimVertices, TrackingParticles, SimClusters, CaloParticles, SimHits,
RecHits) by hand.

## Documentation

The **authoritative, maintained documentation** lives on the project website:

> **http://cms-truth.docs.cern.ch/**

It covers the data model, how to enable and use the graph (navigation API, the
`truth::Branch` view, the hit index, matching reco objects), the physics findings,
validation, pileup, and the roadmap. The MkDocs sources are in the companion
`cms-truth-docs` repository. Start with the **"How to use the graph"** page.

## The three layers

1. **`TruthGraph`** (raw): a compact CSR graph built directly from HepMC +
   `SimTrack`/`SimVertex` by `TruthGraphProducer`.
2. **`truth::Graph`** (logical): a user-facing bipartite Particle<->Vertex graph
   built by `TruthLogicalGraphProducer`; GEN and SIM are merged where robustly
   associated, with navigation (`parents()`, `descendants()`,
   `firstCommonAncestor()`, `hasAncestorPdgId()`, ...) and the `truth::Branch`
   subgraph view + `BranchSelector` selection.
3. **`truth::LogicalGraphHitIndex`**: per-particle direct vs aggregated subgraph
   calorimeter and tracker hits, built by `LogicalGraphHitIndexProducer` (with the
   DetId->RecHit map from `DetIdToRecHitMapProducer`).

Producer chain (order matters): `truthGraphProducer` -> `truthLogicalGraphProducer`
-> `detIdToRecHitMapProducer` -> `truthLogicalGraphHitIndexProducer`. In a release
job these run behind the `enableTruth` process modifier (the
`truthGraphPrevalidation` sequence in `Validation/Configuration`).

## Package layout

- `interface/` + `src/`: the data formats and algorithms: `TruthGraph`,
  `truth::Graph`, `Branch`, `BranchSelector`, `LogicalGraphHitIndex`,
  `BranchHitAssociator`, the `truth::recoHits` adapters (`RecoHitAdapters.h`), and
  `TruthLogicalGraphPostProcessor` (merge/collapse/filter; covered by the cppunit).
- `plugins/`: the producers above, the DOT dumpers, the flat-table producers and
  the pileup `TruthGraphAccumulator`.
- `python/`: the DIGI configuration (`truthGraphMixedDigi_cff`), the reduced variant
  (`mixedTruthGraphCustomize.customiseTruthReduced`), the RECO content customise
  (`customiseTruthMixedReco`) and the event content.
- `test/`: cppunit unit tests and standalone `cmsRun` drivers (graph dumps,
  topology checks).

The DQM validators, their association-map producers, `truthGraphValidation_cff`,
`truthGraphDQMHarvester_cff` and the plot scripts live in `Validation/TruthInfo`.
The reco-to-truth associators with working points live in
`SimGeneral/TruthGraphAssociatorProducers`.

## NanoAOD tables of the graph

`TruthGraphFlatTableProducer` and `TruthGraphHitFlatTableProducer`
(`plugins/TruthGraphNanoTableProducers.cc`) write the graph as NanoAOD flat tables.
`SimGeneral/TruthGraphAssociatorProducers/python/truthGraphNanoTables_cff.py` configures
them: `truthGraphTablesTask` has the graph tables alone, `truthGraphNanoTablesTask` adds the
reco <-> truth association tables and their associators (every module label ends in
`Table`, so `keep nanoaodFlatTable_*Table_*_*` keeps the tables). The tables need the graph
and the hit index of the DIGI step (`truthLogicalGraphProducer`,
`truthLogicalGraphHitIndexProducer`), which the Phase-2 eras (they include the
`enableTruth` modifier) keep in the FEVTDEBUGHLT and RECOSIM event content, not in AODSIM.
The associators also need the reco collections they read; they run unscheduled and fail
on a missing input. The table producers alone write empty tables, with all their columns,
when an input is missing (`skipNonExistingSrc`), so `truthGraphTablesTask` runs on any
input.

- `TruthParticles` (`TruthGraphFlatTableProducer`, module `truthGraphTable`): a
  selection of the particles, one row each in ascending graph id: the signal interaction
  (`keepSignal`), every particle with one of the `keepLevelMask` `truth::LevelFlag` bits,
  pile-up included, every particle a reco object is matched to (`keepMatchedBy`: entries
  of reco -> truth maps; `keepMatchedByTruthToReco`: rows with entries of truth -> reco
  maps), and the ancestors of all of them (`keepAncestors`). Columns: `graphId`, PDG id,
  charge, the momentum as `pt eta phi energy` (the columns of `TruthCB`), generator
  `status` and `statusFlags`, `levelFlags`, `isSignal`, `bunchCrossing`, `eventIndex`,
  `hasGen`, `hasSim`, `role`, the rows of the production vertex and of the endpoint in
  `TruthVertices` (`prodVertexIdx`, `decayVertexIdx`, -1 if none; the production position
  is that vertex row's), and the tracker-calorimeter boundary crossing of the SimTrack
  (`truth::Checkpoint` 0): `hasBoundary` and `boundaryX Y Z T Pt Eta Phi E`, zeros without
  one, as `TruthCB` (`TruthParticleFlatTableProducer`) writes them. The momentum and the
  positions have one representation each (no `px py pz mass`, no copy of the vertex
  position on the particle). The endpoint is the latest in time of the particle's decay or
  interaction vertices (the higher graph id on a tie), i.e. where its last saved
  secondaries were made, not necessarily where it stopped: a Geant4 track can end at
  several, and the children of all of them are in `TruthParticleChildren`. Floats keep 14 mantissa bits for momenta and angles and 16 for
  positions and times (`precision`, `positionPrecision`).
- `TruthVertices` (instance `vertices`): the production vertex and the endpoint of each
  exported particle, and every truth vertex a vertex map matches (`keepVerticesMatchedBy`,
  `keepVerticesMatchedByTruthToReco`), in ascending graph id, with `graphId`, `x y z t`,
  `reason`, `role`, `isSignal`.
- `TruthParticleChildren` (instance `children`): `parentIdx` -> `childIdx` links
  between exported particles through every end vertex of the parent, ordered by parent.
- The graph ids of the particle and vertex rows, in row order (`selected`,
  `selectedVertices`, `std::vector<unsigned int>`).
- Reco <-> truth associations, `TruthBranchAssociationFlatTableProducer` (particle
  maps) and `TruthBranchVertexAssociationFlatTableProducer` (vertex maps), in both
  directions for the collections `recoLabels` names in the domains `tracks`,
  `tracksters`, `pfClustersEcal`, `pfClustersHcal`, `vertices` and `secondaryVertices`
  (not `candidates` and `layerClusters`): `<Reco>ToTruth` (adaptive nominal working
  point, rows parallel to the reco collection) and `TruthTo<Reco>` (rows parallel to
  `TruthParticles` or `TruthVertices`), each with `n<A>Links`/`o<A>Links` columns and a
  `<A>Links` table of `index` (row in the other table), `graphId`, the payload and
  `score`. The payload is `sharedEnergy` (reco -> truth) or `sharedEnergyFraction`
  (truth -> reco) for the calorimeter maps, `nSharedHits` for the track maps (a hit
  count in both directions), `fraction` for the vertex maps. The tables read `selected`
  or `selectedVertices` as their `targets`, so `index` is a row of `TruthParticles` or
  `TruthVertices`. Every map entry has its particle or vertex exported, so no link is
  left out.

The default selection is the signal, the matched particles and their ancestors: at 200
pile-up (TTbar, D104) 49k of the graph's 159k particles. `truthGraphCaloBoundaryTable`,
the opt-in for a hit-level truth, adds every particle that entered the calorimeter
(`keepLevelMask = 1 << 3`): 117k particles, carrying 99.8% of the simulated energy in
the HGCAL rechits. It replaces `truthGraphTable`, so the other tables keep reading it:
`process.truthGraphTable = process.truthGraphCaloBoundaryTable.clone()`. The associators
cost 3.5 s/event at 200 pile-up (tracks 2.1 s, tracksters 1.2 s), the tables 0.06 s;
without associations (`truthGraphTablesTask`) the tables hold the signal and its
ancestors (3.1k particles) and take 3 ms/event.

`TruthGraphHitFlatTableProducer` (`truthGraphHitFlatTable_cfi`) adds the simulated
energy of the exported particles in the hits of any hit table with a uint32 DetId
column (`hitTable`, `hitIdColumn`): one row per (particle, hit) with `particleIdx`,
`hitIdx` and `simEnergy` [GeV], the energy the particle and its descendants left in
the hit's cell, from the `Calo` or `Muon` channel of the hit index. Only exported
particles with one of the `levelMask` bits get rows (default `CaloBoundary`, an
antichain, so no energy is counted twice; use it with `truthGraphCaloBoundaryTable`).
On the `Muon` channel the hit index is keyed by the simulation's unit (DT wire, CSC
layer, RPC roll, GEM eta partition): every rechit of a unit gets the unit's full
energy, so summing over rechits that share a DetId counts it more than once. For
example, on the HGCAL rechit table of this package, which reads the TICL geometry
(`RecoTICL.Geometry.TICLGeom_cff`, made by an alpaka ES producer, so the job also
needs `Configuration.StandardSequences.Accelerators_cff`):

```python
from PhysicsTools.TruthInfo.recHitTable_cfi import recHitTable
from PhysicsTools.TruthInfo.truthGraphHitFlatTable_cfi import truthGraphHitFlatTable
truthGraphHGCalHitTable = truthGraphHitFlatTable.clone(
    hitTable="recHitTable", hitIdColumn="rechit_ID", name="TruthHGCalRecHits")
```

## Build, check, test

All commands assume the CMSSW environment (`cd $CMSSW_BASE/src && cmsenv`):

```bash
scram b -j 8                          # build
scram b code-format code-checks -j 8  # clang-format + clang-tidy (must pass)
scram b runtests                      # cppunit unit tests
```

## Run standalone on a step3.root

```bash
# Dump per-event DOT graphs (options: -n, -m/--merge, -c/--collapse, -o, -t)
cmsRun test/dumpTruthGraphsFromGENSIMRECO_cfg.py path/to/step3.root -n 5

# Branch DQM validators (calo / tracking), from Validation/TruthInfo
cmsRun ../../Validation/TruthInfo/test/validateBranchDQM_cfg.py          path/to/step3.root -n 5
cmsRun ../../Validation/TruthInfo/test/validateBranchTrackingDQM_cfg.py  path/to/step3.root -n 5
```

See the website for the full configuration reference, the navigation API with
examples, and the validation results.
