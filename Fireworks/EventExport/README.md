# Fireworks/EventExport

Exports events for the web event display as standard NanoAOD flat tables (`NanoAODOutputModule`).
Existing CMSSW table producers are reused; this package adds only what CMSSW lacks:

| Plugin | Output table | Notes |
|---|---|---|
| `{Ecal,HBHE,HF,HO,HGC}RecHitTableProducer` | `detId det subdet layer energy time x y z eta phi halfWidth` | One template for every calorimeter; geometry from the EventSetup (`hgcal::RecHitTools`, `CaloGeometry`) |
| `CaloClusterHitTableProducer` | `clusterIdx hitIdx fraction` (one row per link, ordered by cluster) | Any `reco::CaloCluster` collection; hit rows are resolved through the hit table's DetId column, so they always match it |
| `SoATableProducer<Collection>` (`interface/`), typedefs in `plugins/SoATableProducers.cc` | The SoA's columns, named in config | **Any `PortableHostCollection`.** Walks the layout descriptor: Eigen columns flattened (`name_0..`), blocks depth-first, leaves of a different size (e.g. `AssociationMap` offsets/content) in their own tables. A new or changed SoA = typedef + names list. On mismatch it prints the layout it found. Types that exist only in some areas (TICL v6) are guarded by `__has_include` |

SoAs are not persisted in RECO, so SoA tables run **inside the reco job**:
`addEventExport(process, presets=["hgcalSoA", ...])` (HGCal SoA rechits, merged SoA layer clusters, SoA tracksters + LC map).

Reused: `DPGAnalysis/HGCalNanoAOD` (tracksters with layer-cluster lists, TICL candidates, layer clusters,
tracks, sim objects, associations with scores) and NanoAOD's `Simple*FlatTableProducer` (vertices, PF candidates).

## Run

```bash
cmsRun Fireworks/EventExport/test/export_cfg.py inputFiles=file:step3.root outputFile=evd_nano.root \
    presets=calo,hgcal,reco geometry=ExtendedRun4D128 era=Phase2C26I13M9 globalTag=auto:phase2_realistic_T35
python3 Fireworks/EventExport/test/validate_export.py evd_nano.root
```
Inside a reco job: `cmsDriver.py ... --customise Fireworks/EventExport/customise.addEventExport`.
Presets: `calo` (all calorimeter hits), `hgcal` (HGCal hits, layer cluster → hit links, HGCalNanoAOD tables), `reco` (vertices, PF candidates).

## Checks (`validate_export.py`)

Finite positions and positive cell sizes; every cluster → hit link resolves to a hit row; layer-cluster energy
equals Σ fraction × hit energy; trackster → layer-cluster indices in range.
On the D128 CloseBy sample (10 events): all pass, LC energy agreement better than 1e-6.

**SoA, in a TICL v6 GPU step3 job (T4, ttbar, 5 events):** all checks pass:
- SoA layer clusters equal legacy in count and energy;
- SoA tracksters equal the legacy tracksters converted from them (count, raw energy, leading eigenvalue);
- trackster → LC offsets have n+1 entries and per-trackster LC sets equal legacy (order differs);
- principal axes are unit vectors to 2e-7 (zero only for 2-LC tracksters);
- SoA rechits are a subset of legacy (noise threshold: dropped hits ≤ 0.03 GeV) with identical positions.

The axis comparison also confirms that `DPGAnalysis/HGCalNanoAOD` `hgcalTracksters_cfi.py` swaps `eVector0_y` and `eVector0_z`.

## Cost (D128 CloseBy, 0PU, 4 threads)

`hgcalHitTable` 8.5 ms/event for about 20k hits (geometry lookups), `CaloClusterHitTableProducer` 0.3 ms,
NanoAOD output 8 ms. The whole job runs at 0.31 s/event, mostly the reused HGCalNanoAOD track propagation.

In the TICL v6 GPU job (ttbar, about 50k legacy and about 15k SoA HGCal hits per event):

| Module | ms/event |
|---|---|
| SoA tables, all five | 0.29 |
| Legacy `hgcalHitTable` | 18.5 |
| LC → hit links | 1.2 |
| NanoAOD output (ZLIB 4) | 57 |

The reconstruction itself takes about 8.4 s per event per stream, so the export is under 1%.
