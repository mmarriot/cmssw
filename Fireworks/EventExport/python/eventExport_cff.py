"""Event-display export: NanoAOD flat tables for the web viewer.

Reuses existing producers wherever they exist (HGCalNanoAOD for TICL objects, tracks and
associations; NanoAOD's simple table producers for vertices and PF candidates) and adds only
what CMSSW lacks: calorimeter rechits with geometry, and cluster -> hit links.
"""
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, P3Vars


def _hitTable(plugin, name, doc, src):
    return cms.EDProducer(plugin, name=cms.string(name), doc=cms.string(doc),
                          src=cms.VInputTag(*src), skipMissing=cms.bool(True))


# --- calorimeter rechits with geometry (one table per hit type)
hgcalHitTable = _hitTable("HGCRecHitTableProducer", "HGCalHits", "HGCal rechits (CE-E, CE-H silicon, CE-H scintillator)",
                          ["HGCalRecHit:HGCEERecHits", "HGCalRecHit:HGCHEFRecHits", "HGCalRecHit:HGCHEBRecHits"])
ecalHitTable = _hitTable("EcalRecHitTableProducer", "EcalHits", "ECAL barrel rechits", ["ecalRecHit:EcalRecHitsEB"])
hbheHitTable = _hitTable("HBHERecHitTableProducer", "HBHEHits", "HCAL barrel/endcap rechits", ["hbhereco"])
hfHitTable = _hitTable("HFRecHitTableProducer", "HFHits", "HCAL forward rechits", ["hfreco"])
hoHitTable = _hitTable("HORecHitTableProducer", "HOHits", "HCAL outer rechits", ["horeco"])

# --- cluster -> hit links
hgcalLayerClusterHitTable = cms.EDProducer("CaloClusterHitTableProducer",
    name=cms.string("HGCalLayerClusterHits"),
    doc=cms.string("HGCal layer cluster -> HGCalHits links"),
    clusters=cms.InputTag("hgcalMergeLayerClusters"),
    hitTable=cms.InputTag("hgcalHitTable"),
    hitIdColumn=cms.string("detId"),
)

# --- generic reco objects (NanoAOD simple producers)
vertexTable = cms.EDProducer("SimpleVertexFlatTableProducer",
    src=cms.InputTag("offlinePrimaryVertices"), name=cms.string("Vertices"), doc=cms.string("offline primary vertices"),
    cut=cms.string(""), singleton=cms.bool(False), extension=cms.bool(False),
    variables=cms.PSet(
        x=Var("x", float, doc="x [cm]"), y=Var("y", float, doc="y [cm]"), z=Var("z", float, doc="z [cm]"),
        ndof=Var("ndof", float, doc="number of degrees of freedom"),
        chi2=Var("normalizedChi2", float, doc="chi2/ndof"),
        nTracks=Var("tracksSize", "uint16", doc="number of tracks"),
        isFake=Var("isFake", bool, doc="fake vertex"),
    ),
)
pfCandidateTable = cms.EDProducer("SimplePFCandidateFlatTableProducer",
    src=cms.InputTag("particleFlow"), name=cms.string("PFCandidates"), doc=cms.string("particle-flow candidates"),
    cut=cms.string(""), singleton=cms.bool(False), extension=cms.bool(False),
    variables=cms.PSet(P3Vars,
        energy=Var("energy", float, doc="energy [GeV]"),
        pdgId=Var("pdgId", int, doc="PDG id"),
        charge=Var("charge", int, doc="electric charge"),
        vx=Var("vx", float, doc="vertex x [cm]"), vy=Var("vy", float, doc="vertex y [cm]"), vz=Var("vz", float, doc="vertex z [cm]"),
    ),
)

# --- SoA products (in-job only: SoAs are not persisted). Column lists mirror the layouts, one entry per
# leaf layout in declaration order; a mismatch makes the producer print the layout it found.
def _soaTable(plugin, src, name, doc, *blocks):
    return cms.EDProducer(plugin, src=cms.InputTag(src), name=cms.string(name), doc=cms.string(doc),
                          blocks=cms.VPSet(*[cms.PSet(columns=cms.vstring(*cols), table=cms.string(table))
                                             for table, cols in blocks]))

# dim1/dim2 are x/y [cm] for silicon and eta/phi for scintillator (set by the SoA rechit producer).
_hgcalSoARecHitColumns = lambda dim1, dim2: (dim1, dim2, "z:position z [cm]", "layer", "energy:energy [GeV]",
                                             "mipEnergy:energy in MIP where calibrated (DAQ chain); equals energy in the simulation chain", "sigmaNoise:noise, same units as mipEnergy", "recHitIndex", "flags",
                                             "detId:raw DetId", "time:time [ns]", "timeError:time error [ns]")
_silicon, _scintillator = ("x:position x [cm]", "y:position y [cm]"), ("eta:pseudorapidity", "phi:azimuth [rad]")
hgcalSoARecHitTables = {sub: _soaTable("HGCalRecHitSoATableProducer", f"hgcalSoARecHits{sub}", f"HGCalSoARecHits{sub}",
                                       f"HGCal {sub} rechits (SoA)", ("", _hgcalSoARecHitColumns(*dims)))
                        for sub, dims in (("EE", _silicon), ("HSi", _silicon), ("HSci", _scintillator))}
hgcalSoARecHitsEETable, hgcalSoARecHitsHSiTable, hgcalSoARecHitsHSciTable = hgcalSoARecHitTables.values()

hgcalSoALayerClusterTable = _soaTable("CaloClusterSoATableProducer", "hgcalSoAMergeLayerClusters",
    "HGCalSoALayerClusters", "HGCal layer clusters (SoA, merged)",
    ("", ("x:position x [cm]", "y:position y [cm]", "z:position z [cm]", "layer", "cells:number of cells")),
    ("", ("energy:energy [GeV]", "correctedEnergy:[GeV]", "correctedEnergyUncertainty:[GeV]")),
    ("", ("", "algoID", "seedId:seed raw DetId", "flags")),  # caloID skipped: polymorphic type, not exportable
    ("", ("time:[ns]", "timeError:[ns]")))

ticlSoATrackstersTable = _soaTable("TracksterSoATableProducer", "ticlSoATrackstersPID",
    "TICLSoATracksters", "TICL tracksters (SoA)",
    ("", ("rawEnergy:[GeV]", "rawEmEnergy:[GeV]", "regressedEnergy:[GeV]")),
    ("", ("x:barycenter x [cm]", "y:barycenter y [cm]", "z:barycenter z [cm]")),
    ("", ("time:[ns]", "timeError:[ns]")),
    ("", ("eigenvalues:PCA eigenvalues", "eigenvectors:PCA eigenvectors (column-major)", "sigmas", "sigmasPCA")),
    ("", ("idProbabilities:particle-ID probabilities",)),
    ("", ("maxEnergyFraction", "core2EnergyFraction", "core5EnergyFraction", "rmsRadius:[cm]", "minAbsZ:[cm]",
          "lengthZ:[cm]")),
    ("TICLSoATrackstersLCOffsets", ("offset:start of each trackster's layer clusters in TICLSoATrackstersLCs",)),
    ("TICLSoATrackstersLCs", ("layerCluster:index into HGCalSoALayerClusters",)))

# --- presets: compose what you need, e.g. addEventExport(process, presets=["calo", "hgcal", "reco"])
presets = {
    "calo": cms.Task(hgcalHitTable, ecalHitTable, hbheHitTable, hfHitTable, hoHitTable),
    "hgcal": cms.Task(hgcalHitTable, hgcalLayerClusterHitTable),  # + HGCalNanoAOD tables, see customise.py
    "reco": cms.Task(vertexTable, pfCandidateTable),
    "hgcalSoA": cms.Task(*hgcalSoARecHitTables.values(), hgcalSoALayerClusterTable, ticlSoATrackstersTable),
}
