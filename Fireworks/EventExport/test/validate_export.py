#!/usr/bin/env python3
"""Consistency checks on an event-display NanoAOD file. Exit code 1 if any check fails.

python3 validate_export.py evd_nano.root
"""
import sys

import numpy as np
import uproot

failures = []


def check(ok, what):
    print(f"[{'ok' if ok else 'FAIL'}] {what}")
    if not ok:
        failures.append(what)


def main(path):
    events = uproot.open(path)["Events"]
    names = set(events.keys())
    arrays = lambda *b: events.arrays(list(b), library="np")

    for table in ("HGCalHits", "EcalHits", "HBHEHits", "HFHits", "HOHits"):
        if f"n{table}" not in names:
            continue
        a = arrays(f"n{table}", f"{table}_x", f"{table}_y", f"{table}_z", f"{table}_halfWidth", f"{table}_energy")
        n = a[f"n{table}"]
        xyz = np.concatenate([np.concatenate(a[f"{table}_{c}"]) for c in "xyz"]) if n.sum() else np.zeros(0)
        hw = np.concatenate(a[f"{table}_halfWidth"]) if n.sum() else np.zeros(0)
        check(np.isfinite(xyz).all(), f"{table}: {n.sum()} hits over {len(n)} events, finite positions")
        if n.sum():
            print(f"       halfWidth [cm] median {np.median(hw):.2f}, range {hw.min():.2f}-{hw.max():.2f}")
            check((hw > 0).mean() > 0.99, f"{table}: halfWidth > 0 for >99% of hits")

    if "nHGCalLayerClusterHits" in names:
        a = arrays("HGCalLayerClusterHits_clusterIdx", "HGCalLayerClusterHits_hitIdx", "HGCalLayerClusterHits_fraction",
                   "HGCalHits_energy", "HGCalLayerClusters_energy", "nHGCalHits")
        rel, ordered, unmatched, links = [], True, 0, 0
        for ev in range(len(a["nHGCalHits"])):
            cl, hit, frac = (a[f"HGCalLayerClusterHits_{c}"][ev] for c in ("clusterIdx", "hitIdx", "fraction"))
            ordered &= bool(np.all(np.diff(cl) >= 0))
            unmatched += int(np.sum((hit < 0) | (hit >= a["nHGCalHits"][ev])))
            links += len(hit)
            lc_e = a["HGCalLayerClusters_energy"][ev]
            summed = np.bincount(cl, weights=frac * a["HGCalHits_energy"][ev][hit], minlength=len(lc_e))
            rel.append(np.abs(summed - lc_e) / np.maximum(lc_e, 1e-6))
        check(ordered, "cluster -> hit links ordered by cluster")
        check(unmatched == 0, f"all {links} cluster -> hit links resolve to a hit row ({unmatched} unmatched)")
        rel = np.concatenate(rel)
        print(f"       LC energy vs sum(fraction*hit energy): median rel. diff {np.median(rel):.1e}, max {rel.max():.1e}")
        check(np.quantile(rel, 0.99) < 1e-3, "layer-cluster energy = sum of linked hit energies (99% within 0.1%)")

    for ts in ("ticlTrackstersCLUE3DHigh", "ticlTracksterLinks"):
        if f"n{ts}" in names:
            a = arrays(f"{ts}vertices_vertices", "nHGCalLayerClusters")
            ok = all((v < n).all() for v, n in zip(a[f"{ts}vertices_vertices"], a["nHGCalLayerClusters"]))
            check(ok, f"{ts}: trackster -> layer-cluster indices within range")

    if "nTICLSoATracksters" in names:
        a = arrays("nTICLSoATracksters", "TICLSoATrackstersLCOffsets_offset", "TICLSoATrackstersLCs_layerCluster",
                   "nHGCalSoALayerClusters", *[f"TICLSoATracksters_eigenvectors_{k}" for k in range(9)])
        csr_ok, unit_ok = True, True
        for ev, nT in enumerate(a["nTICLSoATracksters"]):
            off = a["TICLSoATrackstersLCOffsets_offset"][ev].astype(np.int64)
            lcs = a["TICLSoATrackstersLCs_layerCluster"][ev]
            csr_ok &= bool(len(off) == nT + 1 and off[0] == 0 and off[-1] == len(lcs) and (np.diff(off) >= 0).all()
                           and (lcs < a["nHGCalSoALayerClusters"][ev]).all())
            axis = np.stack([a[f"TICLSoATracksters_eigenvectors_{k}"][ev] for k in range(3)], axis=1)  # column 0
            norm = np.linalg.norm(axis, axis=1)
            unit_ok &= bool(np.all((norm < 1e-6) | (np.abs(norm - 1) < 1e-5)))  # zero when no PCA (<3 LCs)
        check(csr_ok, "TICLSoATracksters: LC offsets have n+1 entries, monotonic, indices in range")
        check(unit_ok, "TICLSoATracksters: principal axes are unit vectors (or zero)")
        if "nticlTrackstersCLUEstering" in names:  # legacy tracksters converted from the same SoA
            b = arrays("TICLSoATracksters_rawEnergy", "ticlTrackstersCLUEstering_raw_energy")
            same = all(np.allclose(np.sort(s), np.sort(l), rtol=1e-5)
                       for s, l in zip(b["TICLSoATracksters_rawEnergy"], b["ticlTrackstersCLUEstering_raw_energy"]))
            check(same, "TICLSoATracksters: raw energies equal the legacy tracksters converted from the SoA")
    if "nHGCalSoALayerClusters" in names and "nHGCalLayerClusters" in names:
        b = arrays("HGCalSoALayerClusters_energy", "HGCalLayerClusters_energy")
        same = all(np.isclose(s.sum(), l.sum(), rtol=1e-5)
                   for s, l in zip(b["HGCalSoALayerClusters_energy"], b["HGCalLayerClusters_energy"]))
        check(same, "HGCalSoALayerClusters: energy sum equals the legacy layer clusters")

    print("FAILED: " + "; ".join(failures) if failures else "all checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
