#!/usr/bin/env python3
"""The device modules read the HGCAL components in the order of the merged layer clusters of the host, offline and at
HLT: a CLUEstering trackster indexes its layer clusters across the components, so the host merge reads the same SoA
collections in the same order."""

import sys

from RecoTICL.Configuration import presets
from RecoTICL.Configuration.target import TARGETS


def _hlt_module(label):
    # the cfi file of an HLT module is not always named after its label
    import glob, os, re
    import HLTrigger.Configuration.HLT_75e33.modules as modules
    for path in sorted(glob.glob(os.path.join(modules.__path__[0], "*_cfi.py"))):
        with open(path) as f:
            if re.search(r"^%s\s*=" % label, f.read(), re.M):
                name = "HLTrigger.Configuration.HLT_75e33.modules." + os.path.basename(path)[:-3]
                return getattr(__import__(name, fromlist=[label]), label)
    raise LookupError("no HLT_75e33 module cfi defines %s" % label)


def _offline_module(label):
    import RecoTICL.LayerClustering.hgcalLayerClusters_cff as cff
    return getattr(cff, label)


def _label(tag):
    return tag.getModuleLabel() if hasattr(tag, "getModuleLabel") else str(tag).split(":")[0]


def check(name, host_merge, reader, labels):
    host = [_label(tag) for tag in host_merge.layerClustersSoA]
    device = [_label(tag) for tag in labels]
    if host != device:
        print("FAIL [%s]: the host merges %s; %s reads %s" % (name, host, reader, device))
        return 1
    print("OK [%s]: %s reads %s in the order of the host merge" % (name, reader, device))
    return 0


def readers(assembled, target):
    pid = assembled.modules[target.soa_pid_tracksters_label]
    return [(target.cluestering_label, assembled.modules[target.cluestering_label].hgcalLayerClustersSoA),
            (target.soa_tracksters_label, assembled.modules[target.soa_tracksters_label].hgcalLayerClustersSoA),
            (target.soa_pid_tracksters_label, pid.layerClusters)]


def main():
    rc = 0
    for name, host_merge in (("offline", _offline_module("hgcalMergeLayerClusters")),
                             ("HLT", _hlt_module(TARGETS["HLT"].merged_lc))):
        config = presets.v6(target=name)
        config.iterations[0].device_pid = dict(model="CNN")
        for reader, labels in readers(config.assemble(), TARGETS[name]):
            rc |= check(name, host_merge, reader, labels)
    return rc


if __name__ == "__main__":
    sys.exit(main())
