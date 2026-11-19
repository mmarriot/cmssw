# Checks the output of testTrackHitsTable_cfg.py against the expected rows written by TrackHitsTestInputProducer.
import sys
import numpy as np
import uproot

events = uproot.open(sys.argv[1])["Events"]
assert events.num_entries == 1
row = {name: array[0] for name, array in events.arrays(library="np").items()}


def column(name):
    return np.asarray(row[name])


# track 0: a pixel hit, a missing hit (left out), an RPC hit; track 1: no hits; track 2: an RPC hit
assert column("Track_nTrackHits").tolist() == [2, 0, 1], column("Track_nTrackHits")
assert row["nTrackHits"] == 3
assert column("TrackHits_detId").tolist() == column("Expected_detId").tolist()
assert column("TrackHits_layer").tolist() == column("Expected_layer").tolist()
# the pixel hit within 100 um of the cluster barycentre (the CPE corrections, e.g. the Lorentz drift, are smaller),
# the RPC hits at their stored position
for axis in "xyz":
    difference = np.abs(column("TrackHits_" + axis) - column("Expected_" + axis))
    assert difference[0] < 0.01, (axis, difference)
    assert np.all(difference[1:] < 1e-4), (axis, difference)
# with cut = "pt > 15" and maxLen = 1: track 1 only, which has no hits
assert column("CutTrack_nCutTrackHits").tolist() == [0]
assert row["nCutTrackHits"] == 0
print("checkTrackHitsTable.py: OK")
