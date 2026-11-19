#!/bin/bash
# TrackHitsTableProducer on tracks made by TrackHitsTestInputProducer (needs the conditions database)
function die { echo "$1: status $2"; exit $2; }

cmsRun ${SCRAM_TEST_PATH}/testTrackHitsTable_cfg.py || die "cmsRun testTrackHitsTable_cfg.py" $?
python3 ${SCRAM_TEST_PATH}/checkTrackHitsTable.py testTrackHitsTable.root || die "checkTrackHitsTable.py" $?

# tracks without their extras and hits, as in AOD: a clear failure
cmsRun ${SCRAM_TEST_PATH}/testTrackHitsTable_cfg.py aod > testTrackHitsTableAOD.log 2>&1 && die "no failure on AOD" 1
grep -q "this table needs RECO or FEVT" testTrackHitsTableAOD.log || die "no message on AOD" 1
echo "testTrackHitsTable.sh: OK"
