import FWCore.ParameterSet.Config as cms

from ..modules.hltFilteredLayerClustersCLUE3DHigh_cfi import *
from ..modules.hltTiclSeedingGlobal_cfi import *
from ..modules.hltTiclTrackstersCLUE3DHigh_cfi import *

HLTTiclTrackstersCLUE3DHighStepSequence = cms.Sequence(hltFilteredLayerClustersCLUE3DHigh+hltTiclSeedingGlobal+hltTiclTrackstersCLUE3DHigh)

# with ticl_dev the CLUEstering assignment CLUE3DHigh reads is run on device.
from Configuration.ProcessModifiers.ticl_dev import ticl_dev
from ..modules.hltTiclTrackstersCLUEsteringAssignment_cfi import *
ticl_dev.toModify(HLTTiclTrackstersCLUE3DHighStepSequence,
                  func = lambda s : s.associate(cms.Task(hltTiclTrackstersCLUEsteringAssignment)))
