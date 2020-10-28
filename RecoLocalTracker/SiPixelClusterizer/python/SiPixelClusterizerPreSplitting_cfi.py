
import FWCore.ParameterSet.Config as cms

#
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *
from RecoLocalTracker.SiPixelClusterizer.SiPixelClusterizer_cfi import siPixelClusters as _siPixelClusters
siPixelClustersPreSplitting = _siPixelClusters.clone()

from Configuration.Eras.Modifier_run3_common_cff import run3_common
run3_common.toModify(
    siPixelClustersPreSplitting,
    src = cms.InputTag('siPixelDigisMorphed')
)
