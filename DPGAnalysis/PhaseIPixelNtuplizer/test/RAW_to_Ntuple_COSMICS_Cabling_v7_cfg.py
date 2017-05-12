import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PhaseINtuplizer', eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

inputFilesFromReco291877 = cms.untracked.vstring( # with BPix
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/22673ADA-9623-E711-8CE5-02163E0139AC.root",
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/127CC029-9423-E711-978E-02163E0141EA.root",
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/34DACE4A-9923-E711-BF38-02163E0144E5.root",
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/946728E0-AD23-E711-B999-02163E011B43.root",
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/BE55B4D7-9423-E711-A2D6-02163E019C47.root",
"file:/data/veszpv/Cosmics2017/GlobalApr17/RAW_291877/C07E24CD-9523-E711-BE66-02163E011808.root",
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = inputFilesFromReco291877,
    secondaryFileNames = cms.untracked.vstring())

process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('PhaseINtuplizer run on reco'),
    name = cms.untracked.string('PhaseINtuplizer'),
    version = cms.untracked.string('$alpha$'))

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_dataRun2_Express_v2', '')

# ----------------------------------------------------------

process.CablingMapDBReader = cms.ESSource("PoolDBESSource",
  BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0),
    authenticationPath = cms.untracked.string('')
  ),
  ##connect = cms.string('sqlite_file:SiPixelFedCablingMap_phase1_v8.db'), #local DB
  connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'), 
 #SiPixelCabling_Ph0andPilotBlade.db --> PB + BPix + FPix
 #SiPixelCabling_PilotBlade.db --> just PB, no BPix, no FPix
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('SiPixelFedCablingMapRcd'),
      #label = cms.untracked.string('pilotBlade'), 
      #tag = cms.string('SiPixelFedCablingMap_data') #tagname for the local DB
      tag = cms.string('SiPixelFedCablingMap_phase1_v7_May2'),  
      ##tag = cms.string('SiPixelFedCablingMap_phase1_v6.3_May2'),  
    )
  )
)
process.es_prefer_CablingReader = cms.ESPrefer("PoolDBESSource","CablingMapDBReader")

# ----------------------------------------------------------

process.NewAlignment = cms.ESSource("PoolDBESSource",
        BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
        DBParameters = cms.PSet(
            messageLevel = cms.untracked.int32(0),
            authenticationPath = cms.untracked.string('')
        ),
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        ##connect = cms.string('frontier://FrontierPrep/CMS_CONDITIONS'),
        toGet = cms.VPSet(
           cms.PSet(
           record = cms.string('TrackerAlignmentRcd'),
           tag = cms.string('TrackerAlignment_CRUZET17_v1'),
           )
        )
##         connect = cms.string('sqlite_file:myAlignmentFile.db'),
##         toGet = cms.VPSet(
##           cms.PSet(
##               record = cms.string("TrackerAlignmentRcd"),
##               tag = cms.string("Alignments")
##           ),
##           cms.PSet(
##               record = cms.string("TrackerAlignmentErrorExtendedRcd"),
##               tag = cms.string("AlignmentErrorsExtended")
##           ),
##           cms.PSet(
##               record = cms.string("TrackerSurfaceDeformationRcd"),
##               tag = cms.string("AlignmentSurfaceDeformations")
##           ),
##          )
)
process.es_prefer_Alignment = cms.ESPrefer("PoolDBESSource","NewAlignment")

# ----------------------------------------------------------

import HLTrigger.HLTfilters.hltHighLevel_cfi

process.HLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = ['HLT_L1SingleMuCosmics_v*',
                'HLT_L1SingleMuOpen_v*',
                'HLT_L1SingleMuOpen_DT_v*',
                'HLT_L1SingleMu3_v*'],
    throw = False #dont throw except on unknown path names
)

#---------------------------
#  Track Refitter
#---------------------------
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
process.offlineBeamSpot = offlineBeamSpot 

from RecoTracker.TrackProducer.TrackRefitters_cff import *
process.Refitter = RecoTracker.TrackProducer.TrackRefitterP5_cfi.TrackRefitterP5.clone()
process.Refitter.src = 'ctfWithMaterialTracksP5'
process.Refitter.NavigationSchool = ""

process.PhaseIPixelNtuplizerPlugin = cms.EDAnalyzer('PhaseIPixelNtuplizer',
    trajectoryInput = cms.InputTag('Refitter'),
    Cosmics = cms.int32(1),
    fileName = cms.untracked.string("Ntuple.root"),
)

process.hltFilter_step = cms.Sequence(process.HLTFilter)
process.raw2digi_step = cms.Sequence(process.RawToDigi)
process.L1Reco_step = cms.Sequence(process.L1Reco)
process.reconstruction_step = cms.Sequence(process.reconstructionCosmics)
process.PhaseIPixelNtuplizer_step = cms.Path(process.hltFilter_step*process.raw2digi_step*process.L1Reco_step*process.reconstruction_step* process.Refitter * process.PhaseIPixelNtuplizerPlugin)

# do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)
