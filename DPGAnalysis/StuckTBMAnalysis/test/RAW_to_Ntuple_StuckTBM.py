import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('StuckTBMAnalyzer', eras.Run2_2017)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
      "file:/data/veszpv/ZeroBiasRaw2017/300806/02C685A4-567D-E711-9C28-02163E0127B8.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/0A9A5FC6-567D-E711-893A-02163E0134B5.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/109E0792-567D-E711-A08B-02163E0143A3.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/20CD9E24-577D-E711-A0BE-02163E019E54.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/2EF366C4-567D-E711-89FE-02163E01A5C6.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/12D6531A-577D-E711-822D-02163E014282.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/34999391-567D-E711-A1E7-02163E0133E1.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/3AA6B8A4-567D-E711-B975-02163E0144B9.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/3CFCFEB0-567D-E711-9B65-02163E014218.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/46C64CB7-567D-E711-8615-02163E019DD0.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/484E43AF-567D-E711-A73F-02163E012514.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/4A3783A2-567D-E711-8380-02163E019E02.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/4C62C592-567D-E711-8FBA-02163E01390F.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/4EC05BA1-567D-E711-B802-02163E0138DB.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/545D7DD7-567D-E711-BFFA-02163E01374B.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/5C65789F-567D-E711-8BA4-02163E0138DB.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/62028BC1-567D-E711-BA46-02163E011EBA.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/6AAF1DA7-567D-E711-B507-02163E019D12.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/70A203BE-567D-E711-B482-02163E011EF1.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/78335BAA-567D-E711-8038-02163E01A1F2.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/80CF859F-567D-E711-A71D-02163E01A3DE.root",
      "file:/data/veszpv/ZeroBiasRaw2017/300806/8A19A79D-567D-E711-AF63-02163E01A45F.root",
) )

process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('StuckTBMAnalyzer run on reco'),
    name = cms.untracked.string('StuckTBMAnalyzer'),
    version = cms.untracked.string('$alpha$'))

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_dataRun2_Prompt_v1', '')

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:RECOO.root'),
    outputCommands = cms.untracked.vstring('drop *', 'keep *_siPixel*_*_*')
)

process.StuckTBMAnalyzerPlugin = cms.EDAnalyzer('StuckTBMAnalyzer',
    outputFileName = cms.untracked.string("test.root"),
    lumiFileName = cms.untracked.string("lumi300806.cvs"),
    eventSaveDownscaleFactor = cms.untracked.int32(1),
)

process.raw2digi_step = cms.Sequence(process.siPixelDigis)
process.StuckTBMAnalyzer_step = cms.Path(process.raw2digi_step*process.StuckTBMAnalyzerPlugin)

# do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)
