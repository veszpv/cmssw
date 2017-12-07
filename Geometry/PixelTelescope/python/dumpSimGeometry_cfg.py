import FWCore.ParameterSet.Config as cms
import sys, os, operator
import FWCore.ParameterSet.VarParsing as VarParsing
from FWCore.Utilities.Enumerate import Enumerate
from pprint import pprint
from Geometry.PixelTelescope.dictPixelTelescopeGeometry import detectorVersionDict

varType = Enumerate ("Telescope")
defaultVersion=str();

def help():
   print "Usage: cmsRun dumpSimGeometry_cfg.py  tag=TAG version=VERSION "
   print "   tag=tagname"
   print "       indentify geometry scenario "
   print "      ", varType.keys()
   print ""
   print "   version=versionNumber"
   print "       scenario version from PixelTelescope dictionary:"
   print ""
   print "   out=outputFileName"
   print "       default is cmsSimGeom<tag><version>.root"
   print 
   os._exit(1);

def versionCheck(ver):
   if ver == "":
      print "Please, specify PixelTelescope scenario version\n"
      pprint(sorted(detectorVersionDict.items(),key=operator.itemgetter(1)))
      help()

def simGeoLoad(score):
    print "Loading configuration for scenario", options.tag , options.version ,"...\n"
    if score == "Telescope":
       process.load("Geometry.PixelTelescope.Phase2BeamTestXML_cfi")
       
    elif score == "2":
       print "Versioning and dictPixelTelescopeGeometry content is not yet implemented.\n"
       process.load("Geometry.PixelTelescope.Phase2BeamTestXML_cfi")
#       versionCheck(options.version)
#       process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023" + options.version + "XML_cfi")

    else:
      help()

options = VarParsing.VarParsing ()

defaultTag=str();
defaultLevel=14;
defaultOutputFileName="cmsSimGeom-"+ defaultTag +".root"

options.register ('tag',
                  defaultTag, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "info about geometry scenario")
options.register ('version',
                  defaultVersion, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "info about PixelTelescope geometry scenario version")
options.register ('out',
                  defaultOutputFileName, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Output file name")


options.parseArguments()


if (options.out == defaultOutputFileName ):
   options.out = "cmsSimGeom-" + str(options.tag) + str(options.version) + ".root"

process = cms.Process("SIMDUMP")
simGeoLoad(options.tag)

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.add_(cms.ESProducer("TGeoMgrFromDdd",
                            verbose = cms.untracked.bool(False),
                            level = cms.untracked.int32(defaultLevel)
                            ))

process.dump = cms.EDAnalyzer("DumpSimGeometry", 
                              tag = cms.untracked.string(options.tag),
                              outputFileName = cms.untracked.string(options.out))

process.p = cms.Path(process.dump)
