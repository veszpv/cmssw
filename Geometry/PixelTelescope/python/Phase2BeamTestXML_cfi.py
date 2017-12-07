import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring(
         # World volume creation and default CMS materials
        'Geometry/PixelTelescope/data/materials.xml',
        'Geometry/PixelTelescope/data/rotations.xml',
        'Geometry/PixelTelescope/data/extend/cmsextent.xml',
        'Geometry/PixelTelescope/data/cms.xml',
        'Geometry/PixelTelescope/data/cmsMother.xml',
        
         # Define standalone Phase 1 BPIX module and associated materials
        'Geometry/PixelTelescope/data/pixbarmaterial.xml',
        'Geometry/PixelTelescope/data/trackermaterial.xml',   #sometimes included via Run2, sometimes via PhaseI, which one to choose?
        'Geometry/PixelTelescope/data/pixfwdMaterials.xml',
        'Geometry/PixelTelescope/data/Phase1BPIXLayer4Module.xml',
        
         # Define DUT and telescope
        'Geometry/PixelTelescope/data/DUT.xml',
        'Geometry/PixelTelescope/data/telescope.xml',
        
         # Configurable parameters
        'Geometry/PixelTelescope/data/Phase2BeamTestConstants.xml'
      
    ),
    rootNodeName = cms.string('cms:OCMS')
)
