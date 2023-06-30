import FWCore.ParameterSet.Config as cms

XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring(
        'Geometry/CMSCommonData/data/materials/2021/v1/materials.xml',
        'Geometry/CMSCommonData/data/rotations.xml',
        'Geometry/CMSCommonData/data/extend/v2/cmsextent.xml',
        'Geometry/CMSCommonData/data/cavernData/2021/v1/cavernData.xml',
        'Geometry/CMSCommonData/data/cms/2026/v5/cms.xml',
        'Geometry/CMSCommonData/data/cmsMother.xml',
        'Geometry/CMSCommonData/data/eta3/etaMax.xml',
        'Geometry/CMSCommonData/data/cmsTracker.xml',
        'Geometry/CMSCommonData/data/caloBase/2026/v7/caloBase.xml',
        'Geometry/CMSCommonData/data/cmsCalo.xml',
        'Geometry/CMSCommonData/data/muonBase/2026/v5/muonBase.xml',
        'Geometry/CMSCommonData/data/cmsMuon.xml',
        'Geometry/CMSCommonData/data/mgnt.xml',
        'Geometry/CMSCommonData/data/beampipe/2026/v3/beampipe.xml',
        'Geometry/CMSCommonData/data/cmsBeam/2026/v1/cmsBeam.xml',
        'Geometry/CMSCommonData/data/muonMB.xml',
        'Geometry/CMSCommonData/data/muonMagnet.xml',
        'Geometry/CMSCommonData/data/cavern/2021/v1/cavern.xml',
        'Geometry/CMSCommonData/data/cavernFloor/2017/v1/cavernFloor.xml',
        'Geometry/TrackerCommonData/data/PhaseII/trackerParameters.xml',
        'Geometry/TrackerCommonData/data/pixfwdCommon.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_2021_02/pixfwd.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_OT800_IT615_2022_10/pixbar.xml',
        'Geometry/TrackerCommonData/data/trackermaterial.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_2021_02/tracker.xml',
        'Geometry/TrackerCommonData/data/PhaseII/OuterTracker616_2020_04/otst.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_IT702_2021_03/pixel.xml',
        'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker404/trackerbar.xml',
        'Geometry/TrackerCommonData/data/PhaseII/TiltedTracker404/trackerfwd.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_2021_02/trackerStructureTopology.xml',
        'Geometry/TrackerCommonData/data/PhaseII/Tracker_DD4hep_compatible_IT702_2021_03/pixelStructureTopology.xml',
        'Geometry/TrackerSimData/data/PhaseII/Tracker_DD4hep_compatible_2021_02/trackersens.xml',
        'Geometry/TrackerSimData/data/PhaseII/Tracker_DD4hep_compatible_IT702_2021_03/pixelsens.xml',
        'Geometry/TrackerRecoData/data/PhaseII/Tracker_DD4hep_compatible_IT702_2021_03/trackerRecoMaterial.xml',
        'SimTracker/TrackerMaterialAnalysis/data/trackingMaterialGroups_ForPhaseII/v2/trackingMaterialGroups_ForPhaseII.xml',
        'Geometry/TrackerSimData/data/PhaseII/Tracker_DD4hep_compatible_2021_02/trackerProdCuts.xml',
        'Geometry/TrackerSimData/data/PhaseII/Tracker_DD4hep_compatible_IT702_2021_03/pixelProdCuts.xml',
        'Geometry/TrackerSimData/data/trackerProdCutsBEAM.xml',
        'Geometry/CMSCommonData/data/FieldParameters.xml',
    ),
    rootNodeName = cms.string('cms:OCMS')
)
