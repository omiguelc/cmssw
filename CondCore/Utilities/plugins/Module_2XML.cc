#include "CondCore/Utilities/interface/Payload2XMLModule.h"
#include "CondCore/Utilities/src/CondFormats.h"

PAYLOAD_2XML_MODULE(pluginUtilities_payload2xml) {
  m.def("boost_version_label", &cond::boost_version_label, "Get boost version for this release"); \
  PAYLOAD_2XML_CLASS(AlCaRecoTriggerBits);
  PAYLOAD_2XML_CLASS(AlignPCLThresholds);
  PAYLOAD_2XML_CLASS(AlignmentErrors);
  PAYLOAD_2XML_CLASS(AlignmentErrorsExtended);
  PAYLOAD_2XML_CLASS(AlignmentSurfaceDeformations);
  PAYLOAD_2XML_CLASS(Alignments);
  PAYLOAD_2XML_CLASS(BeamSpotObjects);
  PAYLOAD_2XML_CLASS(BeamSpotOnlineObjects);
  PAYLOAD_2XML_CLASS(CSCBadChambers);
  PAYLOAD_2XML_CLASS(CSCBadStrips);
  PAYLOAD_2XML_CLASS(CSCBadWires);
  PAYLOAD_2XML_CLASS(CSCChamberIndex);
  PAYLOAD_2XML_CLASS(CSCChamberMap);
  PAYLOAD_2XML_CLASS(CSCChamberTimeCorrections);
  PAYLOAD_2XML_CLASS(CSCCrateMap);
  PAYLOAD_2XML_CLASS(CSCDBChipSpeedCorrection);
  PAYLOAD_2XML_CLASS(CSCDBCrosstalk);
  PAYLOAD_2XML_CLASS(CSCDBGains);
  PAYLOAD_2XML_CLASS(CSCDBGasGainCorrection);
  PAYLOAD_2XML_CLASS(CSCDBL1TPParameters);
  PAYLOAD_2XML_CLASS(CSCDBNoiseMatrix);
  PAYLOAD_2XML_CLASS(CSCDBPedestals);
  PAYLOAD_2XML_CLASS(CSCDDUMap);
  PAYLOAD_2XML_CLASS(CSCL1TPParameters);
  PAYLOAD_2XML_CLASS(CSCRecoDigiParameters);
  PAYLOAD_2XML_CLASS(CTPPSPixelAnalysisMask);
  PAYLOAD_2XML_CLASS(CTPPSPixelDAQMapping);
  PAYLOAD_2XML_CLASS(CTPPSPixelGainCalibrations);
  PAYLOAD_2XML_CLASS(PPSAlignmentConfig)
  PAYLOAD_2XML_CLASS(PPSAlignmentConfiguration)
  PAYLOAD_2XML_CLASS(PPSAssociationCuts)
  PAYLOAD_2XML_CLASS(CastorChannelQuality);
  PAYLOAD_2XML_CLASS(CastorElectronicsMap);
  PAYLOAD_2XML_CLASS(CastorGainWidths);
  PAYLOAD_2XML_CLASS(CastorGains);
  PAYLOAD_2XML_CLASS(CastorPedestalWidths);
  PAYLOAD_2XML_CLASS(CastorPedestals);
  PAYLOAD_2XML_CLASS(CastorQIEData);
  PAYLOAD_2XML_CLASS(CastorRecoParams);
  PAYLOAD_2XML_CLASS(CastorSaturationCorrs);
  PAYLOAD_2XML_CLASS(CentralityTable);
  PAYLOAD_2XML_CLASS(DTCCBConfig);
  PAYLOAD_2XML_CLASS(DTDeadFlag);
  PAYLOAD_2XML_CLASS(DTHVStatus);
  PAYLOAD_2XML_CLASS(DTKeyedConfig);
  PAYLOAD_2XML_CLASS(DTLVStatus);
  PAYLOAD_2XML_CLASS(DTMtime);
  PAYLOAD_2XML_CLASS(DTReadOutMapping);
  PAYLOAD_2XML_CLASS(DTRecoConditions);
  PAYLOAD_2XML_CLASS(DTRecoUncertainties);
  PAYLOAD_2XML_CLASS(DTStatusFlag);
  PAYLOAD_2XML_CLASS(DTT0);
  PAYLOAD_2XML_CLASS(DTTPGParameters);
  PAYLOAD_2XML_CLASS(DTTtrig);
  PAYLOAD_2XML_CLASS(DYTParamObject);
  PAYLOAD_2XML_CLASS(DYTThrObject);
  PAYLOAD_2XML_CLASS(DropBoxMetadata);
  PAYLOAD_2XML_CLASS(ESCondObjectContainer<ESChannelStatusCode>);
  PAYLOAD_2XML_CLASS(ESCondObjectContainer<ESPedestal>);
  PAYLOAD_2XML_CLASS(ESCondObjectContainer<float>);
  PAYLOAD_2XML_CLASS(ESEEIntercalibConstants);
  PAYLOAD_2XML_CLASS(ESGain);
  PAYLOAD_2XML_CLASS(ESMIPToGeVConstant);
  PAYLOAD_2XML_CLASS(ESMissingEnergyCalibration);
  PAYLOAD_2XML_CLASS(ESRecHitRatioCuts);
  PAYLOAD_2XML_CLASS(ESThresholds);
  PAYLOAD_2XML_CLASS(ESTimeSampleWeights);
  PAYLOAD_2XML_CLASS(EcalADCToGeVConstant);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalChannelStatusCode>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalDQMStatusCode>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalMGPAGainRatio>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalMappingElement>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalPedestal>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalPulseCovariance>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalPulseShape>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalPulseSymmCovariance>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalTPGCrystalStatusCode>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalTPGLinearizationConstant>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalTPGPedestal>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<EcalXtalGroupId>);
  PAYLOAD_2XML_CLASS(EcalCondObjectContainer<float>);
  PAYLOAD_2XML_CLASS(EcalCondTowerObjectContainer<EcalChannelStatusCode>);
  PAYLOAD_2XML_CLASS(EcalCondTowerObjectContainer<EcalDAQStatusCode>);
  PAYLOAD_2XML_CLASS(EcalCondTowerObjectContainer<EcalDQMStatusCode>);
  PAYLOAD_2XML_CLASS(EcalFunParams);
  PAYLOAD_2XML_CLASS(EcalLaserAPDPNRatios);
  PAYLOAD_2XML_CLASS(EcalMustacheSCParameters);
  PAYLOAD_2XML_CLASS(EcalSCDynamicDPhiParameters);
  PAYLOAD_2XML_CLASS(EcalSRSettings);
  PAYLOAD_2XML_CLASS(EcalSampleMask);
  PAYLOAD_2XML_CLASS(EcalSamplesCorrelation);
  PAYLOAD_2XML_CLASS(EcalSimPulseShape);
  PAYLOAD_2XML_CLASS(EcalTBWeights);
  PAYLOAD_2XML_CLASS(EcalTPGFineGrainEBGroup);
  PAYLOAD_2XML_CLASS(EcalTPGFineGrainEBIdMap);
  PAYLOAD_2XML_CLASS(EcalTPGFineGrainStripEE);
  PAYLOAD_2XML_CLASS(EcalTPGFineGrainTowerEE);
  PAYLOAD_2XML_CLASS(EcalTPGLutGroup);
  PAYLOAD_2XML_CLASS(EcalTPGLutIdMap);
  PAYLOAD_2XML_CLASS(EcalTPGPhysicsConst);
  PAYLOAD_2XML_CLASS(EcalTPGSlidingWindow);
  PAYLOAD_2XML_CLASS(EcalTPGSpike);
  PAYLOAD_2XML_CLASS(EcalTPGStripStatus);
  PAYLOAD_2XML_CLASS(EcalTPGTowerStatus);
  PAYLOAD_2XML_CLASS(EcalTPGWeightGroup);
  PAYLOAD_2XML_CLASS(EcalTPGWeightIdMap);
  PAYLOAD_2XML_CLASS(EcalTPGOddWeightGroup);
  PAYLOAD_2XML_CLASS(EcalTPGOddWeightIdMap);
  PAYLOAD_2XML_CLASS(EcalTPGTPMode);
  PAYLOAD_2XML_CLASS(EcalTimeBiasCorrections);
  PAYLOAD_2XML_CLASS(EcalTimeDependentCorrections);
  PAYLOAD_2XML_CLASS(EcalTimeOffsetConstant);
  PAYLOAD_2XML_CLASS(FileBlob);
  PAYLOAD_2XML_CLASS(GBRForest);
  PAYLOAD_2XML_CLASS(GBRForestD);
  //PAYLOAD_2XML_CLASS( HBHENegativeEFilter );
  PAYLOAD_2XML_CLASS(HcalChannelQuality);
  PAYLOAD_2XML_CLASS(HcalDcsValues);
  PAYLOAD_2XML_CLASS(HcalElectronicsMap);
  PAYLOAD_2XML_CLASS(HcalFlagHFDigiTimeParams);
  PAYLOAD_2XML_CLASS(HcalFrontEndMap);
  PAYLOAD_2XML_CLASS(HcalGainWidths);
  PAYLOAD_2XML_CLASS(HcalGains);
  //PAYLOAD_2XML_CLASS( HcalInterpolatedPulseColl );
  //PAYLOAD_2XML_CLASS( HcalItemCollById<HFPhase1PMTData> );
  PAYLOAD_2XML_CLASS(HcalL1TriggerObjects);
  PAYLOAD_2XML_CLASS(HcalLUTCorrs);
  PAYLOAD_2XML_CLASS(HcalLongRecoParams);
  PAYLOAD_2XML_CLASS(HcalLutMetadata);
  PAYLOAD_2XML_CLASS(HcalMCParams);
  //PAYLOAD_2XML_CLASS( HFPhase1PMTParams );
  PAYLOAD_2XML_CLASS(HcalPFCorrs);
  PAYLOAD_2XML_CLASS(HcalParameters);
  PAYLOAD_2XML_CLASS(HcalPedestalWidths);
  PAYLOAD_2XML_CLASS(HcalPedestals);
  PAYLOAD_2XML_CLASS(HcalQIEData);
  PAYLOAD_2XML_CLASS(HcalQIETypes);
  PAYLOAD_2XML_CLASS(HcalRecoParams);
  PAYLOAD_2XML_CLASS(HcalRespCorrs);
  PAYLOAD_2XML_CLASS(HcalSiPMCharacteristics);
  PAYLOAD_2XML_CLASS(HcalSiPMParameters);
  PAYLOAD_2XML_CLASS(HcalTPChannelParameters);
  PAYLOAD_2XML_CLASS(HcalTPParameters);
  PAYLOAD_2XML_CLASS(HcalTimeCorrs);
  PAYLOAD_2XML_CLASS(HcalZDCLowGainFractions);
  PAYLOAD_2XML_CLASS(HcalZSThresholds);
  PAYLOAD_2XML_CLASS(JME::JetResolutionObject);
  PAYLOAD_2XML_CLASS(JetCorrectorParametersCollection);
  PAYLOAD_2XML_CLASS(L1CaloEcalScale);
  PAYLOAD_2XML_CLASS(L1CaloEtScale);
  PAYLOAD_2XML_CLASS(L1CaloGeometry);
  PAYLOAD_2XML_CLASS(L1CaloHcalScale);
  PAYLOAD_2XML_CLASS(L1GctChannelMask);
  PAYLOAD_2XML_CLASS(L1GctJetFinderParams);
  PAYLOAD_2XML_CLASS(L1GtBoardMaps);
  PAYLOAD_2XML_CLASS(L1GtParameters);
  PAYLOAD_2XML_CLASS(L1GtPrescaleFactors);
  PAYLOAD_2XML_CLASS(L1GtPsbSetup);
  PAYLOAD_2XML_CLASS(L1GtStableParameters);
  PAYLOAD_2XML_CLASS(L1GtTriggerMask);
  PAYLOAD_2XML_CLASS(L1GtTriggerMenu);
  PAYLOAD_2XML_CLASS(L1MuCSCPtLut);
  PAYLOAD_2XML_CLASS(L1MuCSCTFAlignment);
  PAYLOAD_2XML_CLASS(L1MuCSCTFConfiguration);
  PAYLOAD_2XML_CLASS(L1MuDTEtaPatternLut);
  PAYLOAD_2XML_CLASS(L1MuDTExtLut);
  PAYLOAD_2XML_CLASS(L1MuDTPhiLut);
  PAYLOAD_2XML_CLASS(L1MuDTPtaLut);
  PAYLOAD_2XML_CLASS(L1MuDTQualPatternLut);
  PAYLOAD_2XML_CLASS(L1MuDTTFMasks);
  PAYLOAD_2XML_CLASS(L1MuDTTFParameters);
  PAYLOAD_2XML_CLASS(L1MuGMTChannelMask);
  PAYLOAD_2XML_CLASS(L1MuGMTParameters);
  PAYLOAD_2XML_CLASS(L1MuGMTScales);
  PAYLOAD_2XML_CLASS(L1MuTriggerPtScale);
  PAYLOAD_2XML_CLASS(L1MuTriggerScales);
  PAYLOAD_2XML_CLASS(L1RCTChannelMask);
  PAYLOAD_2XML_CLASS(L1RCTNoisyChannelMask);
  PAYLOAD_2XML_CLASS(L1RCTParameters);
  PAYLOAD_2XML_CLASS(L1RPCBxOrConfig);
  PAYLOAD_2XML_CLASS(L1RPCConeDefinition);
  PAYLOAD_2XML_CLASS(L1RPCConfig);
  PAYLOAD_2XML_CLASS(L1RPCHsbConfig);
  PAYLOAD_2XML_CLASS(L1RPCHwConfig);
  PAYLOAD_2XML_CLASS(L1TGlobalParameters);
  PAYLOAD_2XML_CLASS(L1TGlobalPrescalesVetos);
  PAYLOAD_2XML_CLASS(L1TGlobalPrescalesVetosFract);
  PAYLOAD_2XML_CLASS(L1TMuonBarrelParams);
  PAYLOAD_2XML_CLASS(L1TMuonEndCapForest);
  PAYLOAD_2XML_CLASS(L1TMuonEndCapParams);
  PAYLOAD_2XML_CLASS(L1TMuonGlobalParams);
  PAYLOAD_2XML_CLASS(L1TMuonOverlapParams);
  PAYLOAD_2XML_CLASS(L1TUtmAlgorithm);
  PAYLOAD_2XML_CLASS(L1TUtmBin);
  PAYLOAD_2XML_CLASS(L1TUtmCondition);
  PAYLOAD_2XML_CLASS(L1TUtmCut);
  PAYLOAD_2XML_CLASS(L1TUtmCutValue);
  PAYLOAD_2XML_CLASS(L1TUtmObject);
  PAYLOAD_2XML_CLASS(L1TUtmScale);
  PAYLOAD_2XML_CLASS(L1TUtmTriggerMenu);
  PAYLOAD_2XML_CLASS(L1TriggerKey);
  PAYLOAD_2XML_CLASS(L1TriggerKeyList);
  PAYLOAD_2XML_CLASS(LHCInfo);
  PAYLOAD_2XML_CLASS(METCorrectorParametersCollection);
  PAYLOAD_2XML_CLASS(MEtXYcorrectParametersCollection);
  PAYLOAD_2XML_CLASS(MagFieldConfig);
  PAYLOAD_2XML_CLASS(MixingModuleConfig);
  PAYLOAD_2XML_CLASS(MuScleFitDBobject);
  PAYLOAD_2XML_CLASS(OOTPileupCorrectionBuffer);
  PAYLOAD_2XML_CLASS(PCaloGeometry);
  PAYLOAD_2XML_CLASS(PGeometricDet);
  PAYLOAD_2XML_CLASS(PHGCalParameters);
  PAYLOAD_2XML_CLASS(PTrackerParameters);
  PAYLOAD_2XML_CLASS(PTrackerAdditionalParametersPerDet);
  PAYLOAD_2XML_CLASS(PerformancePayloadFromBinnedTFormula);
  PAYLOAD_2XML_CLASS(PerformancePayloadFromTFormula);
  PAYLOAD_2XML_CLASS(PerformancePayloadFromTable);
  PAYLOAD_2XML_CLASS(PerformanceWorkingPoint);
  PAYLOAD_2XML_CLASS(PhysicsTFormulaPayload);
  PAYLOAD_2XML_CLASS(PhysicsTGraphPayload);
  //PAYLOAD_2XML_CLASS( PhysicsTools::Calibration::Histogram3D<double, double, double, double> );
  PAYLOAD_2XML_CLASS(PhysicsTools::Calibration::HistogramD3D);
  PAYLOAD_2XML_CLASS(PhysicsTools::Calibration::MVAComputerContainer);
  PAYLOAD_2XML_CLASS(QGLikelihoodCategory);
  PAYLOAD_2XML_CLASS(QGLikelihoodObject);
  PAYLOAD_2XML_CLASS(QGLikelihoodSystematicsObject);
  PAYLOAD_2XML_CLASS(RPCAMCLinkMap);
  PAYLOAD_2XML_CLASS(RPCClusterSize);
  PAYLOAD_2XML_CLASS(RPCDCCLinkMap);
  PAYLOAD_2XML_CLASS(RPCEMap);
  PAYLOAD_2XML_CLASS(RPCLBLinkMap);
  PAYLOAD_2XML_CLASS(RPCObFebmap);
  PAYLOAD_2XML_CLASS(RPCObGas);
  PAYLOAD_2XML_CLASS(RPCObGasMix);
  PAYLOAD_2XML_CLASS(RPCObImon);
  PAYLOAD_2XML_CLASS(RPCObPVSSmap);
  PAYLOAD_2XML_CLASS(RPCObStatus);
  PAYLOAD_2XML_CLASS(RPCObTemp);
  PAYLOAD_2XML_CLASS(RPCObUXC);
  PAYLOAD_2XML_CLASS(RPCObVmon);
  PAYLOAD_2XML_CLASS(RPCStripNoises);
  PAYLOAD_2XML_CLASS(RPFlatParams);
  PAYLOAD_2XML_CLASS(RecoIdealGeometry);
  PAYLOAD_2XML_CLASS(RunInfo);
  PAYLOAD_2XML_CLASS(SiPhase2OuterTrackerLorentzAngle);
  PAYLOAD_2XML_CLASS(SiPixel2DTemplateDBObject);
  PAYLOAD_2XML_CLASS(SiPixelCPEGenericErrorParm);
  PAYLOAD_2XML_CLASS(SiPixelCalibConfiguration);
  PAYLOAD_2XML_CLASS(SiPixelDynamicInefficiency);
  PAYLOAD_2XML_CLASS(SiPixelFedCablingMap);
  PAYLOAD_2XML_CLASS(SiPixelGainCalibrationForHLT);
  PAYLOAD_2XML_CLASS(SiPixelGainCalibrationOffline);
  PAYLOAD_2XML_CLASS(SiPixelGenErrorDBObject);
  PAYLOAD_2XML_CLASS(SiPixelLorentzAngle);
  PAYLOAD_2XML_CLASS(SiPixelQuality);
  PAYLOAD_2XML_CLASS(SiPixelFEDChannelContainer);
  PAYLOAD_2XML_CLASS(SiPixelQualityProbabilities);
  PAYLOAD_2XML_CLASS(SiPixelVCal);
  PAYLOAD_2XML_CLASS(SiPixelTemplateDBObject);
  PAYLOAD_2XML_CLASS(SiStripApvGain);
  PAYLOAD_2XML_CLASS(SiStripApvSimulationParameters);
  PAYLOAD_2XML_CLASS(SiStripBackPlaneCorrection);
  PAYLOAD_2XML_CLASS(SiStripBadStrip);
  PAYLOAD_2XML_CLASS(SiStripConfObject);
  PAYLOAD_2XML_CLASS(SiStripDetVOff);
  PAYLOAD_2XML_CLASS(SiStripFedCabling);
  PAYLOAD_2XML_CLASS(SiStripLatency);
  PAYLOAD_2XML_CLASS(SiStripLorentzAngle);
  PAYLOAD_2XML_CLASS(SiStripNoises);
  PAYLOAD_2XML_CLASS(SiStripPedestals);
  PAYLOAD_2XML_CLASS(SiStripThreshold);
  PAYLOAD_2XML_CLASS(DTCELinkId);
  PAYLOAD_2XML_CLASS(TrackerDetToDTCELinkCablingMap);
  //PAYLOAD_2XML_CLASS( StorableDoubleMap<AbsOOTPileupCorrection> );
  PAYLOAD_2XML_CLASS(TrackProbabilityCalibration);
  PAYLOAD_2XML_CLASS(cond::BaseKeyed);
  PAYLOAD_2XML_CLASS(l1t::CaloConfig);
  PAYLOAD_2XML_CLASS(l1t::CaloParams);
  PAYLOAD_2XML_CLASS(lumi::LumiSectionData);
  PAYLOAD_2XML_CLASS(std::vector<unsigned long long>);
}
