import FWCore.ParameterSet.Config as cms

diamondTimingAnalyzer = cms.EDProducer("DiamondTimingAnalyzer",
    tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
    tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
    tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
    tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"), 
    tagCalibrationFile = cms.string("src/Analyzer/DiamondTimingAnalyzer/python/DiamondCalibration.json"),
    tagValidOOT = cms.int32(-1),
    folder = cms.string("myfolder"),
    Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
    Ntracks_Ucuts = cms.vint32([-1,6,-1,6]),
)