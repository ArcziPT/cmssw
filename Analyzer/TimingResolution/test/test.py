import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *
import string

process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)


#options.register ('slopeAddA',
#				  0,
#				  VarParsing.multiplicity.singleton,
#				  VarParsing.varType.float,
#				  "Slope to be added on first channel")
				  
				  
options.parseArguments()
				  

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
#process.load("DQM.Integration.config.environment_cfi")
#process.dqmEnv.subSystemFolder = "CTPPS"
#process.dqmEnv.eventInfoFolder = "EventInfo"
#process.dqmSaver.path = ""
#process.dqmSaver.tag = "CTPPS"

# process.source = cms.Source("NewEventStreamFileReader",
#     fileNames = cms.untracked.vstring(
#     *(
#     'file:///afs/cern.ch/user/n/nminafra/Work/public/SampicCMS/run319190_ls0130_streamPhysicsTOTEM3_StorageManager.dat',
#     )
#     )
# )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    *(
    #'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/D4CDC51F-9486-E811-8D08-FA163E0E63D7.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/B4F35FF7-9986-E811-AD21-FA163E619822.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/9260F9EF-9986-E811-B6FE-FA163EF4135C.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/76EF31D1-9986-E811-A28F-FA163EEB45BD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/6EB624BF-9086-E811-8DB6-FA163E306BAD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/12EE7D0D-9A86-E811-92E2-FA163E04B59.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/00302FB3-1CB4-E811-96B2-FA163EC9F152.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/0246C0FF-21B4-E811-B82D-FA163E4A6585.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/024F9C3C-1EB4-E811-84C9-FA163EFF8B83.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/0645A771-2AB4-E811-95A8-FA163E4C2520.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/08EB25DD-2EB4-E811-9C78-02163E00B262.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/0CFD4E19-22B4-E811-B2C2-02163E016469.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/10335B7A-2AB4-E811-A125-FA163ED383EE.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/10B9DB48-1EB4-E811-956C-FA163ECAADB8.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/129AA2CD-1CB4-E811-9E4F-FA163EDA78EE.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/16ED4249-18B4-E811-8446-02163E01A068.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/1867BEB4-B9B4-E811-BD84-02163E010D99.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/1A21EF36-1EB4-E811-9871-FA163ED8A558.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/1AA81129-1EB4-E811-BCB5-FA163E0ED030.root'
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/20499182-0C78-E811-AC77-A4BF0112BD52.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/6E5F383F-0D78-E811-A8E6-FA163ECDEF72.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/BE61E668-0B78-E811-BD3D-FA163E60E3B4.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/EE9F884F-0A78-E811-8B9E-FA163EA96792.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/AA21D02F-A784-E811-8530-A4BF01277792.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/AA07FE1A-B284-E811-B673-FA163EBC0B5E.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A84D3A74-B084-E811-862D-FA163E51AFD3.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A6FCDC7C-AD84-E811-B784-02163E01A129.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A66307F6-AC84-E811-8F61-FA163E4C942A.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A65749B5-AD84-E811-AB6F-FA163EA6A331.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A6010C62-AD84-E811-B972-FA163E373CC4.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A4202676-AD84-E811-8E86-FA163E49296B.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A259320A-A984-E811-B227-FA163E1051E6.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A22F0E55-AC84-E811-B791-FA163EA8044F.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A04C1629-AB84-E811-A323-FA163ED286F5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A02A9F7A-B184-E811-8DD1-FA163EFD7D19.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9EBD908F-AB84-E811-B035-02163E019ED0.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9EB80A05-AD84-E811-A216-FA163EE0CB15.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9E92D3C0-AD84-E811-93BB-FA163EA1E4FD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9E7EAFDF-B284-E811-985A-02163E010C30.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9C7858F6-AF84-E811-8DE6-02163E01507D.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9C78511E-AD84-E811-BC8D-FA163E126D94.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/989BC9F0-AC84-E811-A4D8-02163E01800D.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/98963900-AB84-E811-A42E-FA163E89DF43.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/968BA940-AD84-E811-B922-FA163ECE132A.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/94CEC45E-B084-E811-96EE-FA163EB83185.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/90C9F87E-B084-E811-8883-02163E019F71.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/90A68047-AB84-E811-A7D7-FA163ECEF311.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/904287C1-A784-E811-A52C-FA163E6369F5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8EC1E27E-AD84-E811-8A42-FA163E862C57.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8CCDB182-AD84-E811-A059-FA163EF4E283.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8CB1EC5E-AD84-E811-B45E-FA163E73F4B9.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C8271AE-AD84-E811-A5F3-02163E00BF95.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C6ABA46-AD84-E811-A048-FA163E7B6908.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C004A69-B084-E811-8A4C-FA163E55E9B1.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8ADD616A-AD84-E811-B94C-FA163E52AA20.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/86736060-AE84-E811-9E72-FA163EFF9AD1.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/862AD75A-AE84-E811-8BC5-FA163E376AE7.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8409D191-AF84-E811-8E16-FA163EFE0B4E.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/82741E0F-A984-E811-913C-FA163ECFE6A5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/820B7559-AD84-E811-BA28-FA163EF96BC3.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/00183235-C47C-E811-ADB7-FA163E06876F.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/02C21249-BF7C-E811-9D6F-FA163E370861.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/0668D8CA-C37C-E811-A1D1-FA163E399762.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/08DBB77F-C97C-E811-92DB-FA163E4A4299.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/0C70B3EF-C47C-E811-991F-FA163E12E196.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/10792C91-BB7C-E811-8FB3-FA163EC10AF6.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/2488BF83-C77C-E811-B6EE-FA163EFDDA6B.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3696483B-C67C-E811-A7B2-FA163E26BFDA.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3C9CB4AE-C17C-E811-9902-02163E019FE8.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3E1F414D-C67C-E811-90D7-FA163E160041.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/425EC7AA-C77C-E811-9458-02163E019F9D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/44111C17-C57C-E811-8E15-FA163EBD447D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/4ABA2487-C07C-E811-84E7-02163E01A179.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/4ED1812B-BD7C-E811-A217-FA163EA5977C.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/5088CA45-B07C-E811-BEE6-FA163EB83F5F.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/5E1F4173-C17C-E811-81F2-FA163E14F11E.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6236C7F1-BF7C-E811-A4A4-02163E01A07C.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6294EA44-C37C-E811-A096-FA163EFE36BC.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6A092290-B97C-E811-B848-FA163E3A303D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6AB60FBD-C47C-E811-AB7A-FA163EA2F3FC.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/70B5D057-C37C-E811-A904-FA163E160041.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/72454809-C27C-E811-8CEE-FA163EDA08D4.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/7681FF63-C57C-E811-B0DC-FA163E850BF4.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/786B16E7-BC7C-E811-99B7-FA163EC1D4C0.root',
    )
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_HLT_v7', '')



# rechits production
process.load('Geometry.VeryForwardGeometry.geometryRPFromDD_2018_cfi')
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff')

process.diamondTimingAnalyzer = cms.EDAnalyzer('DiamondTimingAnalyzer',
 tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
 tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
 tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"),
 Ntracks_Lcuts = cms.vint32([-1,-1,-1,-1]),
 Ntracks_Ucuts = cms.vint32([-1,-1,-1,-1]),
 tagCalibrationFile = cms.string('src/Analyzer/Calibrations/DiamondCalibration.json'),
 tagValidOOT = cms.int32(-1)
) 

 #tagTimingLocalTrack = cms.InputTag("ctppsDiamondLocalTrack"),
 #TTracksMux_Lcuts = cms.vint32([-1,-1]),
 #TTracksMux_Ucuts = cms.vint32([-1,-1]), 

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)


process.path = cms.Path(process.diamondTimingAnalyzer)
process.schedule = cms.Schedule(process.path)