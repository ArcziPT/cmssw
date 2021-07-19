// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingAnalyzerDQM
// Class:      DiamondTimingAnalyzerDQM
//
/**\class DiamondTimingAnalyzerDQM DiamondTimingAnalyzerDQM.cc Analyzer/DiamondTimingAnalyzerDQM/plugins/DiamondTimingAnalyzerDQM.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arkadiusz Wolk
//         Created:  Fri, 16 Jul 2021 14:08:13 GMT
//
//

#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"

#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DiamondDetectorClass.h"

#include <map>
#include <cmath>

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"


//
// class declaration
//

struct Histograms_DiamondTimingAnalyzerDQM {
  TH1F* run_number;

  std::vector<TH1F*>  trk_t_SPC;
	std::vector<TH1F*>  trk_res;
		
	std::map<ChannelKey, double> res_L2_val;
	std::map<ChannelKey, double> res_3p_val;
	std::map<ChannelKey, dqm::reco::MonitorElement*>  res_L2;
	std::map<ChannelKey, dqm::reco::MonitorElement*>  res_L2_3p;
	std::map<ChannelKey, dqm::reco::MonitorElement*>  trk_t_L2; //tracks time
	std::map<ChannelKey, dqm::reco::MonitorElement*>  trk_t_L2_3p;
	std::map<ChannelKey, dqm::reco::MonitorElement*>  trk_exp_res_L2; //tracks expected resolution
	std::map<ChannelKey, dqm::reco::MonitorElement*>  trk_exp_res_L2_3p; //tracks expected resolution (3 planes)
	std::map<std::pair<int,int>, TGraph*>  res_L2_g; //graph
	std::map<std::pair<int,int>, TGraph*>  res_L2_3p_g; //graph
	std::map<std::pair<int,int>, TGraph*>  imp_res_L2_g; //improved resolution L2 graph
	std::map<std::pair<int,int>, TGraph*>  imp_res_L2_3p_g; //improved resolution L2 (3 planes) graph

  std::vector<dqm::reco::MonitorElement*>  tracks_t_SPC_vs_BX;  //<sector>
	std::vector<dqm::reco::MonitorElement*>  tracks_t_SPC_vs_LS;  //<sector>	  
	std::map<int,TProfile*> tracks_t_SPC_vs_BX_prof; //profile
	std::map<int,TProfile*> tracks_t_SPC_vs_LS_prof; //profile

  std::map<ChannelKey, dqm::reco::MonitorElement*> t;
	std::map<ChannelKey, dqm::reco::MonitorElement*> tot;
	std::map<ChannelKey, dqm::reco::MonitorElement*> valid_t;
	std::map<ChannelKey, dqm::reco::MonitorElement*> valid_tot;
	std::map<ChannelKey, dqm::reco::MonitorElement*> tot_vs_t;
	std::map<ChannelKey, dqm::reco::MonitorElement*> pad_tom_220; //pad tomography
	std::map<ChannelKey, dqm::reco::MonitorElement*> pad_tom_210;
};

class DiamondTimingAnalyzerDQM : public DQMGlobalEDAnalyzer<Histograms_DiamondTimingAnalyzerDQM> {
public:
  explicit DiamondTimingAnalyzerDQM(const edm::ParameterSet&);
  ~DiamondTimingAnalyzerDQM() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&,
                      edm::Run const&,
                      edm::EventSetup const&,
                      Histograms_DiamondTimingAnalyzerDQM&) const override;

  void dqmAnalyze(edm::Event const&, edm::EventSetup const&, Histograms_DiamondTimingAnalyzerDQM const&) const override;

  // ------------ member data ------------
  edm::Service<TFileService> fs_;

  // ---------- objects to retrieve ---------------------------
	edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondDigi> > tokenDigi_;
	edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_;
	edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_;
	edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack> > tokenPixelLocalTrack_;
  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;

  int cal_channel_;
	int cal_plane_;
	int cal_sector_;
	double cal_par_0_;
	double cal_par_1_;
	double cal_par_2_;
	double cal_par_3_;
	double resolution_L2_;
	double resolution_L2_3p_;

  //external 
	DiamondDetectorClass DiamondDet;
	int valid_OOT_;

  std::map< std::pair< int , int >, std::pair< int , int > > Ntracks_cuts_map_; //arm, station ,, Lcut,Ucut
};

//
// constants, enums and typedefs
//
static const int CHANNELS_X_PLANE  = 12;
static const int PLANES_X_DETECTOR = 4;
static const int MAX_SECTOR_NUMBER = 2;

enum Sector_id_{
	SECTOR_45_ID,
	SECTOR_56_ID
};
	  
enum Plane_id_{
	PLANE_0_ID,
	PLANE_1_ID,
	PLANE_2_ID,
	PLANE_3_ID
};
	  
enum Station_id_{
	STATION_210_M_ID,
	STATION_TIMING_ID,
	STATION_220_M_ID
};


//
// constructors and destructor
//
DiamondTimingAnalyzerDQM::DiamondTimingAnalyzerDQM(const edm::ParameterSet& iConfig)
  : 
  tokenDigi_(consumes<edm::DetSetVector<CTPPSDiamondDigi>>(iConfig.getParameter<edm::InputTag>("tagDigi"))),
  tokenRecHit_(consumes<edm::DetSetVector<CTPPSDiamondRecHit>>(iConfig.getParameter<edm::InputTag>("tagRecHit"))),
  tokenLocalTrack_(consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagLocalTrack"))),
  tokenPixelLocalTrack_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagPixelLocalTrack"))),
  DiamondDet(iConfig,tokenRecHit_,tokenLocalTrack_),
  valid_OOT_(iConfig.getParameter<int>("tagValidOOT")) {
  
  //usesResource("TFileService"); 
  
	Ntracks_cuts_map_[std::make_pair(SECTOR_45_ID,STATION_210_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[0],
																					 iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[0]);
	Ntracks_cuts_map_[std::make_pair(SECTOR_45_ID,STATION_220_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[1],
																					 iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[1]);
	Ntracks_cuts_map_[std::make_pair(SECTOR_56_ID,STATION_210_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[2],
																					 iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[2]);
	Ntracks_cuts_map_[std::make_pair(SECTOR_56_ID,STATION_220_M_ID)] = std::make_pair(iConfig.getParameter< std::vector <int> >( "Ntracks_Lcuts" )[3],
																					 iConfig.getParameter< std::vector <int> >( "Ntracks_Ucuts" )[3]);

  }

//
// member functions
//

// ------------ method called for each event  ------------

void DiamondTimingAnalyzerDQM::dqmAnalyze(edm::Event const& iEvent,
                           edm::EventSetup const& iSetup,
                           Histograms_DiamondTimingAnalyzerDQM const& histos) const {
}


void DiamondTimingAnalyzerDQM::bookHistograms(DQMStore::IBooker& iBooker,
                               edm::Run const& run,
                               edm::EventSetup const& iSetup,
                               Histograms_DiamondTimingAnalyzerDQM& histos) const {
  iBooker.cd();
  iBooker.setCurrentFolder(".");
  std::string ch_name;

  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (!CTPPSDiamondDetId::check(it->first))
      continue;

    const CTPPSDiamondDetId detid(it->first);

    if (detid.station() != 1)
      continue;

    detid.channelName(ch_name);
    histos.t[detid.rawId()] = iBooker.book1D("t_" + ch_name, ch_name + ";t (ns);Entries", 1200, -60., 60.);
    histos.valid_t[detid.rawId()] = iBooker.book1D("valid_t_" + ch_name, ch_name + ";t (ns);Entries", 1200, -60., 60.);
    histos.tot[detid.rawId()] = iBooker.book1D("tot_" + ch_name, ch_name + ";ToT (ns);Entries", 100, -20., 20.);
    histos.valid_tot[detid.rawId()] = iBooker.book1D("valid_tot_" + ch_name, ch_name + ";ToT (ns);Entries", 100, -20., 20.);
    histos.tot_vs_t[detid.rawId()] =
        iBooker.book2D("tot_vs_t_" + ch_name, ch_name + ";ToT (ns);t (ns)", 240, 0., 60., 450, -20., 25.);
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingAnalyzerDQM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folder", "MY_FOLDER");
  descriptions.add("diamondtiminganalyzerdqm", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingAnalyzerDQM);
