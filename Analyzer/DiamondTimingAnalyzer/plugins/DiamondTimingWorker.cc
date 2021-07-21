// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingWorker
// Class:      DiamondTimingWorker
//
/**\class DiamondTimingWorker DiamondTimingWorker.cc Analyzer/DiamondTimingWorker/plugins/DiamondTimingWorker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arkadiusz Wolk
//         Created:  Mon, 19 Jul 2021 13:22:33 GMT
//
//

#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DiamondDetectorClass.h"

//
// class declaration
//
class DiamondTimingWorker : public DQMEDAnalyzer {
public:
  explicit DiamondTimingWorker(const edm::ParameterSet&);
  ~DiamondTimingWorker() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker &iBooker, edm::Run const &, edm::EventSetup const &iSetup) override;
  void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) override;

  // ---------- objects to retrieve ---------------------------
  edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondDigi>> tokenDigi_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondRecHit>> tokenRecHit_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack>> tokenLocalTrack_;
  edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>> tokenPixelLocalTrack_;
  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;

  // ------------ member data ------------
  struct Histograms_DiamondTiming {
    std::map<uint32_t, MonitorElement*> t;
    std::map<uint32_t, MonitorElement*> valid_t;
    std::map<uint32_t, MonitorElement*> tot;
    std::map<uint32_t, MonitorElement*> valid_tot;
    std::map<uint32_t, MonitorElement*> t_vs_tot;
  };
  Histograms_DiamondTiming histos;

  DiamondDetectorClass DiamondDet;
  int validOOT;
  std::map< std::pair< int , int >, std::pair< int , int > > Ntracks_cuts_map_; //arm, station ,, Lcut,Ucut
};

//
// constants, enums and typedefs
//
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
// static data member definitions
//

//
// constructors and destructor
//
DiamondTimingWorker::DiamondTimingWorker(const edm::ParameterSet& iConfig)
  :
  tokenDigi_(consumes<edm::DetSetVector<CTPPSDiamondDigi>>(iConfig.getParameter<edm::InputTag>("tagDigi"))),
  tokenRecHit_(consumes<edm::DetSetVector<CTPPSDiamondRecHit>>(iConfig.getParameter<edm::InputTag>("tagRecHit"))),
  tokenLocalTrack_(consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagLocalTrack"))),
  tokenPixelLocalTrack_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagPixelLocalTrack"))),
  geomEsToken_(esConsumes<edm::Transition::BeginRun>()),
  DiamondDet(iConfig,tokenRecHit_,tokenLocalTrack_),
  validOOT(iConfig.getParameter<int>("tagValidOOT")){
  
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

void DiamondTimingWorker::analyze(const edm::Event &iEvent,
                                  const edm::EventSetup &iSetup){
  using namespace edm;
  
  //retrieve data
  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > timingRecHit;
  edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelLocalTrack;

  iEvent.getByToken(tokenRecHit_, timingRecHit );
  iEvent.getByToken(tokenPixelLocalTrack_, pixelLocalTrack );
  

  ////////////////////////////////////////////////////////////////
  //
  //		EXTRACT PIXELS TRACK NUMBER
  //      Will be used for sector independent event selection
  //
  ///////////////////////////////////////////////////////////////// 
  std::map< std::pair< int , int >, int> Pixel_Mux_map_; //arm, station
  Pixel_Mux_map_.clear();
  std::vector<bool> Sector_TBA(2,true);
  
  for(const auto& RP_trks : *pixelLocalTrack){ //array of tracks
    const CTPPSDetId detId( RP_trks.detId() );
	  //std::cout << "Tracks in arm " << detId.arm() << ", station " << detId.station() << ", rp " << detId.rp() << std::endl;
      
	  for(const auto& trk : RP_trks) {
		  if(!trk.isValid()) continue;
		    Pixel_Mux_map_[ std::make_pair(detId.arm(), detId.station()) ]++;
    }	 
	} 
	
  for (const auto& Ntracks_cuts_iter_ :  Ntracks_cuts_map_){
    if((Ntracks_cuts_iter_.second.first < 0) || (Ntracks_cuts_iter_.second.second < 0)) continue; // don't care condition
    if((Pixel_Mux_map_[Ntracks_cuts_iter_.first] < Ntracks_cuts_iter_.second.first) ||
       (Pixel_Mux_map_[Ntracks_cuts_iter_.first] > Ntracks_cuts_iter_.second.second)) //condition violated
      Sector_TBA[Ntracks_cuts_iter_.first.first] = false;
  } 
	
  if(!(Sector_TBA[0] || Sector_TBA[1])) return;

  
  ////////////////////////////////////////////////////////////////
  //
  //		EXTRACT Dimoand detector info
  //
  ///////////////////////////////////////////////////////////////// 
  DiamondDet.ExtractData(iEvent);
 

  ////////////////////////////////////////////////////////////////
  //
  //		control over PCL calibration quality
  //
  /////////////////////////////////////////////////////////////////  
  for (const auto& recHits : *timingRecHit){ //rechits = array of hits in one channel
    const CTPPSDiamondDetId detId(recHits.detId());
	  if(!(Sector_TBA[detId.arm()])) continue;
	
	  // Perform channel histogram
    for (const auto& recHit : recHits){ //rechit
		  //std::cout << "Hits in channel " << detId.channel() << ", plane " << detId.plane() << ", rp " << detId.rp()<< ", station " << detId.station()<< 
	    //", arm " << detId.arm() << std::endl;

		  if (((recHit.ootIndex() !=0) && validOOT != -1) || recHit.multipleHits()) 
        continue;
	
      //T,TOT and OOT for all hits, important for monitoring the calibration	
		  histos.t[detId.rawId()]-> Fill(recHit.time());
		  histos.tot[detId.rawId()]-> Fill(recHit.toT());

      // T,TOT and OOT complete hits (T and TOT available), important for monitoring the calibration
      if(DiamondDet.PadActive(detId.arm(), detId.plane(),detId.channel())){
        histos.valid_tot[detId.rawId()]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()));
        histos.t_vs_tot[detId.rawId()]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()));
        histos.valid_t[detId.rawId()]-> Fill(DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()));
      }
    }
  }
}


void DiamondTimingWorker::bookHistograms(DQMStore::IBooker& iBooker,
                                           edm::Run const& run,
                                           edm::EventSetup const& iSetup) {
  
  std::string ch_name, ch_path;
  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (!CTPPSDiamondDetId::check(it->first))
      continue;
    
    const CTPPSDiamondDetId detid(it->first);

    // if(detid.station() != 1)
    //   continue;

    detid.channelName(ch_name);
    detid.channelName(ch_path, CTPPSDiamondDetId::nPath);

    iBooker.setCurrentFolder(ch_path);
    
    histos.t[detid.rawId()] = iBooker.book1D("t_" + ch_name, ch_name + ";t (ns);Entries", 1200, -60., 60.);
    histos.valid_t[detid.rawId()] = iBooker.book1D("valid_t_" + ch_name, ch_name + ";t (ns);Entries", 1200, -60., 60.);
    histos.tot[detid.rawId()] = iBooker.book1D("tot_" + ch_name, ch_name + ";ToT (ns);Entries", 100, -20., 20.);
    histos.valid_tot[detid.rawId()] = iBooker.book1D("valid_tot_" + ch_name, ch_name + ";ToT (ns);Entries", 100, -20., 20.);
    histos.t_vs_tot[detid.rawId()] =
        iBooker.book2D("t_vs_tot_" + ch_name, ch_name + ";ToT (ns);t (ns)", 240, 0., 60., 450, -20., 25.);
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingWorker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setAllowAnything();
  // desc.add<std::string>("folder", "MY_FOLDER");
  // descriptions.add("DiamondTimingWorker", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingWorker);
