// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingAnalyzer
// Class:      DiamondTimingAnalyzer
//
/**\class DiamondTimingAnalyzer DiamondTimingAnalyzer.cc Analyzer/DiamondTimingAnalyzer/plugins/DiamondTimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arkadiusz Wolk
//         Created:  Fri, 16 Jul 2021 08:40:49 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class DiamondTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit DiamondTimingAnalyzer(const edm::ParameterSet&);
  ~DiamondTimingAnalyzer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void initHistograms(const CTPPSDiamondDetId&);

  //
  // constants, enums and typedefs
  //
  typedef std::map<std::pair<int,int>, TH2F*> TH2F_map;
	typedef std::map<std::pair<int,int>, TH1F*> TH1F_map;

  static const int CHANNELS_X_PLANE  = 12;
  static const int PLANES_X_DETECTOR = 4;
  static const int MAX_SECTOR_NUMBER = 2;
	edm::Service<TFileService> fs_;

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

  // ---------- objects to retrieve ---------------------------
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondDigi> > tokenDigi_;
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_;
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_;
	edm::EDGetTokenT< edm::DetSetVector<CTPPSPixelLocalTrack> > tokenPixelLocalTrack_;

	// ---------- directories ---------------------------
	std::map< CTPPSDiamondDetId, TFileDirectory > mainDir_map_;
	std::map< CTPPSDiamondDetId, TFileDirectory > planeDir_map_;
	std::map< ChannelKey, TFileDirectory > channelDir_map_;
		
	std::vector<TFileDirectory> dir_cyl_V ;	
	std::vector < std::map<int, TFileDirectory > > dir_cyl_L2_3p_Res_V ;
	std::vector < std::map<int, TFileDirectory > > dir_cyl_L2Res_V ;
		

	// ---------- global histograms ---------------------------
	TH1F* RunNumber_Histo_;
		
	std::map< ChannelKey, TH1F*> T_Histo_map_;
	std::map< ChannelKey, TH1F*> TOT_Histo_map_;
	std::map< ChannelKey, TH1F*> ValidT_Histo_map_;
	std::map< ChannelKey, TH1F*> ValidTOT_Histo_map_;
	std::map< ChannelKey, TH2F*> TOTvsT_Histo_map_;
	std::map< ChannelKey, TH2F*> Pad_tomography_220_Hmap_;
	std::map< ChannelKey, TH2F*> Pad_tomography_210_Hmap_;
		
		
	// ---------- pixel mux map ---------------------------
	std::map< std::pair< int , int >, int> Pixel_Mux_map_; //arm, station
		

	// ---------- channel graphs ---------------------------	  
	std::vector<TH1F*>  Tracks_time_SPC_histo_V_;
	std::vector<TH1F*>  Tracks_resolution_histo_V_;
		
	std::map<ChannelKey, double> Resolution_L2_map_;
	std::map<ChannelKey, double> Resolution_L2_3p_map_;
	std::map<ChannelKey, TH1F*>  Resolution_L2_histo_map_;
	std::map<ChannelKey, TH1F*>  Resolution_L2_3p_histo_map_;
	std::map<ChannelKey, TH1F*>  TrkTime_L2_histo_map_;
	std::map<ChannelKey, TH1F*>  TrkTime_L2_3p_histo_map_;
	std::map<ChannelKey, TH1F*>  TrkExpectedRes_L2_histo_map_;
	std::map<ChannelKey, TH1F*>  TrkExpectedRes_L2_3p_histo_map_;
	std::map<std::pair<int,int>, TGraph*>  Resolution_L2_graph_map_;
	std::map<std::pair<int,int>, TGraph*>  Resolution_L2_3p_graph_map_;
	std::map<std::pair<int,int>, TGraph*>  ImprouvedRes_L2_graph_map_;
	std::map<std::pair<int,int>, TGraph*>  ImprouvedRes_L2_3p_graph_map_;
		
		
	std::vector<TH2F*>  Tracks_time_SPC_BX_Histo_V_;  //<sector>
	std::vector<TH2F*>  Tracks_time_SPC_LS_Histo_V_;  //<sector>	  
	std::map<int,TProfile*> Tracks_time_SPC_BX_Profile_map_;
	std::map<int,TProfile*> Tracks_time_SPC_LS_Profile_map_;
		
	TTree  *TCalibration_;
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
		
	// selection parameters
	std::map< std::pair< int , int >, std::pair< int , int > > Ntracks_cuts_map_; //arm, station ,, Lcut,Ucut
};
	

//
// constructors and destructor
//
DiamondTimingAnalyzer::DiamondTimingAnalyzer(const edm::ParameterSet& iConfig)
  : 
  tokenDigi_(consumes<edm::DetSetVector<CTPPSDiamondDigi>>(iConfig.getParameter<edm::InputTag>("tagDigi"))),
  tokenRecHit_(consumes<edm::DetSetVector<CTPPSDiamondRecHit>>(iConfig.getParameter<edm::InputTag>("tagRecHit"))),
  tokenLocalTrack_(consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagLocalTrack"))),
  tokenPixelLocalTrack_(consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>(iConfig.getParameter<edm::InputTag>("tagPixelLocalTrack"))),
  Tracks_time_SPC_histo_V_(2),
  Tracks_resolution_histo_V_(2),
  Tracks_time_SPC_BX_Histo_V_(2),
  Tracks_time_SPC_LS_Histo_V_(2), 
  DiamondDet(iConfig,tokenRecHit_,tokenLocalTrack_),
  valid_OOT_(iConfig.getParameter<int>("tagValidOOT")) {
  
  usesResource("TFileService"); 
  
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

void DiamondTimingAnalyzer::initHistograms(const CTPPSDiamondDetId& detId){
  ChannelKey recHitKey(detId.arm(),detId.plane(),detId.channel());
  //now do what ever initialization is needed
  if (mainDir_map_.find(detId.rpId()) == mainDir_map_.end()){
	  std::cout << "RP is new, creating directory and histos" << std::endl;
    std::string dirName;
    detId.rpName(dirName, CTPPSDiamondDetId::nPath);
    std::string rpName;
    detId.rpName(rpName, CTPPSDiamondDetId::nFull);

    // create directory for the detector, if not already done
    mainDir_map_[ detId.rpId() ] = fs_->mkdir( dirName );
	}

  if (planeDir_map_.find( detId.planeId() ) == planeDir_map_.end()){
	  std::cout << "Plane is new, creating directory and histos" << std::endl;
    std::string dirName;
    detId.planeName(dirName, CTPPSDiamondDetId::nPath);
    std::string plName;
    detId.planeName(plName, CTPPSDiamondDetId::nFull);

    // create directory for the plane, if not already done
    planeDir_map_[ detId.planeId() ] = fs_->mkdir( dirName );
  } 

  //create channels histograms
  if (channelDir_map_.find(recHitKey) == channelDir_map_.end()){
	  std::cout << "Channels is new, creating directory and histos" << std::endl;
    std::string dirName;
    detId.channelName(dirName, CTPPSDiamondDetId::nPath);
    std::string chName;
    detId.channelName(chName, CTPPSDiamondDetId::nFull);

    // create directory for the channel, if not already done
    channelDir_map_[recHitKey] = fs_->mkdir( dirName );
	  
	  
	  
	  // create all histograms
    std::string T_Histo_name(chName);
    T_Histo_name.insert(0, "T_Distribution_");
    T_Histo_map_[recHitKey] = channelDir_map_[recHitKey].make<TH1F>(T_Histo_name.c_str(), T_Histo_name.c_str(), 1200, -60, 60 );
    std::string ValidT_Histo_name(chName);
    ValidT_Histo_name.insert(0, "VALID_T_Distribution_");
    ValidT_Histo_map_[recHitKey] = channelDir_map_[recHitKey].make<TH1F>(ValidT_Histo_name.c_str(), ValidT_Histo_name.c_str(), 1200, -60, 60 );
    
	
    std::string TOT_Histo_name(chName);
    TOT_Histo_name.insert(0, "TOT_Distribution_");
    TOT_Histo_map_[recHitKey] = channelDir_map_[recHitKey].make<TH1F>(TOT_Histo_name.c_str(), TOT_Histo_name.c_str(), 100, -20, 20 );
    std::string ValidTOT_Histo_name(chName);
    ValidTOT_Histo_name.insert(0, "ValidTOT_Distribution_");
    ValidTOT_Histo_map_[recHitKey] = channelDir_map_[recHitKey].make<TH1F>(ValidTOT_Histo_name.c_str(), ValidTOT_Histo_name.c_str(), 100, -20, 20 );

    std::string TOTvsT_Histo_name(chName);
    TOTvsT_Histo_name.insert(0, "TvsTOT_Distribution_");
    TOTvsT_Histo_map_[recHitKey] = channelDir_map_[recHitKey].make<TH2F>(TOTvsT_Histo_name.c_str(), TOTvsT_Histo_name.c_str(), 240, 0, 60, 450, -20, 25 );

    std::string PadTomo220_Histo_name(chName);
    PadTomo220_Histo_name.insert(0, "Pad tomography 220");
    Pad_tomography_220_Hmap_[recHitKey] = channelDir_map_[recHitKey].make<TH2F>(PadTomo220_Histo_name.c_str(), PadTomo220_Histo_name.c_str(), 1000, -100, 100, 1000, -100, 100);

    std::string PadTomo210_Histo_name(chName);
    PadTomo210_Histo_name.insert(0, "Pad tomography 210");
    Pad_tomography_210_Hmap_[recHitKey] = channelDir_map_[recHitKey].make<TH2F>(PadTomo210_Histo_name.c_str(), PadTomo210_Histo_name.c_str(), 1000, -100, 100, 1000, -100, 100);
	} 
}

// ------------ method called for each event  ------------
void DiamondTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  //retrieve data
  edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > timingRecHit;
  edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > pixelLocalTrack;
  iEvent.getByToken( tokenRecHit_, timingRecHit );
  iEvent.getByToken( tokenPixelLocalTrack_, pixelLocalTrack );
  

  ////////////////////////////////////////////////////////////////
  //
  //		EXTRACT PIXELS TRACK NUMBER
  //      Will be used for sector independent event selection
  //
  ///////////////////////////////////////////////////////////////// 
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


  RunNumber_Histo_->Fill(iEvent.id().run());

  
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
	
	  ChannelKey recHitKey(detId.arm(),detId.plane(),detId.channel());
	
    if(channelDir_map_.find(recHitKey) == channelDir_map_.end() && recHits.size() > 0){  
		  std::cout << "Found new channel with data! Creating Histograms..." << std::endl;
		  initHistograms( detId );
	  }
	
	  // Perform channel histogram
    for (const auto& recHit : recHits){ //rechit
		  //std::cout << "Hits in channel " << detId.channel() << ", plane " << detId.plane() << ", rp " << detId.rp()<< ", station " << detId.station()<< 
	    //", arm " << detId.arm() << std::endl;

		  if (((recHit.ootIndex() !=0 ) &&  valid_OOT_!=-1) ||  recHit.multipleHits()) continue;
	
      //T,TOT and OOT for all hits, important for monitoring the calibration	
		  T_Histo_map_[recHitKey]-> Fill( recHit.time() );
		  TOT_Histo_map_[recHitKey]-> Fill( recHit.toT() );

      // T,TOT and OOT complete hits (T and TOT available), important for monitoring the calibration
      if (DiamondDet.PadActive(detId.arm(), detId.plane(),detId.channel())){
        ValidTOT_Histo_map_[recHitKey]-> Fill( DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()) );
        TOTvsT_Histo_map_[recHitKey]-> Fill(DiamondDet.GetToT(detId.arm(), detId.plane(),detId.channel()), DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()));
        ValidT_Histo_map_[recHitKey]-> Fill( DiamondDet.GetTime(detId.arm(), detId.plane(),detId.channel()) );
      }
    }
  }



  ////////////////////////////////////////////////////////////////
  //
  //		RESOLUTION STUDIES
  //
  ///////////////////////////////////////////////////////////////// 
	

  ////////////////////////////////////////////////////////////////
  //
  //		4 planes!!!!!
  //
  ///////////////////////////////////////////////////////////////// 

  for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map()){ // loop on predigested tracks
	  int sec_number = LocalTrack_mapIter.first.z0() > 0.0 ? SECTOR_45_ID : SECTOR_56_ID;
    if (!(Sector_TBA[sec_number])) continue;

    bool mark_tag=false;

    //std::cout << "check for mark" << std::endl;	
	  if(DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)==1 && 
		  DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)==1 && 
		  DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)==1 && 
		  DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)==1 &&
		  DiamondDet.GetTrackMuxInSector(sec_number)==1) mark_tag=true;
	
	  std::vector<ChannelKey>  hit_selected(4);
	 
	  //int NumPlanes=0;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)>0) NumPlanes++;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)>0) NumPlanes++;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)>0) NumPlanes++;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)>0) NumPlanes++;
	  //if (NumPlanes<4) continue;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)>1) continue;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)>1) continue;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)>1) continue;
	  //if (DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)>1) continue;
	  //if (DiamondDet.GetTrackMuxInSector(sec_number)!=1) continue;
    
	  std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>::const_iterator  hit_iter;

	  //std::cout << "check for size" << std::endl;		
	  if(LocalTrack_mapIter.second.size() == 0) continue;
	
	  double Track_time_SPC = 100.0;
	  double Track_precision_SPC = 100.0;

	  //std::cout << "start loop on hit" << std::endl;			
	  for (hit_iter=LocalTrack_mapIter.second.begin();hit_iter < LocalTrack_mapIter.second.end(); hit_iter++ ) // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> > 
	  {
	    //std::cout << "looping on hits" << std::endl;
		  double hit_time_SPC= DiamondDet.GetTime((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
		  double hit_prec_SPC= DiamondDet.GetPadPrecision((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
		  double hit_weig_SPC= DiamondDet.GetPadWeight((*hit_iter).first.sector,(*hit_iter).first.plane,(*hit_iter).first.channel);
		
		  if (mark_tag){
				//std::cout << "marking hits" << std::endl;
				hit_selected[(*hit_iter).first.plane] = (*hit_iter).first;  // save for resolution reco
		  }
		
		  if (hit_iter == LocalTrack_mapIter.second.begin()){
			
			  //	std::cout << "setting starting timing" << std::endl;
			  Track_time_SPC = hit_time_SPC;
			  Track_precision_SPC =  hit_prec_SPC;
		  }else{
				//std::cout << "updating  timing" << std::endl;
			  Track_time_SPC = (Track_time_SPC*pow(Track_precision_SPC,-2) + hit_time_SPC*hit_weig_SPC)/(pow(Track_precision_SPC,-2)+hit_weig_SPC);
			  Track_precision_SPC = pow((pow(Track_precision_SPC,-2)+hit_weig_SPC),-0.5);	
		  }
		
		  //std::cout<<"T="<<hit_time_SPC<<", P="<<hit_prec_SPC<<", W="<<hit_weig_SPC<<", Precision="<<Track_precision_SPC<<std::endl;
	  }
			
	  Tracks_time_SPC_histo_V_[ sec_number ]->Fill(Track_time_SPC);
	  Tracks_resolution_histo_V_[ sec_number ]->Fill(Track_precision_SPC);
	
	  Tracks_time_SPC_BX_Histo_V_[ sec_number ]->Fill(iEvent.bunchCrossing(), Track_time_SPC );  //<sector>
	  Tracks_time_SPC_LS_Histo_V_[ sec_number ]->Fill(iEvent.luminosityBlock(), Track_time_SPC );  //<sector>	  

	  if(mark_tag){
		  for (int pl_mark = 0 ; pl_mark < PLANES_X_DETECTOR; pl_mark++){
			  double Marked_track_time = 12.5;
			  double Marked_track_precision = 25.0;
			  double Marked_hit_time=0.0;
			  int Marked_hit_channel=-1;
			
			  for (int pl_loop = 0 ; pl_loop < PLANES_X_DETECTOR; pl_loop++){
				  if (pl_loop == pl_mark){
					  Marked_hit_time= DiamondDet.GetTime(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
					  Marked_hit_channel = hit_selected[pl_loop].channel;
					  continue;
				  }

				  double Others_hit_time= DiamondDet.GetTime(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
				  double Others_hit_prec= DiamondDet.GetPadPrecision(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
				  double Others_hit_weig= DiamondDet.GetPadWeight(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
				 
			    if (hit_iter == LocalTrack_mapIter.second.begin()){
					  Marked_track_time = Others_hit_time;
					  Marked_track_precision = Others_hit_prec;
          }else{	
					  Marked_track_time = (Marked_track_time*pow(Marked_track_precision,-2) + Others_hit_time*Others_hit_weig)/(pow(Marked_track_precision,-2)+Others_hit_weig);
					  Marked_track_precision = pow((pow(Marked_track_precision,-2)+Others_hit_weig),-0.5);	
				  }
			  }
			
			  double Marked_hit_difference = Marked_hit_time-Marked_track_time;
				
        // creta pair histos if do not exist
			  if (Resolution_L2_histo_map_.find(ChannelKey(sec_number,pl_mark,Marked_hit_channel))==Resolution_L2_histo_map_.end()){
				  std::string Histo_name = "L2_resolution_DT_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
											std::to_string(Marked_hit_channel);
				  Resolution_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
					  dir_cyl_L2Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );
					
				  Histo_name = "L2_partialTrk_Time_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
											std::to_string(Marked_hit_channel);
				  TrkTime_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
				  	dir_cyl_L2Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );		
					
				  Histo_name = "L2_partialTrk_ExpectedRes_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
											std::to_string(Marked_hit_channel);
				  TrkExpectedRes_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
					  dir_cyl_L2Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1000, 0, 2 );					
			  }
			
			  Resolution_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_hit_difference);
			  TrkTime_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_time);
			  TrkExpectedRes_L2_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_precision);
					
		  }
	  }
  }


  ////////////////////////////////////////////////////////////////
  //
  //		3 planes!!!!!
  //
  ///////////////////////////////////////////////////////////////// 

    
  for (const auto& LocalTrack_mapIter : DiamondDet.GetDiamondTrack_map()){ // loop on predigested tracks
    int sec_number = LocalTrack_mapIter.first.z0() > 0.0 ? SECTOR_45_ID : SECTOR_56_ID;
    if (!(Sector_TBA[sec_number])) continue;
    
    std::map<int,int> active_pl_map;
    if (sec_number==0){
      active_pl_map[0]=0;
      active_pl_map[1]=1;
      active_pl_map[2]=1;
      active_pl_map[3]=1;
    }else{
      active_pl_map[0]=1;
      active_pl_map[1]=1;
      active_pl_map[2]=1;
      active_pl_map[3]=0;
    }

    //std::cout << "check for mark" << std::endl;	
    if (!((DiamondDet.GetMuxInTrack(sec_number,PLANE_0_ID)==1 || active_pl_map[PLANE_0_ID]==0 )&& 
      (DiamondDet.GetMuxInTrack(sec_number,PLANE_1_ID)==1 || active_pl_map[PLANE_1_ID]==0 )&& 
      (DiamondDet.GetMuxInTrack(sec_number,PLANE_2_ID)==1 || active_pl_map[PLANE_2_ID]==0 )&& 
      (DiamondDet.GetMuxInTrack(sec_number,PLANE_3_ID)==1 || active_pl_map[PLANE_3_ID]==0 )&&
      DiamondDet.GetTrackMuxInSector(sec_number)==1)) continue;
    
    std::map<int,ChannelKey>  hit_selected;
    std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>::const_iterator  hit_iter;

    //std::cout << "check for size" << std::endl;		
    if(LocalTrack_mapIter.second.size() == 0) continue; //check if it no hit

    //std::cout << "start loop on hit" << std::endl;			
    
    // hits in track loop LocalTrack_mapIter.second = std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> > 
    for (hit_iter=LocalTrack_mapIter.second.begin();hit_iter < LocalTrack_mapIter.second.end(); hit_iter++ ){
      if (active_pl_map[(*hit_iter).first.plane] ==1)
        hit_selected[(*hit_iter).first.plane] = (*hit_iter).first;  // save for resolution reco
    }
    
    for (int pl_mark = 0 ; pl_mark < PLANES_X_DETECTOR; pl_mark++){
      if (active_pl_map[pl_mark] ==0) continue;
    
      double Marked_track_time = 12.5;
      double Marked_track_precision = 25.0;
      double Marked_hit_time=0.0;
      int Marked_hit_channel=-1;
      
      for (int pl_loop = 0 ; pl_loop < PLANES_X_DETECTOR; pl_loop++){
        if (active_pl_map[pl_loop] == 0) continue;
        if (pl_loop == pl_mark){
          Marked_hit_time= DiamondDet.GetTime(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
          Marked_hit_channel = hit_selected[pl_loop].channel;
          continue;
        }

        double Others_hit_time= DiamondDet.GetTime(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
        double Others_hit_prec= DiamondDet.GetPadPrecision(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
        double Others_hit_weig= DiamondDet.GetPadWeight(hit_selected[pl_loop].sector,hit_selected[pl_loop].plane,hit_selected[pl_loop].channel);
        
        if (hit_iter == LocalTrack_mapIter.second.begin()){
          
          Marked_track_time = Others_hit_time;
          Marked_track_precision = Others_hit_prec;
        }else{	
          Marked_track_time = (Marked_track_time*pow(Marked_track_precision,-2) + Others_hit_time*Others_hit_weig)/(pow(Marked_track_precision,-2)+Others_hit_weig);
          Marked_track_precision = pow((pow(Marked_track_precision,-2)+Others_hit_weig),-0.5);	
        }
      }
      
      double Marked_hit_difference = Marked_hit_time-Marked_track_time; 
        
      // creta pair histos if do not exist
      if (Resolution_L2_3p_histo_map_.find(ChannelKey(sec_number,pl_mark,Marked_hit_channel))==Resolution_L2_3p_histo_map_.end()){
        std::string Histo_name = "L2_3p_resolution_DT_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
                      std::to_string(Marked_hit_channel);
        Resolution_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
          dir_cyl_L2_3p_Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );
          
        Histo_name = "L2_3p_partialTrk_Time_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
                      std::to_string(Marked_hit_channel);
        TrkTime_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
          dir_cyl_L2_3p_Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1200, -60, 60 );		
          
        Histo_name = "L2_3p_partialTrk_ExpectedRes_sec_"+std::to_string(sec_number)+"_plane_"+std::to_string(pl_mark)+"_channel_"+
                      std::to_string(Marked_hit_channel);
        TrkExpectedRes_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)] = 
          dir_cyl_L2_3p_Res_V[sec_number][pl_mark].make<TH1F>(Histo_name.c_str(), Histo_name.c_str(), 1000, 0, 2 );					
      }
      
      Resolution_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_hit_difference);
      TrkTime_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_time);
      TrkExpectedRes_L2_3p_histo_map_[ChannelKey(sec_number,pl_mark,Marked_hit_channel)]-> Fill(Marked_track_precision);
          
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void DiamondTimingAnalyzer::beginJob() {
  dir_cyl_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr" ));
	dir_cyl_V.push_back(fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr" ));
	
	dir_cyl_L2Res_V.resize(2);	
	dir_cyl_L2Res_V[0][0]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2Res/plane0" );
	dir_cyl_L2Res_V[0][1]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2Res/plane1" );
	dir_cyl_L2Res_V[0][2]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2Res/plane2" );
	dir_cyl_L2Res_V[0][3]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2Res/plane3" );
	dir_cyl_L2Res_V[1][0]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2Res/plane0" );
	dir_cyl_L2Res_V[1][1]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2Res/plane1" );
	dir_cyl_L2Res_V[1][2]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2Res/plane2" );
	dir_cyl_L2Res_V[1][3]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2Res/plane3" );

	dir_cyl_L2_3p_Res_V.resize(2);	
	dir_cyl_L2_3p_Res_V[0][0]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2_3p_Res/plane0" );
	dir_cyl_L2_3p_Res_V[0][1]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2_3p_Res/plane1" );
	dir_cyl_L2_3p_Res_V[0][2]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2_3p_Res/plane2" );
	dir_cyl_L2_3p_Res_V[0][3]=fs_->mkdir( "CTPPS/TimingDiamond/sector 45/station 220cyl/cyl_hr/L2_3p_Res/plane3" );
	dir_cyl_L2_3p_Res_V[1][0]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2_3p_Res/plane0" );
	dir_cyl_L2_3p_Res_V[1][1]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2_3p_Res/plane1" );
	dir_cyl_L2_3p_Res_V[1][2]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2_3p_Res/plane2" );
	dir_cyl_L2_3p_Res_V[1][3]=fs_->mkdir( "CTPPS/TimingDiamond/sector 56/station 220cyl/cyl_hr/L2_3p_Res/plane3" );

	
	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++)
	{
		for(int pl_number=0; pl_number < PLANES_X_DETECTOR; pl_number++){
			Resolution_L2_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2Res_V[sec_number][pl_number].make<TGraph>(12);
			std::string Hname="L2 resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			Resolution_L2_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());
		
			Resolution_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2_3p_Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="L2_3p resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			Resolution_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());
		
			ImprouvedRes_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2_3p_Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="Improuved L2_3p resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			ImprouvedRes_L2_3p_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());	
		
			ImprouvedRes_L2_graph_map_[ std::make_pair(sec_number,pl_number)] = dir_cyl_L2Res_V[sec_number][pl_number].make<TGraph>(12);
			Hname="Improuved L2 resolution sector "+std::to_string(sec_number)+" plane "+std::to_string(pl_number);
			ImprouvedRes_L2_graph_map_[ std::make_pair(sec_number,pl_number)]  -> SetNameTitle(Hname.c_str(),Hname.c_str());			
		}
	}	
	
	for (int sec_number=0; sec_number < MAX_SECTOR_NUMBER; sec_number++){
	  std::string name="Timing track time SPC sector "+std::to_string(sec_number);
		Tracks_time_SPC_histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), 1200, -60, 60);			
	  
    name="Timing track resolution sector "+std::to_string(sec_number);
		Tracks_resolution_histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH1F>(name.c_str(), name.c_str(), 1000, 0, 1);

		name="Timing track time SPC Vs BX sector "+std::to_string(sec_number);
		Tracks_time_SPC_BX_Histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH2F>(name.c_str(), name.c_str(), 4000, 0, 4000, 500, -5 , 5 );
		
		name="Timing track time SPC Vs LS sector "+std::to_string(sec_number);
		Tracks_time_SPC_LS_Histo_V_[ sec_number ] = dir_cyl_V[sec_number].make<TH2F>(name.c_str(), name.c_str(), 4000, 0, 4000, 500, -5 , 5 );
	}
	
	TFileDirectory BaseDir = fs_->mkdir( "Calib/" );
	
	TCalibration_ = BaseDir.make<TTree>("CalibTree","rootuple");
	
	TCalibration_ -> Branch("Sector"  	, &(this->cal_sector_), "Sector/I");
	TCalibration_ -> Branch("Plane"		, &(this->cal_plane_), "Plane/I");
	TCalibration_ -> Branch("Channel"  	, &(this->cal_channel_), "Channel/I");
	TCalibration_ -> Branch("Par_0"		, &(this->cal_par_0_), "Par_0/D" );
	TCalibration_ -> Branch("Par_1"		, &(this->cal_par_1_), "Par_1/D" );
	TCalibration_ -> Branch("Par_2"		, &(this->cal_par_2_), "Par_2/D" );
	TCalibration_ -> Branch("Par_3"		, &(this->cal_par_3_), "Par_3/D" );
	TCalibration_ -> Branch("Res_L2"		, &(this->resolution_L2_), "Res_L2/D" );
	TCalibration_ -> Branch("Res_L2_3p"		, &(this->resolution_L2_3p_), "Res_L2_3p/D" );
	
	
	TFileDirectory  GlobalInfoDir = fs_->mkdir( "GlobalInfo/" );
	RunNumber_Histo_ = GlobalInfoDir.make<TH1F>("RunNumber", "RunNumber", 300000, 300000, 330000);
}

// ------------ method called once each job just after ending the event loop  ------------
void DiamondTimingAnalyzer::endJob() {
  std::cout << "Starting end job task" << std::endl;

  ///////////////////////////////////////////
  // deriving full track based L2 resolution	
  ///////////////////////////////////////////
	
	for (int sec_id=0; sec_id < MAX_SECTOR_NUMBER; sec_id ++){
		for (int pl_id=0; pl_id < PLANES_X_DETECTOR; pl_id ++){
			for (int ch_id=0; ch_id < CHANNELS_X_PLANE; ch_id ++){
				ChannelKey histo_key(sec_id,pl_id,ch_id);
				if (Resolution_L2_histo_map_.find(histo_key)!=Resolution_L2_histo_map_.end()){
					if (Resolution_L2_histo_map_[histo_key]->GetEntries() > 100){
						Resolution_L2_histo_map_[histo_key]->Fit("gaus","+Q","",-10,10);
			
						if (Resolution_L2_histo_map_[histo_key]->GetFunction("gaus")!= NULL){				
							double ResL2_mean = Resolution_L2_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(1);
							double ResL2_sigma = Resolution_L2_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(2);
							Resolution_L2_histo_map_[histo_key]->Fit("gaus","","",ResL2_mean-(2.2*ResL2_sigma),ResL2_mean+(2.2*ResL2_sigma));
							ResL2_sigma = Resolution_L2_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(2);
							double Exp_sigma =  TrkExpectedRes_L2_histo_map_[histo_key]-> GetMean();
							Resolution_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, ResL2_sigma);
							if (ResL2_sigma > Exp_sigma){
								ImprouvedRes_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2),0.5));
								Resolution_L2_map_[histo_key] = pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2),0.5);
							}else{
								ImprouvedRes_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0.05);
								Resolution_L2_map_[histo_key] = 0.05;
							}
						}else{
							Resolution_L2_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0);
							Resolution_L2_map_[histo_key] = 0.400; 
						}	
					}
				}
				if (Resolution_L2_3p_histo_map_.find(histo_key)!=Resolution_L2_3p_histo_map_.end()){
					if (Resolution_L2_3p_histo_map_[histo_key]->GetEntries() > 100){
						Resolution_L2_3p_histo_map_[histo_key]->Fit("gaus","+Q","",-10,10);
			
						if (Resolution_L2_3p_histo_map_[histo_key]->GetFunction("gaus")!= NULL){				
							double ResL2_3p_mean = Resolution_L2_3p_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(1);
							double ResL2_3p_sigma = Resolution_L2_3p_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(2);
							Resolution_L2_3p_histo_map_[histo_key]->Fit("gaus","","",ResL2_3p_mean-(2.2*ResL2_3p_sigma),ResL2_3p_mean+(2.2*ResL2_3p_sigma));
							ResL2_3p_sigma = Resolution_L2_3p_histo_map_[histo_key]->GetFunction("gaus")->GetParameter(2);
							double Exp_sigma =  TrkExpectedRes_L2_3p_histo_map_[histo_key]-> GetMean();
							Resolution_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, ResL2_3p_sigma);
							if (ResL2_3p_sigma > Exp_sigma){
								ImprouvedRes_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, pow(pow(ResL2_3p_sigma,2)-pow(Exp_sigma,2),0.5));
								Resolution_L2_3p_map_[histo_key] = pow(pow(ResL2_3p_sigma,2)-pow(Exp_sigma,2),0.5);
							}
							else{
								ImprouvedRes_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0.05);
								Resolution_L2_3p_map_[histo_key] = 0.05;
							}
						}
						else{
							std::cout << "WARNING: L2_3p resolution fit unseccessfull for sector  " << sec_id <<" plane " << pl_id <<" channel " << ch_id << std::endl;
							Resolution_L2_3p_graph_map_[std::make_pair(sec_id,pl_id)]->SetPoint(ch_id+1, ch_id, 0);
							Resolution_L2_3p_map_[histo_key] = 0.400; 
						}	
					}
				}
				
				// Merging  ch. calibration with updated resolution in output calibration tree				
				cal_sector_    = sec_id;
				cal_plane_  = pl_id;
				cal_channel_= ch_id;
				cal_par_0_  = DiamondDet.GetSPCMap()[histo_key].params[0];
				cal_par_1_  = DiamondDet.GetSPCMap()[histo_key].params[1];
				cal_par_2_  = DiamondDet.GetSPCMap()[histo_key].params[2];
				cal_par_3_  = DiamondDet.GetSPCMap()[histo_key].params[3];
				
				if (Resolution_L2_map_.find(histo_key) != Resolution_L2_map_.end() )
					resolution_L2_ = Resolution_L2_map_[histo_key];
				else
					resolution_L2_= 0.999;
				
				if (Resolution_L2_3p_map_.find(histo_key) != Resolution_L2_3p_map_.end() )
					resolution_L2_3p_ = Resolution_L2_3p_map_[histo_key];
				else
					resolution_L2_3p_= 0.999;
				
				TCalibration_ ->Fill();
			}						
		}
	}		

	//////////////////////////////////////////////////
  // profiles for SPC average trackTime Vs BX/LS	
  //////////////////////////////////////////////////
	for (int sec_id=0; sec_id < MAX_SECTOR_NUMBER; sec_id ++){
		Tracks_time_SPC_BX_Profile_map_[ sec_id ] = dir_cyl_V[ sec_id ].make<TProfile>(*Tracks_time_SPC_BX_Histo_V_[ sec_id ]->ProfileX("_BX_pfx",1,-1));
		Tracks_time_SPC_LS_Profile_map_[ sec_id ] = dir_cyl_V[ sec_id ].make<TProfile>(*Tracks_time_SPC_LS_Histo_V_[ sec_id ]->ProfileX("_LS_pfx",1,-1));
	}
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingAnalyzer);
