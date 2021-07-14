#ifndef BOSS__DIAMONDDETECTORCLASS__LIBGUARD
#define BOSS__DIAMONDDETECTORCLASS__LIBGUARD


// system include files
#include <memory> 

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"



#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include "TF1.h"
 
#include <map>
#include <math.h>


		
		
struct ChannelKey
	{
		int sector, plane, channel;
	
		bool operator < (const ChannelKey &rhs) const {
			return std::tie(this->sector, this->plane, this->channel) < std::tie(rhs.sector, rhs.plane, rhs.channel);
		}
	
		ChannelKey(const int sector=-1, const int plane=-1, const int channel=-1)
		: sector(sector), plane(plane), channel(channel)
		{};
	};
	
struct Calib_par
	{
		double offset, precision;
		std::vector<double> params;
	
		Calib_par(const double offset=0, const double precision=200, const std::vector<double> params={0.0,0.0,0.0,0.0})
		: offset(offset), precision(precision), params(params)
		{};
	
	};

	
	
const std::vector<int> Ch_position {0,2,4,6,8,10,11,9,7,5,3,1};
typedef std::map< CTPPSDiamondLocalTrack , std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit> > > LocalTrack_map_type;


		
class DiamondDetectorClass  {
	
	
   public:
      explicit DiamondDetectorClass(const edm::ParameterSet&, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> >, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > );
      ~ DiamondDetectorClass();
	  
	  void ExtractData(const edm::Event&);
	  
	  void ExtractCalibrations();
	  std::map < ChannelKey , Calib_par >  GetDefaultCalibrations();
	  
	  inline int GetMux(int sector, int plane) 			{return Mux_map_[std::make_pair(sector,plane)];}
	  inline int GetMuxValidT(int sector, int plane) 	{return Mux_validT_map_[std::make_pair(sector,plane)];}
	  inline int GetMuxInTrack(int sector, int plane) 	{return Mux_inTrack_map_[std::make_pair(sector,plane)];}
	  inline bool SecIsSaturated(int sector) 			{return saturationV_[sector]!=0;}
	  
	  inline bool PairActive(int sector, int planeA, int channelA, int planeB, int channelB)   			
		{ return (RecHit_map_.find(ChannelKey(sector,planeA,channelA)) != RecHit_map_.end() && RecHit_map_.find(ChannelKey(sector,planeB,channelB)) != RecHit_map_.end());}
		
	  inline bool PadActive(int sector, int plane, int channel)   			
		{ return (RecHit_map_.find(ChannelKey(sector,plane,channel)) != RecHit_map_.end());}
		
	  inline double GetTime(int sector, int plane, int channel)
		{return RecHit_map_[ChannelKey(sector,plane,channel)].at(0).getT();}
		
	  inline double GetPadPrecision(int sector, int plane, int channel)
		{return SPC_map_[ChannelKey(sector,plane,channel)].precision;}
		
	  inline double GetPadWeight(int sector, int plane, int channel)
		{return (pow(SPC_map_[ChannelKey(sector,plane,channel)].precision,-2));}
		
	  inline double GetPadOOT(int sector, int plane, int channel)
		{return RecHit_map_[ChannelKey(sector,plane,channel)].at(0).getOOTIndex();}
		
	  double GetTime_SPC(int, int , int );	
		
	  inline double GetToT(int sector, int plane, int channel)
		{return RecHit_map_[ChannelKey(sector,plane,channel)].at(0).getToT();}	
		
	  inline int GetMH(int sector, int plane, int channel)
		{return RecHit_map_[ChannelKey(sector,plane,channel)].at(0).getMultipleHits() == false ? 0 : 1;}	
		
		
	  inline std::vector<CTPPSDiamondRecHit>  GetPad(int sector, int plane, int channel) {return RecHit_map_[ChannelKey(sector,plane,channel)];}
	  inline std::map< ChannelKey, Calib_par> GetSPCMap() {return SPC_map_;}
		
	  
	  int GetSpread(int, int);
	  int GetTrackSpread(int, int);
	  
	  inline LocalTrack_map_type GetDiamondTrack_map (){return LocalTrack_map_;}	  
	  int  GetTrackMuxInSector(int sector);
	
	
	  void DumpEvent(void);

   private:
   
   	enum Sector_id_
	{
		SECTOR_45_ID,
		SECTOR_56_ID
    };
	  
	enum Plane_id_
	{
		PLANE_0_ID,
		PLANE_1_ID,
		PLANE_2_ID,
		PLANE_3_ID	
    };
	
    static const int CHANNELS_X_PLANE=12;
    static const int PLANES_X_DETECTOR=4;
    static const int MAX_SECTOR_NUMBER=2;
	
	std::map< ChannelKey, Calib_par>    SPC_map_;

	// ---------- objects to retrieve ---------------------------
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_;
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_;

	// ---------- data extracted---------------------------

	std::map< ChannelKey, std::vector<CTPPSDiamondRecHit> > RecHit_map_;
	LocalTrack_map_type  LocalTrack_map_;
	
	
	// ---------- mux map ---------------------------
	
	std::map< std::pair< int , int >, int> Mux_map_;   //arm, plane
	std::map< std::pair< int , int >, int> Mux_validT_map_;  //arm, plane
	std::map< std::pair< int , int >, int> Mux_inTrack_map_;  //arm, plane
	
	//  event status flags
	  
	std::vector<int> saturationV_;
	
	TF1 *CalibFunc;
	std::string Calib_input_path_;
	int valid_OOT_;
	
};





#endif
