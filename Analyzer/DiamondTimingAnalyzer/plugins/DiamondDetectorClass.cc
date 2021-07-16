
#include "DiamondDetectorClass.h"



namespace pt = boost::property_tree;

DiamondDetectorClass::DiamondDetectorClass(const edm::ParameterSet& iConfig, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_input, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_input)
:

tokenRecHit_          ( tokenRecHit_input),
tokenLocalTrack_          ( tokenLocalTrack_input),
saturationV_ (2),
Calib_input_path_ (iConfig.getParameter< std::string >( "tagCalibrationFile" )),
valid_OOT_ (iConfig.getParameter< int >( "tagValidOOT" ))
{
	CalibFunc = new TF1("myFermi","[0]/(exp((x-[1])/[2])+1)+[3]",0,30);
	
	SPC_map_= this-> GetDefaultCalibrations();
	this-> ExtractCalibrations();

}

void DiamondDetectorClass::ExtractCalibrations()
{
	

	// Create a root
	pt::ptree root;
	// Load the json file in this ptree
	pt::read_json(Calib_input_path_.c_str(), root);

	std::cout << "reading calibration from "<< Calib_input_path_.c_str()<< std::endl;
	
	SPC_map_.clear();
   
   // pt::write_json(std::cout, root);	
	for (const auto& sec_iter : root.get_child("Parameters.Sectors") )
	{
		pt::ptree sec_node = sec_iter.second;
		int sec_id = sec_node.get<int>("sector");
		for (const auto& pl_iter : sec_node.get_child("Planes") )
		{
			pt::ptree pl_node = pl_iter.second;
			int pl_id = pl_node.get<int>("plane");
			for (const auto& ch_iter : pl_node.get_child("Channels") )
			{
				pt::ptree ch_node = ch_iter.second;
				int ch_id = ch_node.get<int>("channel");
				std::vector<double> par_vec;
				for (const auto& params_iter : ch_node.get_child("param") )
				{
					par_vec.push_back(params_iter.second.get_value<double>());
				}
				Calib_par ch_calibration(ch_node.get<double>("time_offset", 0.0), ch_node.get<double>("time_precision"), par_vec);
				SPC_map_[ChannelKey(sec_id,pl_id,ch_id)]=ch_calibration;
				std::cout << "Calib found (sec "<< sec_id <<" pl "<< pl_id << " ch "<< ch_id << " : " << " off "<< ch_node.get<double>("time_offset", 0.0)  
							<< " pre "<< ch_node.get<double>("time_precision")  << " p0 "<< par_vec[0]<< " p1 " << par_vec[1]
							<< " p2 "<< par_vec[2]<< " p3 "<< par_vec[3] << std::endl;
			}
		}
	}
	

}



DiamondDetectorClass::~ DiamondDetectorClass()
{
}


void DiamondDetectorClass::ExtractData(const edm::Event& iEvent)
{
	
	
	edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > timingRecHit;
	edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > LocalTracks;
	iEvent.getByToken( tokenRecHit_, timingRecHit );
	iEvent.getByToken( tokenLocalTrack_, LocalTracks );
	
	
	
	RecHit_map_.clear();
	RecHit_map_.clear();
	LocalTrack_map_.clear();
	LocalTrack_map_.clear();
	saturationV_[0]=0;
	saturationV_[1]=0;
	Mux_map_.clear();
	Mux_validT_map_.clear();
	Mux_inTrack_map_.clear();
	
	
	
	
    //extract timing  tracks
	
	for (const auto& locTracks : *LocalTracks) //rechits = array of hits in one channel
	{
		for (const auto& locTrack : locTracks)
		{
			if((valid_OOT_==2) && (locTrack.ootIndex() != valid_OOT_) && (locTrack.ootIndex() != valid_OOT_+1)) continue;
			if (valid_OOT_!=2 && (locTrack.ootIndex() != valid_OOT_) &&  valid_OOT_!=-1 ) continue;
			
			LocalTrack_map_.insert(std::make_pair(locTrack,std::vector<std::pair<ChannelKey,CTPPSDiamondRecHit>>()));
				
			
		}
	}
				//	std::cout <<  " track found in 45 = " << LocalTrack_mapV_[0].size() << std::endl;
  
    // extract reco hit
  
	for (const auto& recHits : *timingRecHit) //rechits = array of hits in one channel
	{
		const CTPPSDiamondDetId detId( recHits.detId() );
		
		// retrieve and order all events in map. 
		
		for (const auto& recHit : recHits ) //rechit
		{
			//std::cout << "Hits in channel " << detId.channel() << ", plane " << detId.plane() << ", rp " << detId.rp()<< ", station " << detId.station()<< ", arm " << detId.arm() << std::endl;
					    
			
			if (((recHit.ootIndex() != (int)(SPC_map_[ChannelKey(detId.arm(),detId.plane(),detId.channel())].offset/25))  &&  valid_OOT_!=-1) ||  recHit.multipleHits()) continue;
			
			
			Mux_map_[std::make_pair(detId.arm(),detId.plane())]++;
			
			//put in hit map ("select hit with valid leading time and TOT")
			
			if(recHit.time()!=0.0 && recHit.toT()> 0.0) 
			{
				RecHit_map_[ChannelKey(detId.arm(),detId.plane(),detId.channel())].push_back(recHit);
				Mux_validT_map_[std::make_pair(detId.arm(),detId.plane())]++;
				
				bool counted = false;
			
				for (auto& locTrack_mapIter : LocalTrack_map_) 
				{
					if (locTrack_mapIter.first.containsHit(recHit, 0.1))
					{
						locTrack_mapIter.second.push_back(std::make_pair(ChannelKey(detId.arm(),detId.plane(),detId.channel()),recHit));
						if (!counted)
						{
							Mux_inTrack_map_[std::make_pair(detId.arm(),detId.plane())]++;
							counted=true;
						}
						//std::cout << "hit of sector "<<detId.arm() << " assigned to track with z= " << locTrack_mapIter.first.getZ0() << std::endl;
					}
				}
			
			}	
		
		}
	}
				//	std::cout <<  " track found in 45 = " << LocalTrack_mapV_[0].size() << std::endl;
  
	//check for saturation
  
	for (int sector_id=0; sector_id <2; sector_id++)
	{		
		for  (int ch_number=0; ch_number < CHANNELS_X_PLANE; ch_number++)
		{
			if (RecHit_map_.find(ChannelKey(sector_id,PLANE_2_ID,ch_number)) != RecHit_map_.end() && RecHit_map_.find(ChannelKey(sector_id,PLANE_3_ID,ch_number)) != RecHit_map_.end())
			{	 
				if (RecHit_map_[ChannelKey(sector_id,PLANE_2_ID,ch_number)].at(0).toT() > 15 && RecHit_map_[ChannelKey(sector_id,PLANE_3_ID,ch_number)].at(0).toT() > 15) // double saturated
				{	
					saturationV_[sector_id] = 1;
					//std::cout << "saturated!!" << std::endl;
				}
					
			} 
		}
	
	//std::cout << RecHit_mapV_[sector_id].size() << std::endl;
	
	}   
}

double DiamondDetectorClass::GetTime_SPC(int sector, int plane, int channel)
{
	double Time_raw = RecHit_map_[ChannelKey(sector,plane,channel)].at(0).time();
	double ToT = RecHit_map_[ChannelKey(sector,plane,channel)].at(0).toT();
	
	//TF1 *myfermi = new TF1("myFermi","[0]/(exp((x-[1])/[2])+1)+[3]",0,30);
	
	double P0 = SPC_map_[ChannelKey(sector,plane,channel)].params[0];
	double P1 = SPC_map_[ChannelKey(sector,plane,channel)].params[1];
	double P2 = SPC_map_[ChannelKey(sector,plane,channel)].params[2];
	double P3 = SPC_map_[ChannelKey(sector,plane,channel)].params[3];
	
	CalibFunc -> SetParameters(P0,P1,P2,P3);
	
	double Time_SPC = Time_raw - ( CalibFunc-> Eval(ToT) ); 
	
	//double Time_SPC = Time_raw - (C0 + C1*ToT + C2*ToT*ToT + C3*ToT*ToT*ToT);
 
	return Time_SPC;
	
	
}

int DiamondDetectorClass::GetSpread(int sector, int plane)
{
	
	int min=13, max=-1;
	for (const auto& rechit : RecHit_map_)
	{
		int rh_plane = rechit.first.plane;
		int rh_sector = rechit.first.sector;
		if ((rh_plane == plane) && (rh_sector == sector))
		{
			int rh_channel = Ch_position[rechit.first.channel]; 
			if (min > rh_channel) min = rh_channel;
			if (max < rh_channel) max = rh_channel;
		}
	}
	
	
	if ((max-min)>=0) return max-min+1;
	else return 0;
	
	
}

int DiamondDetectorClass::GetTrackSpread(int sector, int plane)
{
	
	int min=13, max=-1;
	for (const auto& locTrack_mapIter : LocalTrack_map_) //loop on track	
	{	
		for (const auto&  HitInfo_Hit_pair_iter : locTrack_mapIter.second)
		{
			int rh_plane = HitInfo_Hit_pair_iter.first.plane;
			int rh_sector = HitInfo_Hit_pair_iter.first.sector;
			if ((rh_plane == plane) && (rh_sector == sector))
			{
				int rh_channel = Ch_position[HitInfo_Hit_pair_iter.first.channel]; 
				if (min > rh_channel) min = rh_channel;
				if (max < rh_channel) max = rh_channel;
			}
		}
	}
	
	
	if ((max-min)>=0) return max-min+1;
	else return 0;
	
	
}


int DiamondDetectorClass::GetTrackMuxInSector(int sector)
{
	int mux=0;
	for (const auto& locTrack_mapIter : LocalTrack_map_) //loop on track	
	{	
		if (((locTrack_mapIter.first.z0() > 0.0 ) && sector==SECTOR_45_ID) ||
		    ((locTrack_mapIter.first.z0() < 0.0 ) && sector==SECTOR_56_ID))
			mux++;		
	}

	return mux;	
	
}


std::map < ChannelKey, Calib_par >  DiamondDetectorClass::GetDefaultCalibrations ()
        {
          std::map< ChannelKey,Calib_par  > m;
		  m[ChannelKey(0,0,0) ] = Calib_par(0.0,0.400,  std::vector<double> {0.0,0.0,0.0,0.0});
		  m[ChannelKey(0,0,1) ] = Calib_par(0.0,0.400,  std::vector<double> {0.0,0.0,0.0,0.0});
		  m[ChannelKey(0,0,2) ] = Calib_par(0.0,0.200,  std::vector<double> {2.006,14.34,1.505,4.63});
		  m[ChannelKey(0,0,3) ] = Calib_par(0.0,0.200,  std::vector<double> {1.127,11.55,0.6548,5.272});
		  m[ChannelKey(0,0,4) ] = Calib_par(0.0,0.200,  std::vector<double> {1.191,11.86,0.7645,5.306});
		  m[ChannelKey(0,0,5) ] = Calib_par(0.0,0.200,  std::vector<double> {1.557,14.28,1.022,4.963});
		  m[ChannelKey(0,0,6) ] = Calib_par(0.0,0.200,  std::vector<double> {2.067,14.33,1.362,4.674});
		  m[ChannelKey(0,0,7) ] = Calib_par(0.0,0.200,  std::vector<double> {1.178,12.04,0.6804,5.165});
		  m[ChannelKey(0,0,8) ] = Calib_par(0.0,0.200,  std::vector<double> {0.8463,11.67,0.5013,5.311});
		  m[ChannelKey(0,0,9) ] = Calib_par(0.0,0.200,  std::vector<double> {1.027,13.92,0.8086,4.901});
		  m[ChannelKey(0,0,10)] = Calib_par(0.0,0.200,  std::vector<double> {1.649,11.59,1.155,5.103});
		  m[ChannelKey(0,0,11)] = Calib_par(0.0,0.200,  std::vector<double> {2.308,11.02,0.8602,5.158});
		  m[ChannelKey(0,1,0) ] = Calib_par(0.0,0.200,  std::vector<double> {1.915,12.3,1.26,5.161});
		  m[ChannelKey(0,1,1) ] = Calib_par(0.0,0.200,  std::vector<double> {1.03,11.87,0.5981,5.621});
		  m[ChannelKey(0,1,2) ] = Calib_par(0.0,0.200,  std::vector<double> {0.9698,11.79,0.5261,5.741});
		  m[ChannelKey(0,1,3) ] = Calib_par(0.0,0.200,  std::vector<double> {1.57,11.72,0.7466,4.77});
		  m[ChannelKey(0,1,4) ] = Calib_par(0.0,0.200,  std::vector<double> {0.9819,12.21,0.5467,4.93});
		  m[ChannelKey(0,1,5) ] = Calib_par(0.0,0.200,  std::vector<double> {2.326,11.65,0.8615,5.738});
		  m[ChannelKey(0,1,6) ] = Calib_par(0.0,0.200,  std::vector<double> {1.1815,12.65,0.6459,5.525});
		  m[ChannelKey(0,1,7) ] = Calib_par(0.0,0.200,  std::vector<double> {1.523,12.42,0.8805,4.582});
		  m[ChannelKey(0,1,8) ] = Calib_par(0.0,0.200,  std::vector<double> {2.455,11.17,0.8015,4.774});
		  m[ChannelKey(0,1,9) ] = Calib_par(0.0,0.200,  std::vector<double> {1.415,11.61,0.7076,5.533});
		  m[ChannelKey(0,1,10)] = Calib_par(0.0,0.200,  std::vector<double> {12950,-1.953,1.384,4.817});
		  m[ChannelKey(0,1,11)] = Calib_par(0.0,0.150,  std::vector<double> {1.512,11.76,1.039,4.592});
		  m[ChannelKey(0,2,0) ] = Calib_par(0.0,0.150,  std::vector<double> {2.34,11.23,0.8874,4.048});
		  m[ChannelKey(0,2,1) ] = Calib_par(0.0,0.150,  std::vector<double> {2.717,10.23,1.356,4.343});
		  m[ChannelKey(0,2,2) ] = Calib_par(0.0,0.150,  std::vector<double> {4.953,9.694,0.8014,4.372});
		  m[ChannelKey(0,2,3) ] = Calib_par(0.0,0.150,  std::vector<double> {1.807,11.57,0.6542,4.783});
		  m[ChannelKey(0,2,4) ] = Calib_par(0.0,0.150,  std::vector<double> {1.823,11.97,1.143,4.638});
		  m[ChannelKey(0,2,5) ] = Calib_par(0.0,0.150,  std::vector<double> {1.987,11.47,0.6343,4.526});
		  m[ChannelKey(0,2,6) ] = Calib_par(0.0,0.150,  std::vector<double> {1.555,12.04,0.6252,4.535});
		  m[ChannelKey(0,2,7) ] = Calib_par(0.0,0.150,  std::vector<double> {3.029,11.03,1.408,4.506});
		  m[ChannelKey(0,2,8) ] = Calib_par(0.0,0.150,  std::vector<double> {1.282,12.04,0.6962,4.713});
		  m[ChannelKey(0,2,9) ] = Calib_par(0.0,0.150,  std::vector<double> {2.091,10.8,0.6909,4.295});
		  m[ChannelKey(0,2,10)] = Calib_par(0.0,0.150,  std::vector<double> {1.706,11.44,1.209,4.801});
		  m[ChannelKey(0,2,11)] = Calib_par(0.0,0.150,  std::vector<double> {2.433,11.73,1.38,4.326});
		  m[ChannelKey(0,3,0) ] = Calib_par(0.0,0.150,  std::vector<double> {2.612,11.8,1.677,4.24});
		  m[ChannelKey(0,3,1) ] = Calib_par(0.0,0.150,  std::vector<double> {2.253,11.38,1.792,4.782});
		  m[ChannelKey(0,3,2) ] = Calib_par(0.0,0.150,  std::vector<double> {1.425,11.3,0.8377,4.982});
		  m[ChannelKey(0,3,3) ] = Calib_par(0.0,0.150,  std::vector<double> {1.557,11.37,0.6042,4.456});
		  m[ChannelKey(0,3,4) ] = Calib_par(0.0,0.150,  std::vector<double> {0.9531,11.64,0.6303,4.378});
		  m[ChannelKey(0,3,5) ] = Calib_par(0.0,0.150,  std::vector<double> {2.481,11.12,0.9857,5.059});
		  m[ChannelKey(0,3,6) ] = Calib_par(0.0,0.150,  std::vector<double> {2.889,11.6,0.9576,4.846});
		  m[ChannelKey(0,3,7) ] = Calib_par(0.0,0.150,  std::vector<double> {0.9654,11.7,0.6441,4.266});
		  m[ChannelKey(0,3,8) ] = Calib_par(0.0,0.150,  std::vector<double> {2.077,11.16,0.805,4.314});
		  m[ChannelKey(0,3,9) ] = Calib_par(0.0,0.150,  std::vector<double> {1.228,11.38,0.7877,5.022});
		  m[ChannelKey(0,3,10)] = Calib_par(0.0,0.150,  std::vector<double> {1.933,11.41,1.387,4.281});
		  m[ChannelKey(0,3,11)] = Calib_par(0.0,0.150,  std::vector<double> {2.013,12.17,1.365,3.863});         
		  m[ChannelKey(1,0,0) ] = Calib_par(0.0,0.200,  std::vector<double> {1.392,12.19,0.4392,5.645});
		  m[ChannelKey(1,0,1) ] = Calib_par(0.0,0.200,  std::vector<double> {2.423,11.74,1.435,5.555});
		  m[ChannelKey(1,0,2) ] = Calib_par(0.0,0.200,  std::vector<double> {1.751,11.89,1.006,5.747});
		  m[ChannelKey(1,0,3) ] = Calib_par(0.0,0.200,  std::vector<double> {1.14,12.19,0.6755,6.112});
		  m[ChannelKey(1,0,4) ] = Calib_par(0.0,0.200,  std::vector<double> {1.057,12.43,0.6998,6.249});
		  m[ChannelKey(1,0,5) ] = Calib_par(0.0,0.200,  std::vector<double> {2.052,11.96,0.9994,5.843});
		  m[ChannelKey(1,0,6) ] = Calib_par(0.0,0.200,  std::vector<double> {2.763,11.54,1.02,5.777});
		  m[ChannelKey(1,0,7) ] = Calib_par(0.0,0.200,  std::vector<double> {1.262,12.56,0.7395,6.036});
		  m[ChannelKey(1,0,8) ] = Calib_par(0.0,0.200,  std::vector<double> {0.9912,12.43,0.5865,6.053});
		  m[ChannelKey(1,0,9) ] = Calib_par(0.0,0.200,  std::vector<double> {1.115,11.83,0.7117,5.765});
		  m[ChannelKey(1,0,10)] = Calib_par(0.0,0.200,  std::vector<double> {1.354,12.06,0.9533,6.269});
		  m[ChannelKey(1,0,11)] = Calib_par(0.0,0.200,  std::vector<double> {1.644,12.38,0.6117,6.051});
		  m[ChannelKey(1,1,0) ] = Calib_par(0.0,0.200,  std::vector<double> {2.511,11.97,1.671,5.887});
		  m[ChannelKey(1,1,1) ] = Calib_par(0.0,0.200,  std::vector<double> {1.196,12.28,0.8028,6.532});
		  m[ChannelKey(1,1,2) ] = Calib_par(0.0,0.200,  std::vector<double> {0.9922,12.59,0.5932,6.385});
		  m[ChannelKey(1,1,3) ] = Calib_par(0.0,0.200,  std::vector<double> {1.399,11.81,0.6517,5.725});
		  m[ChannelKey(1,1,4) ] = Calib_par(0.0,0.200,  std::vector<double> {1.082,12.1,0.5651,5.861});
		  m[ChannelKey(1,1,5) ] = Calib_par(0.0,0.200,  std::vector<double> {1.698,13.32,0.593,6.537});
		  m[ChannelKey(1,1,6) ] = Calib_par(0.0,0.200,  std::vector<double> {1.572,13.49,0.5928,6.291});
		  m[ChannelKey(1,1,7) ] = Calib_par(0.0,0.200,  std::vector<double> {1.267,12.33,0.7168,5.698});
		  m[ChannelKey(1,1,8) ] = Calib_par(0.0,0.200,  std::vector<double> {1.375,11.89,0.6036,5.688});
		  m[ChannelKey(1,1,9) ] = Calib_par(0.0,0.200,  std::vector<double> {1.064,12.61,0.5809,6.294});
		  m[ChannelKey(1,1,10)] = Calib_par(0.0,0.200,  std::vector<double> {1.02,11.74,0.6375,5.83});
		  m[ChannelKey(1,1,11)] = Calib_par(0.0,0.200,  std::vector<double> {1.465,11.5,0.9603,5.735});
		  m[ChannelKey(1,2,0) ] = Calib_par(0.0,0.150,  std::vector<double> {2.203,12.73,1.584,3.654});
		  m[ChannelKey(1,2,1) ] = Calib_par(0.0,0.150,  std::vector<double> {1.028,12.75,0.7025,4.12});
		  m[ChannelKey(1,2,2) ] = Calib_par(0.0,0.150,  std::vector<double> {1.686,12.08,0.6042,4.085});
		  m[ChannelKey(1,2,3) ] = Calib_par(0.0,0.150,  std::vector<double> {1.636,11.42,0.878,4.235});
		  m[ChannelKey(1,2,4) ] = Calib_par(0.0,0.150,  std::vector<double> {1.556,11.9,0.7281,4.336});
		  m[ChannelKey(1,2,5) ] = Calib_par(0.0,0.150,  std::vector<double> {1.825,12.91,0.8837,4.246});
		  m[ChannelKey(1,2,6) ] = Calib_par(0.0,0.150,  std::vector<double> {1.528,13.29,0.7223,4.293});
		  m[ChannelKey(1,2,7) ] = Calib_par(0.0,0.150,  std::vector<double> {1.426,12,0.7137,4.331});
		  m[ChannelKey(1,2,8) ] = Calib_par(0.0,0.150,  std::vector<double> {1.38,11.77,0.6664,4.22});
		  m[ChannelKey(1,2,9) ] = Calib_par(0.0,0.150,  std::vector<double> {2.399,11.49,1.2,3.85});
		  m[ChannelKey(1,2,10)] = Calib_par(0.0,0.150,  std::vector<double> {1.068,11.95,0.7024,4.264});
		  m[ChannelKey(1,2,11)] = Calib_par(0.0,0.150,  std::vector<double> {1.517,11.81,1.075,4.122});
		  m[ChannelKey(1,3,0) ] = Calib_par(0.0,0.400,  std::vector<double> {0.0,0.0,0.0,0.0});
		  m[ChannelKey(1,3,1) ] = Calib_par(0.0,0.150,  std::vector<double> {1.355,12.92,1.07,4.287});
		  m[ChannelKey(1,3,2) ] = Calib_par(0.0,0.150,  std::vector<double> {1.294,12.3,0.7493,4.428});
		  m[ChannelKey(1,3,3) ] = Calib_par(0.0,0.150,  std::vector<double> {1.037,12.21,0.6056,4.223});
		  m[ChannelKey(1,3,4) ] = Calib_par(0.0,0.150,  std::vector<double> {1.376,12.38,0.9060,4.241});
		  m[ChannelKey(1,3,5) ] = Calib_par(0.0,0.150,  std::vector<double> {1.742,12.81,0.8493,4.551});
		  m[ChannelKey(1,3,6) ] = Calib_par(0.0,0.150,  std::vector<double> {1.131,13.61,0.65,4.537});
		  m[ChannelKey(1,3,7) ] = Calib_par(0.0,0.150,  std::vector<double> {2.029,11.88,1.15,4.054});
		  m[ChannelKey(1,3,8) ] = Calib_par(0.0,0.150,  std::vector<double> {1.237,12.42,0.844,4.055});
		  m[ChannelKey(1,3,9) ] = Calib_par(0.0,0.150,  std::vector<double> {1.163,12.63,0.7625,4.278});
		  m[ChannelKey(1,3,10)] = Calib_par(0.0,0.150,  std::vector<double> {1.157,12.66,0.8589,3.944});
		  m[ChannelKey(1,3,11)] = Calib_par(0.0,0.150,  std::vector<double> {1.571,12.81,1.208,3.787});
          return m;
}

void DiamondDetectorClass::DumpEvent()
{
	for (const auto& recHits : RecHit_map_) //rechits = array of hits in one channel
	{
		std::cout << "Hit in sector " << recHits.first.sector << ", plane " << recHits.first.plane << ", channel " << recHits.first.channel << 
					 " : T " << recHits.second.at(0).time() << ", TOT " << recHits.second.at(0).toT() << std::endl;
	}
}

