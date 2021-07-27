// -*- C++ -*-
//
// Package:    Analyzer/DiamondTimingHarvester
// Class:      DiamondTimingHarvester
//
/**\class DiamondTimingHarvester DiamondTimingHarvester.cc Analyzer/DiamondTimingHarvester/plugins/DiamondTimingHarvester.cc

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
#include <bitset>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DiamondDetectorClass.h"

#include "CondFormats/PPSObjects/interface/PPSTimingCalibration.h"
#include "CondFormats/DataRecord/interface/PPSTimingCalibrationRcd.h"

#include "DiamondTimingCalibration.h"
#include "JSONProducer.h"

//
// class declaration
//

class DiamondTimingHarvester : public DQMEDHarvester{
public:
  explicit DiamondTimingHarvester(const edm::ParameterSet&);
  ~DiamondTimingHarvester() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) override;
  void dqmEndRun(DQMStore::IBooker &iBooker,
                 DQMStore::IGetter &iGetter,
                 edm::Run const &iRun,
                 edm::EventSetup const &iSetup) override;

  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
  edm::ESGetToken<PPSTimingCalibration, PPSTimingCalibrationRcd> calibEsToken_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiamondTimingHarvester::DiamondTimingHarvester(const edm::ParameterSet& iConfig)
  : 
  geomEsToken_(esConsumes<edm::Transition::EndRun>()),
  calibEsToken_(esConsumes<edm::Transition::EndRun>())
  {
  // now do what ever initialization is needed
}

//
// member functions
//

// ------------ method called for each event  ------------

void DiamondTimingHarvester::dqmEndJob(DQMStore::IBooker &iBooker, DQMStore::IGetter &iGetter) {
  std::cout<<"######## EndJob ########"<<std::endl;
}

void DiamondTimingHarvester::dqmEndRun(DQMStore::IBooker &iBooker,
               DQMStore::IGetter &iGetter,
               edm::Run const &iRun,
               edm::EventSetup const &iSetup) {
  
  std::map<ChannelKey, double> Resolution_L2_map_;
  
  ///////////////////////////////////////////
  // deriving full track based L2 resolution	
  ///////////////////////////////////////////
  std::cout<<"####### EndRun #######"<<std::endl;

  std::string ch_name, ch_path;
  const auto& geom = iSetup.getData(geomEsToken_);
  const auto& calib = DiamondTimingCalibration(iSetup.getData(calibEsToken_));
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (!CTPPSDiamondDetId::check(it->first))
      continue;
    
    const CTPPSDiamondDetId detid(it->first);
    ChannelKey histo_key(detid);
    detid.channelName(ch_name);
    detid.channelName(ch_path, CTPPSDiamondDetId::nPath);

    auto* l2_res = iGetter.get(ch_path + "/" + "l2_res_" + ch_name);
    auto* expected_trk_time = iGetter.get(ch_path + "/" + "expected_trk_time_res_" + ch_name);
		if(l2_res->getEntries() > 100){
			l2_res->getTH1F()->Fit("gaus","+Q","",-10,10);
			
			if(l2_res->getTH1F()->GetFunction("gaus") != NULL){				
				double ResL2_mean = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(1);
				double ResL2_sigma = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(2);
				l2_res->getTH1F()->Fit("gaus","","",ResL2_mean-(2.2*ResL2_sigma),ResL2_mean+(2.2*ResL2_sigma));
				ResL2_sigma = l2_res->getTH1F()->GetFunction("gaus")->GetParameter(2);

				double Exp_sigma = expected_trk_time->getTH1F()->GetMean();
				if (ResL2_sigma > Exp_sigma)
					Resolution_L2_map_[histo_key] = pow(pow(ResL2_sigma,2)-pow(Exp_sigma,2), 0.5)*1000;
				else
					Resolution_L2_map_[histo_key] = 50;
			}else
				Resolution_L2_map_[histo_key] = 400;
		}else
      Resolution_L2_map_[histo_key] = calib.timePrecision(histo_key);
  }

  // std::cout<<"####### Calib #######"<<std::endl;
  // const auto& calib = iSetup.getData(calibEsToken_);
  // for(int sec=0; sec<2; sec++){
  //   for(int st=0; st<1; st++){
  //     for(int plane=0; plane<4; plane++){
  //       for(int ch=0; ch<12; ch++){
  //         std::cout<<"("<<sec<<", "<<st<<", "<<plane<<", "<<ch<<")"<<" = "<<calib.timePrecision(sec, st, plane, ch)<<std::endl;
  //       }
  //     }
  //   }
  // }
  // DiamondTimingCalibration c(calib);
  // std::cout<<c<<std::endl;

  // const auto& geom = iSetup.getData(geomEsToken_);
  // for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
  //   if (!CTPPSDiamondDetId::check(it->first))
  //     continue;

  //   CTPPSDiamondDetId id((*it).first);
  //   std::cout<<id.arm()<<", "<<id.station()<<", "<<id.plane()<<", "<<id.channel()<<std::endl;
  // }

  iBooker.setCurrentFolder("/");
  for(auto e : Resolution_L2_map_){
    auto* mEl = iBooker.book1D("res_" + std::to_string(e.first.planeKey.sector) + "_" + std::to_string(e.first.planeKey.station) + "_" + std::to_string(e.first.planeKey.plane) + "_" + std::to_string(e.first.channel), "title;x;y", 1200, -60., 60.);
    mEl->Fill(e.second);
  }

  JSON::save(geom, calib, Resolution_L2_map_, "DiamondCalibrationOut.json");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingHarvester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setAllowAnything();
  // desc.add<std::string>("folder", "MY_FOLDER");
  // descriptions.add("DiamondTimingHarvester", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingHarvester);
