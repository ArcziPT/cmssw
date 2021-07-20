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

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"

//
// class declaration
//

struct Histograms_DiamondTiming {
 dqm::reco::MonitorElement* histo_;
};

class DiamondTimingWorker : public DQMEDAnalyzer {
public:
  explicit DiamondTimingWorker(const edm::ParameterSet&);
  ~DiamondTimingWorker() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker &iBooker, edm::Run const &, edm::EventSetup const &iSetup) override;
  void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) override;

  // ------------ member data ------------
  std::string folder_ = "data";

  Histograms_DiamondTiming histos;
  edm::ESGetToken<CTPPSGeometry, VeryForwardRealGeometryRecord> geomEsToken_;
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
DiamondTimingWorker::DiamondTimingWorker(const edm::ParameterSet& iConfig)
  : geomEsToken_(esConsumes<edm::Transition::BeginRun>()) {
  // now do what ever initialization is needed
}

//
// member functions
//

// ------------ method called for each event  ------------

void DiamondTimingWorker::analyze(const edm::Event &iEvent,
                                  const edm::EventSetup &iSetup){
  histos.histo_->Fill(1.);
}


void DiamondTimingWorker::bookHistograms(DQMStore::IBooker& iBooker,
                                           edm::Run const& run,
                                           edm::EventSetup const& iSetup) {
  iBooker.setCurrentFolder(folder_);
  histos.histo_ = iBooker.book1D("EXAMPLE", "EXAMPLE", 10, 0., 10.);

  const auto& geom = iSetup.getData(geomEsToken_);
  for (auto it = geom.beginSensor(); it != geom.endSensor(); ++it) {
    if (!CTPPSDiamondDetId::check(it->first))
      continue;
    
    const CTPPSDiamondDetId detid(it->first);
    
    if (detid.station() != 1)
      continue;
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
