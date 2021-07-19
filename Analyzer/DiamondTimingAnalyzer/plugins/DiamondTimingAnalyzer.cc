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
//         Created:  Mon, 19 Jul 2021 13:22:33 GMT
//
//

#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

struct Histograms_DiamondTimingAnalyzer {
 dqm::reco::MonitorElement* histo_;
};

class DiamondTimingAnalyzer : public DQMGlobalEDAnalyzer<Histograms_DiamondTimingAnalyzer> {
public:
  explicit DiamondTimingAnalyzer(const edm::ParameterSet&);
  ~DiamondTimingAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&,
                      edm::Run const&,
                      edm::EventSetup const&,
                      Histograms_DiamondTimingAnalyzer&) const override;

  void dqmAnalyze(edm::Event const&, edm::EventSetup const&, Histograms_DiamondTimingAnalyzer const&) const override;

  // ------------ member data ------------
  std::string folder_;
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
DiamondTimingAnalyzer::DiamondTimingAnalyzer(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")) {
  // now do what ever initialization is needed
}

DiamondTimingAnalyzer::~DiamondTimingAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------

void DiamondTimingAnalyzer::dqmAnalyze(edm::Event const& iEvent,
                           edm::EventSetup const& iSetup,
                           Histograms_DiamondTimingAnalyzer const& histos) const {
  histos.histo_->Fill(1.);
}


void DiamondTimingAnalyzer::bookHistograms(DQMStore::IBooker& ibook,
                               edm::Run const& run,
                               edm::EventSetup const& iSetup,
                               Histograms_DiamondTimingAnalyzer& histos) const {
  ibook.setCurrentFolder(folder_);
  histos.histo_ = ibook.book1D("EXAMPLE", "EXAMPLE", 10, 0., 10.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiamondTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setAllowAnything();
  // desc.add<std::string>("folder", "MY_FOLDER");
  // descriptions.add("diamondtiminganalyzer", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DiamondTimingAnalyzer);
