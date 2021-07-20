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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDHarvester.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


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
DiamondTimingHarvester::DiamondTimingHarvester(const edm::ParameterSet& iConfig){
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
    std::cout<<"######## EndRun ########"<<std::endl;
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
