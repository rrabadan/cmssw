// -*- C++ -*-
//
// Package:    MyTools/HLTbits
// Class:      HLTbits
// 
/**\class HLTbits HLTbits.cc MyTools/HLTbits/plugins/HLTbits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raul Iraq Rabadan Trejo
//         Created:  Thu, 15 Jun 2017 21:49:21 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HLTbits : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HLTbits(const edm::ParameterSet&);
      ~HLTbits();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
      //std::vector<std::string> triggNames_;

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
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
HLTbits::HLTbits(const edm::ParameterSet& iConfig):
  triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerLabel")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


HLTbits::~HLTbits()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HLTbits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 UInt_t trigger=0;
 edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
 iEvent.getByToken(triggerToken_,hltTriggerResultHandle);

 if(hltTriggerResultHandle.isValid())
 {
   const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hltTriggerResultHandle);
   if ( triggerNames.size()>0 ) for (unsigned int l=0; l < triggerNames.size() ;l++) std::cout << triggerNames.triggerName(l) << std::endl;
  std::cout << trigger << std::endl;
 }
}

// ------------ method called once each job just before starting event loop  ------------
void 
HLTbits::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTbits::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTbits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTbits);
