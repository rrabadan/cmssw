#ifndef HLTDisplacedtktkV0VtxProducer_h
#define HLTDisplacedtktkV0VtxProducer_h

/** \class HLTDisplacedtktkVtxProducer_h
 *
 *  
 *  produces kalman vertices from di-track
 *  takes track candidates as input
 *  configurable cuts on pt, eta, pair pt, inv. mass
 *  the two tracks have to be both tracks or both muons
 *
 *  \author Alexander.Schmidt@cern.ch
 *  \date   3. Feb. 2011
 *  \adapted for D Gian.Michele.Innocenti@cern.ch
 *  \date   5. Aug. 2015
 *
 */



#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"
#include <vector>

namespace edm {
  class ConfigurationDescriptions;
}

class HLTDisplacedtktkV0VtxProducer : public edm::stream::EDProducer<> {
 public:
  explicit HLTDisplacedtktkV0VtxProducer(const edm::ParameterSet&);
  ~HLTDisplacedtktkV0VtxProducer();
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);  
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:  

  bool overlap(const reco::TrackRef& trackref1, const reco::TrackRef& trackref2);
  bool checkPreviousCand(const reco::TrackRef& trackref, std::vector<reco::RecoChargedCandidateRef>& ref2);

  const edm::InputTag                                          muCandTag_;
  const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> muCandToken_;
  const edm::InputTag                                          trkCandTag_;
  const edm::EDGetTokenT<reco::RecoChargedCandidateCollection> trkCandToken_;
  const edm::InputTag                                          beamSpotTag_;
  const edm::EDGetTokenT<reco::BeamSpot>                       beamSpotToken_;
  const edm::InputTag                                          previousCandTag_;
  const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> previousCandToken_;
  const double maxEta_;
  const double minPt_;
  const double tkChi2Cut_;
  const int    tkNHitsCut_;
  const double tkIPsigXYcut_;
  //const double tkIPsigZcut_;
  const double minPtPair_;
  const double minInvMass_;
  const double maxInvMass_;
  const double kShortMassCut_;
  const double lambdaMassCut_;
  const double tkDCACut_;
  const double overlapDR_;
  //const int triggerTypeDaughters_;  

};

#endif
