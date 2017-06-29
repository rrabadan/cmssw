#ifndef HLTDisplacedtktkV0Filter_h
#define HLTDisplacedtktkV0Filter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
namespace edm {
  class ConfigurationDescriptions;
}
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

class HLTDisplacedtktkV0Filter : public HLTFilter {

  public:
    explicit HLTDisplacedtktkV0Filter(const edm::ParameterSet&);
    ~HLTDisplacedtktkV0Filter();
    static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);	
    virtual void beginJob() override ;
    virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct) const override;
    virtual void endJob() override ;

  private:
    bool fastAccept_;
    double minLxySignificance_;
    double maxLxySignificance_;
    double maxNormalisedChi2_;
    double minVtxProbability_;
    double minCosinePointingAngle_;

    edm::InputTag                            DisplacedVertexTag_;
    edm::EDGetTokenT<reco::VertexCollection> DisplacedVertexToken_;
    edm::InputTag                            beamSpotTag_;
    edm::EDGetTokenT<reco::BeamSpot>         beamSpotToken_;
    edm::InputTag                                          TrackTag_;
    edm::EDGetTokenT<reco::RecoChargedCandidateCollection> TrackToken_;

};
#endif
