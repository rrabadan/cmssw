#ifndef HLTmumuV0Filter_h
#define HLTmumuV0Filter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
namespace edm {
  class ConfigurationDescriptions;
}
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

// ----------------------------------------------------------------------

class HLTmumuV0Filter : public HLTFilter {

 public:
  explicit HLTmumuV0Filter(const edm::ParameterSet&);
  ~HLTmumuV0Filter();
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
  virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct) const override;

 private:

  edm::InputTag                                          muCandTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> muCandToken_;
  edm::InputTag                                          trkCandTag_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> trkCandToken_;
  edm::InputTag                                          MuMuVertexTag_;
  edm::EDGetTokenT<reco::VertexCollection>               MuMuVertexToken_;
  edm::InputTag                                          V0VertexTag_;
  edm::EDGetTokenT<reco::VertexCollection>               V0VertexToken_;

  const double maxDelR_;
  const double minCosinePointingAngle_;
  const double displacement_;

  static bool triggerdByPreviousLevel(const reco::RecoChargedCandidateRef &, const std::vector<reco::RecoChargedCandidateRef> &);

};
#endif
