#ifndef __OniaV0Producer_h_
#define __OniaV0Producer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing Onia(chi,psi,etc.) and a track (K, pi, etc.)
 */

class OniaV0Producer : public edm::EDProducer {

 public:
  explicit OniaV0Producer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);
  
  virtual void endJob();
  edm::EDGetTokenT<pat::CompositeCandidateCollection> OniaCollection_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> V0Collection_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  std::vector<double> OniaMassCuts_;
  std::vector<double> OniaV0MassCuts_;
  std::vector<double> V0MassCuts_;
  std::vector<double> MassTraks_;
  double hits_;
  double ptcut_;
  double chi2cut_;
  double IPsigXYcut_;
  bool OnlyBest_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool  IsTheSame(const reco::Candidate& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeOniaV0Candidate(const pat::CompositeCandidate& onia, 
                const pat::CompositeCandidate& V0);
  int candidates;
  int nevents;
  int nonia;
  int nreco;
};

#endif // __OniaV0Producer_h_
