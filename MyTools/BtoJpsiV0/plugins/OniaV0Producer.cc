#include "MyTools/BtoJpsiV0/plugins/OniaV0Producer.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

OniaV0Producer::OniaV0Producer(const edm::ParameterSet& ps):
  OniaCollection_ (consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("Onia"))),
  V0Collection_   (consumes<pat::CompositeCandidateCollection>(ps.getParameter<edm::InputTag>("V0Col"))),
  thebeamspot_    (consumes<reco::BeamSpot>(ps.getParameter<edm::InputTag>("beamSpotTag"))),
  OniaMassCuts_   (ps.getParameter<std::vector<double>>("OniaMassCuts")),
  OniaV0MassCuts_ (ps.getParameter<std::vector<double>>("OniaV0MassCuts")),
  V0MassCuts_     (ps.getParameter<std::vector<double>>("V0MassCuts")),
  MassTraks_      (ps.getParameter<std::vector<double>>("MassTracks")),
  hits_           (ps.getParameter<bool>("hits")),
  ptcut_          (ps.getParameter<bool>("ptcut")),
  chi2cut_        (ps.getParameter<bool>("chi2cut")),
  IPsigXYcut_     (ps.getParameter<bool>("IPsigXYcutcut")),
  OnlyBest_       (ps.getParameter<bool>("OnlyBest"))
{
  produces<pat::CompositeCandidateCollection>("OniaV0Candidates");
  candidates = 0;
  nevents = 0;
  nonia = 0;
  nreco = 0;
}
 
void OniaV0Producer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::unique_ptr<pat::CompositeCandidateCollection> OniaV0CandColl(new pat::CompositeCandidateCollection);

  // the beamspot 
  reco::Vertex theBeamSpotV;
  edm::Handle<reco::BeamSpot> theBeamSpot;
  event.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  edm::Handle<pat::CompositeCandidateCollection> onia;
  event.getByToken(OniaCollection_,onia);

  edm::Handle<pat::CompositeCandidateCollection> theV0Handle;
  event.getByToken(V0Collection_, theV0Handle);

  uint ncombo = 0;
  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];
  float OniaV0MassMax_ = OniaV0MassCuts_[1];
  float OniaV0MassMin_ = OniaV0MassCuts_[0];
  float V0MassMax_ = V0MassCuts_[1];
  float V0MassMin_ = V0MassCuts_[0];

  int ionia = -1;
  // Note: Dimuon cand are sorted by decreasing vertex probability then first chi is associated with "best" dimuon 
  for (pat::CompositeCandidateCollection::const_iterator oniaCand = onia->begin(); oniaCand != onia->end(); ++oniaCand){
    ionia++;
    // If J/psi use reco mass otherwise use mQ
    float oniaM = oniaCand->mass();
    if ( oniaM < OniaMassMax_  && oniaM > OniaMassMin_ ) {
      
      const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon1"));
      const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon2"));

      if ( !(pmu1->isSoftMuon(theBeamSpotV)) || fabs(pmu1->eta()) > 2.4 ) continue;
      if ( !(pmu2->isSoftMuon(theBeamSpotV)) || fabs(pmu2->eta()) > 2.4 ) continue;
      
      for (pat::CompositeCandidateCollection::const_iterator itV0 = theV0Handle->begin(); itV0 != theV0Handle->end(); ++itV0 ){

        float V0mass = oniaCand->mass();
        if ( V0mass > V0MassMax_  && V0mass < V0MassMin_ ) continue;
        
        if ( IsTheSame(*itV0->daughter(0), *pmu1) || IsTheSame(*itV0->daughter(0), *pmu2) ) continue;
        if ( IsTheSame(*itV0->daughter(1), *pmu1) || IsTheSame(*itV0->daughter(1), *pmu2) ) continue;

        reco::TrackRef reftk1 = (dynamic_cast<const reco::RecoCandidate*>(itV0->daughter("V0")->daughter("track0")))->track();
        reco::TrackRef reftk2 = (dynamic_cast<const reco::RecoCandidate*>(itV0->daughter("V0")->daughter("track1")))->track();

        if (!reftk1->quality(reco::TrackBase::highPurity)) continue;
        if (!reftk2->quality(reco::TrackBase::highPurity)) continue;

        if ( reftk1->numberOfValidHits() < hits_ || reftk1->numberOfValidHits() < hits_ ) continue;
        if ( reftk1->pt() < ptcut_ || reftk2->pt() < ptcut_) continue;
        if ( reftk1->normalizedChi2() < chi2cut_ || reftk2->normalizedChi2() < chi2cut_) continue;
        
        float dxy   = reftk1->dxy(bs.position());
        float d0err = reftk1->d0Error();
        if ( dxy/d0err < IPsigXYcut_ ) continue;
        dxy   = reftk2->dxy(bs.position());
        d0err = reftk2->d0Error();
        if ( dxy/d0err < IPsigXYcut_ ) continue;

        pat::CompositeCandidate OniaV0Cand =  makeOniaV0Candidate(*oniaCand,*itV0);
        
        if ( OniaV0Cand.mass() < OniaV0MassMax_ && OniaV0Cand.mass() > OniaV0MassMin_){
          OniaV0CandColl->push_back(OniaV0Cand);
          candidates++;  
          ncombo++;
        }
      }
    }
     if (OnlyBest_) break;
  }
  if ( ncombo != OniaV0CandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != OniaT ("<<OniaV0CandColl->size()<<")"<< std::endl;
  if ( onia->size() > 0 )  nonia++;
  if ( ncombo > 0 ) nreco++; 
  event.put(std::move(OniaV0CandColl),"OniaV0Candidates"); 
  nevents++;
}

void OniaV0Producer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "OniaTrak Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with Onia candidates " << nonia << std::endl;
  std::cout << "Events with OniaTrak candidates " << nreco << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " OniaTrak candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

bool OniaV0Producer::IsTheSame(const reco::Candidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

const pat::CompositeCandidate
OniaV0Producer::makeOniaV0Candidate(const pat::CompositeCandidate& onia, const pat::CompositeCandidate& v0){

  pat::CompositeCandidate Cand;
  Cand.addDaughter(onia,"onia");
  Cand.addDaughter(v0,"V0");
  Cand.setCharge(onia.charge()+ v0.charge());

  reco::Candidate::LorentzVector vcand = onia.p4() + v0.p4();
  Cand.setP4(vcand);

  return Cand;
}

reco::Candidate::LorentzVector OniaV0Producer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(OniaV0Producer);
