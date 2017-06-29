#include <algorithm>
#include <cmath>
#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "HLTrigger/btau/src/HLTmumuV0Filter.h"
#include "TMath.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;
namespace{
  const double Pi = 3.14159265359;
}

// ----------------------------------------------------------------------
HLTmumuV0Filter::HLTmumuV0Filter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig),
  muCandTag_  (iConfig.getParameter<edm::InputTag>("MuonTag")),
  muCandToken_(consumes<reco::RecoChargedCandidateCollection>(muCandTag_)),
  trkCandTag_  (iConfig.getParameter<edm::InputTag>("TrackTag")),
  trkCandToken_(consumes<reco::RecoChargedCandidateCollection>(trkCandTag_)),
  MuMuVertexTag_  (iConfig.getParameter<edm::InputTag>("MuMuVertexTag")),
  MuMuVertexToken_(consumes<reco::VertexCollection>(MuMuVertexTag_)),
  V0VertexTag_  (iConfig.getParameter<edm::InputTag>("V0VertexTag")),
  V0VertexToken_(consumes<reco::VertexCollection>(V0VertexTag_)),
  maxDelR_(iConfig.getParameter<double>("MaxDeltaR")),
  minCosinePointingAngle_(iConfig.getParameter<double>("MinCosinePointingAngle")),
  displacement_(iConfig.getParameter<double>("Displacement"))
{
}

// ----------------------------------------------------------------------
HLTmumuV0Filter::~HLTmumuV0Filter() = default;

void
HLTmumuV0Filter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<double>("MaxDeltaR",2.5);
  desc.add<double>("Displacement" ,0.0);
  desc.add<double>("MinCosinePointingAngle" ,0.0);
  desc.add<edm::InputTag>("MuonTag",edm::InputTag("hltL3MuonCandidates"));
  desc.add<edm::InputTag>("TrackTag",edm::InputTag("hltMumukAllConeTracks"));
  desc.add<edm::InputTag>("MuMuVertexTag",edm::InputTag("hltDisplacedmumuVtxProducerDoubleMu4Jpsi"));
  desc.add<edm::InputTag>("V0VertexTag",edm::InputTag("hltDisplacedtktkV0Vtx"));
  descriptions.add("HLTmumuV0Filter",desc);
}

// ----------------------------------------------------------------------
bool HLTmumuV0Filter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct) const {

  //get the beamspot position
  //edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  //iEvent.getByToken(beamSpotToken_,recoBeamSpotHandle);
  //const reco::BeamSpot& vertexBeamSpot = *recoBeamSpotHandle;

  // get vertices
  reco::VertexCollection displacedMuMuVertexColl;
  edm::Handle<reco::VertexCollection> displacedMuMuVertexCollHandle;
  bool foundMuMuVertexColl = iEvent.getByToken(MuMuVertexToken_, displacedMuMuVertexCollHandle);
  if(foundMuMuVertexColl) displacedMuMuVertexColl = *displacedMuMuVertexCollHandle;

  reco::VertexCollection displacedV0VertexColl;
  edm::Handle<reco::VertexCollection> displacedV0VertexCollHandle;
  bool foundV0VertexColl = iEvent.getByToken(V0VertexToken_, displacedV0VertexCollHandle);
  if(foundV0VertexColl) displacedV0VertexColl = *displacedV0VertexCollHandle;

  // get muon collection
  Handle<RecoChargedCandidateCollection> mucands;
  iEvent.getByToken(muCandToken_,mucands);

  // get track candidates around displaced muons
  Handle<RecoChargedCandidateCollection> trkcands;
  iEvent.getByToken(trkCandToken_,trkcands);

  // Ref to Candidate object to be recorded in filter object
  RecoChargedCandidateRef refMu1;
  RecoChargedCandidateRef refMu2;
  RecoChargedCandidateRef refTrk1;
  RecoChargedCandidateRef refTrk2;
    
  if (saveTags()) {
    filterproduct.addCollectionTag(muCandTag_);
    filterproduct.addCollectionTag(trkCandTag_);
  }

  bool triggered = false;

  // loop over vertex collection
  reco::VertexCollection::iterator it;
  for(it = displacedMuMuVertexColl.begin(); it!= displacedMuMuVertexColl.end(); it++){
    reco::Vertex displacedMuMuVertex = *it;

    // check if the vertex actually consists of exactly two muon + 1 track, throw exception if not
    if(displacedMuMuVertex.tracksSize() != 2) throw cms::Exception("BadLogic") << "HLTmumuV0Filter: ERROR: the Jpsi vertex must have " 
                                                                           << "exactly two muons by definition. It now has n trakcs = "
                                                                        << displacedMuMuVertex.tracksSize() << std::endl;

    // get the two tracks from the vertex
    auto trackIt =  displacedMuMuVertex.tracks_begin();
    reco::TrackRef vertextkRef1 =  (*trackIt).castTo<reco::TrackRef>() ;
    trackIt++;
    reco::TrackRef vertextkRef2 =  (*trackIt).castTo<reco::TrackRef>();
    //trackIt++;
    //reco::TrackRef vertextkRef3 =  (*trackIt).castTo<reco::TrackRef>();
    
    reco::VertexCollection::iterator it2;
    for(it2 = displacedV0VertexColl.begin(); it2!= displacedV0VertexColl.end(); it2++){
      reco::Vertex displacedV0Vertex = *it2;

      // check if the vertex actually consists of exactly two muon + 1 track, throw exception if not
      if(displacedV0Vertex.tracksSize() != 2) throw cms::Exception("BadLogic") << "HLTmumuV0Filter: ERROR: the V0 vertex must have " 
                                                                           << "exactly two tracks by definition. It now has n trakcs = "
                                                                           << displacedV0Vertex.tracksSize() << std::endl;
    
      
      auto trackIt2 =  displacedV0Vertex.tracks_begin();
      reco::TrackRef vertextkRef3 =  (*trackIt2).castTo<reco::TrackRef>() ;
      trackIt2++;
      reco::TrackRef vertextkRef4 =  (*trackIt2).castTo<reco::TrackRef>();
      
      // first find the two muon tracks in the muon collection
      reco::RecoChargedCandidateCollection::const_iterator mucand1;
      reco::RecoChargedCandidateCollection::const_iterator mucand2;    
      reco::RecoChargedCandidateCollection::const_iterator tkcand1;    
      reco::RecoChargedCandidateCollection::const_iterator tkcand2;
      
      int iFoundRefs = 0;
      bool threeMuons = false;
      for (auto cand=mucands->begin(); cand!=mucands->end(); cand++) {
        reco::TrackRef tkRef = cand->get<reco::TrackRef>();
        if     (tkRef == vertextkRef1 && iFoundRefs==0) {mucand1 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef1 && iFoundRefs==1) {mucand2 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef1 && iFoundRefs==2) {threeMuons = true;}
        if     (tkRef == vertextkRef2 && iFoundRefs==0) {mucand1 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef2 && iFoundRefs==1) {mucand2 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef2 && iFoundRefs==2) {threeMuons = true;}
        if     (tkRef == vertextkRef3 && iFoundRefs==0) {mucand1 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef3 && iFoundRefs==1) {mucand2 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef3 && iFoundRefs==2) {threeMuons = true;}
        if     (tkRef == vertextkRef3 && iFoundRefs==0) {mucand1 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef3 && iFoundRefs==1) {mucand2 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef3 && iFoundRefs==2) {threeMuons = true;}
        if     (tkRef == vertextkRef4 && iFoundRefs==0) {mucand1 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef4 && iFoundRefs==1) {mucand2 = cand; iFoundRefs++;}
        else if(tkRef == vertextkRef4 && iFoundRefs==2) {threeMuons = true;}
      }
      if(threeMuons) throw cms::Exception("BadLogic") << "HLTmumuV0Filterr: ERROR: the vertex must have "
                                                    << " exactly two muons by definition."  << std::endl;

      bool threeTrks = false;
      int iTrkFoundRefs = 0;
      for (auto cand=trkcands->begin(); cand!=trkcands->end(); cand++) {
        reco::TrackRef tkRef = cand->get<reco::TrackRef>();
        if     (tkRef == vertextkRef1 && iTrkFoundRefs==0) {tkcand1 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef1 && iTrkFoundRefs==1) {tkcand2 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef1 && iTrkFoundRefs==2) {threeTrks = true;}
        if     (tkRef == vertextkRef2 && iTrkFoundRefs==0) {tkcand1 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef2 && iTrkFoundRefs==1) {tkcand2 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef2 && iTrkFoundRefs==2) {threeTrks = true;}
        if     (tkRef == vertextkRef3 && iTrkFoundRefs==0) {tkcand1 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef3 && iTrkFoundRefs==1) {tkcand2 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef3 && iTrkFoundRefs==2) {threeTrks = true;}
        if     (tkRef == vertextkRef4 && iTrkFoundRefs==0) {tkcand1 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef4 && iTrkFoundRefs==1) {tkcand2 = cand; iTrkFoundRefs++;}
        else if(tkRef == vertextkRef4 && iTrkFoundRefs==2) {threeTrks = true;}
      }
      if(threeTrks) throw cms::Exception("BadLogic") << "HLTmumuV0Filterr: ERROR: the vertex must have "
                                                 << " exactly two track by definition."  << std::endl;

      
      // calculate two-track transverse momentum
      math::XYZVector pperp(tkcand1->px() + tkcand2->px(),
                          tkcand1->py() + tkcand2->py(), 
                          0.);
            
      // get vertex position and error to calculate the decay length significance
      reco::Vertex::Point refpoint=displacedMuMuVertex.position();
      GlobalPoint refVertex (refpoint.x(), refpoint.y(), refpoint.z());

      reco::Vertex::Point vpoint=displacedV0Vertex.position();
      GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());

      float deta = vpoint.Eta() - refpoint.Eta();
      float dphi = fabs(vpoint.Phi() - refpoint.Phi()) > Pi ? vpoint.Phi() - refpoint.Phi() : 2*Pi - fabs(vpoint.Phi() - refpoint.Phi()); 

      float dR = sqrt(deta*deta +dphi*dphi);

      GlobalPoint displacementFromRef( -1*((refVertex.x() - secondaryVertex.x())), 
                                          -1*((refVertex.y() - secondaryVertex.y())), 
                                           0 );
      float displacement = std::abs(displacementFromRef.perp());

      //calculate the angle between the decay length and the mumu momentum
      Vertex::Point vperp(displacementFromRef.x(),displacementFromRef.y(),0.);
      float cosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());

      //if (pperp.R()  < minPt_                 ) continue;
      //if ( dR           > maxDelR_               ) continue;
      if ( displacement < displacement_          ) continue;
      if ( cosAlpha     < minCosinePointingAngle_) {}//continue;
      
      triggered = true;
          
      refMu1=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (mucands,distance(mucands->begin(), mucand1)));
      filterproduct.addObject(TriggerMuon,refMu1);
      refMu2=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (mucands,distance(mucands->begin(), mucand2)));
      filterproduct.addObject(TriggerMuon,refMu2);
      refTrk1=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (trkcands,distance(trkcands->begin(),tkcand1)));
      filterproduct.addObject(TriggerTrack,refTrk1);
      refTrk2=RecoChargedCandidateRef( Ref<RecoChargedCandidateCollection> (trkcands,distance(trkcands->begin(),tkcand2)));
      filterproduct.addObject(TriggerTrack,refTrk2);
    }//end loop V0 vertices
  }//end loop Jpsi vertices

  LogDebug("HLTDisplacedMumuTrkFilter") << " >>>>> Result of HLTDisplacedMuMuTrkFilter is "<< triggered;
  return triggered;
}



bool HLTmumuV0Filter::triggerdByPreviousLevel(const reco::RecoChargedCandidateRef & candref, const std::vector<reco::RecoChargedCandidateRef>& vcands){
  unsigned int i=0;
  unsigned int i_max=vcands.size();
  for (;i!=i_max;++i){
    if (candref == vcands[i]) return true;
  }
  return false;
}
