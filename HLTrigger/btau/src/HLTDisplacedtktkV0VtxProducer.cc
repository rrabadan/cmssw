#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "HLTrigger/btau/src/HLTDisplacedtktkV0VtxProducer.h"
#include <DataFormats/Math/interface/deltaR.h>

using namespace edm;
using namespace reco;
using namespace std; 
using namespace trigger;

namespace {
  const double piMass = 0.13957018;
  const double piMassSquared = piMass*piMass;
  const double protonMass = 0.938272046;
  const double protonMassSquared = protonMass*protonMass;
  const double kShortMass = 0.497614;
  const double lambdaMass = 1.115683;
}
//
// constructors and destructor
//
HLTDisplacedtktkV0VtxProducer::HLTDisplacedtktkV0VtxProducer(const edm::ParameterSet& iConfig):	
	muCandTag_ (iConfig.getParameter<edm::InputTag>("MuCand")),
  muCandToken_ (consumes<reco::RecoChargedCandidateCollection>(muCandTag_)),
  trkCandTag_ (iConfig.getParameter<edm::InputTag>("TrackCand")),
	trkCandToken_(consumes<reco::RecoChargedCandidateCollection>(trkCandTag_)),
  beamSpotTag_ (iConfig.getParameter<edm::InputTag> ("BeamSpotTag")),
  beamSpotToken_(consumes<reco::BeamSpot>(beamSpotTag_)),
  previousCandTag_(iConfig.getParameter<edm::InputTag>("PreviousCandTag")),
	previousCandToken_(consumes<trigger::TriggerFilterObjectWithRefs>(previousCandTag_)),
	maxEta_ (iConfig.getParameter<double>("MaxEta")),
	minPt_ (iConfig.getParameter<double>("MinPt")),
  tkChi2Cut_ (iConfig.getParameter<double>("TkChi2Cut")),
	tkNHitsCut_ (iConfig.getParameter<int>("TkNHitsCut")), 
  tkIPsigXYcut_(iConfig.getParameter<double>("tkIPsigXYcut")),
  //tkIPsigZcut_(iConfig.getParameter<double>("tkIPsigZcut")),
  minPtPair_ (iConfig.getParameter<double>("MinPtPair")),
	minInvMass_ (iConfig.getParameter<double>("MinInvMass")),
	maxInvMass_ (iConfig.getParameter<double>("MaxInvMass")),
  kShortMassCut_ (iConfig.getParameter<double>("KShortMassCut")),
  lambdaMassCut_ (iConfig.getParameter<double>("LambdaMassCut")),
  tkDCACut_ (iConfig.getParameter<double>("TkDCACut")),
  overlapDR_(iConfig.getParameter<double>("OverlapDR"))
  //triggerTypeDaughters_(iConfig.getParameter<int>("triggerTypeDaughters"))
{
	produces<VertexCollection>();
}


HLTDisplacedtktkV0VtxProducer::~HLTDisplacedtktkV0VtxProducer() = default;

void
HLTDisplacedtktkV0VtxProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("MuCand",edm::InputTag("hltL3MuonCandidates"));
  desc.add<edm::InputTag>("TrackCand",edm::InputTag("hltJpsiTkAllConeTracksIter"));
  desc.add<edm::InputTag>("BeamSpotTag",edm::InputTag("hltOnlineBeamSpot"));
  desc.add<edm::InputTag>("PreviousCandTag",edm::InputTag("hltDisplacedmumuFilterDoubleMu4Jpsi"));
  desc.add<double>("MaxEta",2.5);
  desc.add<double>("MinPt",0.35);
  desc.add<double>("TkChi2Cut",10.0);
  desc.add<int>("TkNHitsCut",7);
  desc.add<double>("tkIPsigXYcut",2);
  //desc.add<double>("tkIPsigZcut",-1);
  desc.add<double>("MinPtPair",0.0);
  desc.add<double>("MinInvMass",0.0);
  desc.add<double>("MaxInvMass",1.5);
  desc.add<int>("ChargeOpt",-1);
  desc.add<double>("KShortMassCut",1.30);
  desc.add<double>("LambdaMassCut",0.30);
  desc.add<double>("TkDCACut",2.0);
  desc.add<double>("OverlapDR",1.44e-4);
  //desc.add<int>("triggerTypeDaughters",0);
  descriptions.add("hltDisplacedtktkV0VtxProducer", desc);
}

// ------------ method called once each job just before starting event loop  ------------
void HLTDisplacedtktkV0VtxProducer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void HLTDisplacedtktkV0VtxProducer::endJob() 
{
 	
}

// ------------ method called on each new Event  ------------
void HLTDisplacedtktkV0VtxProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	/*double const firstTrackMass = massParticle1_;
	double const firstTrackMass2 = firstTrackMass*firstTrackMass;
	double const secondTrackMass = massParticle2_;
	double const secondTrackMass2 = secondTrackMass*secondTrackMass;
  */
  //get hold of muon tracks
  Handle<RecoChargedCandidateCollection> mucands;
  iEvent.getByToken(muCandToken_,mucands);
  //std::cout << "+++++ HLTDisplacedtktkV0VtxProducer " << std::endl;
	//get the transient track builder:
	edm::ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::ESHandle<MagneticField> theMagneticFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(theMagneticFieldHandle);
  const MagneticField* theMagneticField = theMagneticFieldHandle.product();

	// get track candidate around displaced muon
	Handle<RecoChargedCandidateCollection> trackcands;
	iEvent.getByToken(trkCandToken_,trackcands);

  reco::BeamSpot vertexBeamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_,recoBeamSpotHandle);
  vertexBeamSpot = *recoBeamSpotHandle;

  std::unique_ptr<VertexCollection> vertexCollection(new VertexCollection());

	// look at all trackcands,  check cuts and make vertices
	double e1,e2;
	Particle::LorentzVector p,p1,p2;

  //RecoChargedCandidateCollection::const_iterator mucand1;
  //RecoChargedCandidateCollection::const_iterator mucand2;
	RecoChargedCandidateCollection::const_iterator cand1;
	RecoChargedCandidateCollection::const_iterator cand2;
	
	// get the objects passing the previous filter
	Handle<TriggerFilterObjectWithRefs> previousCands;
	iEvent.getByToken(previousCandToken_,previousCands);

	vector<RecoChargedCandidateRef> vPrevCands;
	previousCands->getObjects(TriggerMuon,vPrevCands);
	
	for (cand1=trackcands->begin(); cand1!=trackcands->end(); cand1++) {
      TrackRef tk1 = cand1->get<TrackRef>();
	    if ( tk1->numberOfValidHits() < tkNHitsCut_ || tk1->normalizedChi2() > tkChi2Cut_ ) continue;
      double ipsigXY = std::abs(tk1->dxy(vertexBeamSpot)/tk1->dxyError());
      //double ipsigXY = std::abs(tk1->dz(vertexBeamSpot)/tk1->dzError());
      LogDebug("HLTDisplacedtktkV0VtxProducer") << " 1st track in loop: q*pt= " << cand1->charge()*cand1->pt() << ", eta= " << cand1->eta() << ", hits= " << tk1->numberOfValidHits();
      
      // cuts
      if ( ipsigXY < tkIPsigXYcut_ ) continue;
      if (abs(cand1->eta())>maxEta_) continue;
	    if (cand1->pt() < minPt_) continue;
	    
 	
	    cand2=trackcands->begin();
	    //if(massParticle1_==massParticle2_){cand2 = cand1+1;}
      cand2 = cand1+1;

	    for (; cand2!=trackcands->end(); cand2++) {
        TrackRef tk2 = cand2->get<TrackRef>();
		    if(tk1==tk2) continue;
        if ( tk2->numberOfValidHits() < tkNHitsCut_ || tk2->normalizedChi2() > tkChi2Cut_ ) continue;
        ipsigXY = std::abs(tk2->dxy(vertexBeamSpot)/tk2->dxyError());
        //double ipsigXY = std::abs(tk2->dz(vertexBeamSpot)/tk2->dzError());
        LogDebug("HLTDisplacedtktkV0VtxProducer") << " 2nd track in loop: q*pt= " << cand2->charge()*cand2->pt() << ", eta= " << cand2->eta() << ", hits= " << tk2->numberOfValidHits() << ", d0= " << tk2->d0();
        
        // cuts
        if ( ipsigXY < tkIPsigXYcut_ ) continue;
        if (abs(cand2->eta())>maxEta_) continue;
			  if (cand2->pt() < minPt_) continue;
			 
        // opposite sign or same sign
			  //if (chargeOpt_<0) {
			  if (cand1->charge()*cand2->charge()>0) continue;
			  //} else if (chargeOpt_>0) {
			  //  if (cand1->charge()*cand2->charge()<0) continue;
			  //}

        bool Overlap1 = false;
        bool Overlap2 = false;

        for ( RecoChargedCandidateCollection::const_iterator mucand = mucands->begin(); mucand != mucands->end(); ++mucand ){
          TrackRef mutrk = mucand->get<TrackRef>();
          //first check if this muon passed the previous filter
          if( ! checkPreviousCand( mutrk, vPrevCands) ) continue;
        
          // eta and pt cut
          if (fabs(mutrk->eta()) > maxEta_)    continue;
          if (mutrk->pt()        < minPt_ )    continue;
          Overlap1 = overlap( mutrk, tk1);
          Overlap2 = overlap( mutrk, tk2);
          if( Overlap1 || Overlap2 ) break;
        }

        if ( Overlap2 ) continue;
        if ( Overlap1 ) break;
        // measure distance between tracks at their closest approach
        reco::TransientTrack* negTransTkPtr = nullptr;
        reco::TransientTrack* posTransTkPtr = nullptr;
        if ( tk1->charge() < 0 && tk2->charge() > 0){
          negTransTkPtr = new TransientTrack(*tk1, theMagneticField);
          posTransTkPtr = new TransientTrack(*tk2, theMagneticField);
        }
        else if ( tk1->charge() > 0 && tk2->charge() < 0){
          posTransTkPtr = new TransientTrack(*tk1, theMagneticField);
          negTransTkPtr = new TransientTrack(*tk2, theMagneticField);
        }
        else 
          continue;
      
        if ( !negTransTkPtr->impactPointTSCP().isValid() || !posTransTkPtr->impactPointTSCP().isValid() ) continue;
        FreeTrajectoryState const & state1 = negTransTkPtr->impactPointTSCP().theState();
        FreeTrajectoryState const & state2 = posTransTkPtr->impactPointTSCP().theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(state1, state2);
        
        if (!cApp.status()) continue;
        float dca = std::abs(cApp.distance());
        if (dca > tkDCACut_) continue;
        
        
			  // Combined ditrack system
			  e1 = sqrt(cand1->momentum().Mag2()+piMassSquared);
			  e2 = sqrt(cand2->momentum().Mag2()+piMassSquared);
			  p1 = Particle::LorentzVector(cand1->px(),cand1->py(),cand1->pz(),e1);
			  p2 = Particle::LorentzVector(cand2->px(),cand2->py(),cand2->pz(),e2);
			  p = p1+p2;
			  if (p.pt()<minPtPair_) continue;
			 
			  double invmass = abs(p.mass());
			  LogDebug("HLTDisplacedtktkV0VtxProducer") << " ... 1-2 invmass= " << invmass;
			 
			  if (invmass<minInvMass_) continue;
			  if (invmass>maxInvMass_) continue;
			  
			  // do the vertex fit
			  vector<TransientTrack> t_tks;
			  TransientTrack ttkp1 = (*theB).build(&tk1);
			  TransientTrack ttkp2 = (*theB).build(&tk2);
			  t_tks.push_back(ttkp1);
			  t_tks.push_back(ttkp2);
			 
			  KalmanVertexFitter kvf;
			  TransientVertex tv = kvf.vertex(t_tks);

			  if (!tv.isValid()) continue;

        Vertex vertex = tv;
        GlobalPoint vtxPos(vertex.x(), vertex.y(), vertex.z());

        std::auto_ptr<TrajectoryStateClosestToPoint> trajPlus;
        std::auto_ptr<TrajectoryStateClosestToPoint> trajMins;
        std::vector<reco::TransientTrack> theRefTracks;
        if (tv.hasRefittedTracks()) {
          theRefTracks = tv.refittedTracks();
        } 
        if ( theRefTracks.size() > 1) {
          reco::TransientTrack* thePositiveRefTrack = 0;
          reco::TransientTrack* theNegativeRefTrack = 0;
          for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) {
            if (iTrack->track().charge() > 0.) {
              thePositiveRefTrack = &*iTrack;
            } else if (iTrack->track().charge() < 0.) {
              theNegativeRefTrack = &*iTrack;
            }
          }
          if (thePositiveRefTrack == 0 || theNegativeRefTrack == 0) continue;
          trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
          trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
        } 
        else {
          trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
          trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
        }
    
        if (trajPlus.get() == 0 || trajMins.get() == 0 || !trajPlus->isValid() || !trajMins->isValid()) continue;
    
        GlobalVector positiveP(trajPlus->momentum());
        GlobalVector negativeP(trajMins->momentum());
        GlobalVector totalP(positiveP + negativeP);
    
        // calculate total energy of V0 3 ways: assume it's a kShort, a Lambda, or a LambdaBar.
        double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
        double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
        double protonE = sqrt(positiveP.mag2() + protonMassSquared);
        double antiProtonE = sqrt(negativeP.mag2() + protonMassSquared);
        double kShortETot = piPlusE + piMinusE;
        double lambdaEtot = protonE + piMinusE;
        double lambdaBarEtot = antiProtonE + piPlusE;

        // Create momentum 4-vectors for the 3 candidate types
        const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);
        const reco::Particle::LorentzVector lambdaP4(totalP.x(), totalP.y(), totalP.z(), lambdaEtot);
        const reco::Particle::LorentzVector lambdaBarP4(totalP.x(), totalP.y(), totalP.z(), lambdaBarEtot);

        if ( kShortP4.M() > kShortMass + kShortMassCut_ && kShortP4.M() < kShortMass - kShortMassCut_) continue;
        if ( lambdaP4.M() > lambdaMass + lambdaMassCut_ && lambdaP4.M() < lambdaMass - lambdaMassCut_) continue;
        if ( lambdaBarP4.M() > lambdaMass + lambdaMassCut_ && lambdaBarP4.M() < lambdaMass - lambdaMassCut_) continue;
        
			  // put vertex in the event
			  vertexCollection->push_back(vertex);
        //delete negTransTkPtr;
        //delete posTransTkPtr;
        //negTransTkPtr = posTransTkPtr = nullptr;
      }
    }
   	iEvent.put(std::move(vertexCollection));
}

bool HLTDisplacedtktkV0VtxProducer::overlap(const TrackRef& trackref1, const TrackRef& trackref2){
  if (deltaR(trackref1->eta(), trackref1->phi(),trackref2->eta(), trackref2->phi()) < overlapDR_) return 1;
  return 0;
}

bool HLTDisplacedtktkV0VtxProducer::checkPreviousCand(const TrackRef& trackref, vector<RecoChargedCandidateRef> & refVect){
  bool ok=false;
  for (auto & i : refVect) {
    if ( i->get<TrackRef>() == trackref ) {
      ok=true;
      break;
    }
  }
  return ok;
}
