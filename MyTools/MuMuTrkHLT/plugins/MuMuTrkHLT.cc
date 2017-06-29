// -*- C++ -*-
//
// Package:    MyTools/MuMuTrkHLT
// Class:      MuMuTrkHLT
// 
/**\class MuMuTrkHLT MuMuTrkHLT.cc MyTools/MuMuTrkHLT/plugins/MuMuTrkHLT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raul Iraq Rabadan Trejo
//         Created:  Fri, 16 Jun 2017 09:46:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"

#include "MyTools/MuMuTrkHLT/src/ntupleTree.h"
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuTrkHLT : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTrkHLT(const edm::ParameterSet&);
      ~MuMuTrkHLT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginEvent();
      
      void fillHlt        (const edm::Handle<edm::TriggerResults> &, 
                           const edm::Handle<trigger::TriggerEvent> &,
                           const edm::TriggerNames &,
                           const edm::Event &,
                           bool       isTag);
      void fillGen        (const edm::Handle<reco::GenParticleCollection> & genParticles,
                           const edm::Event                               & event );
      void fillHltMuons   (const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                           const edm::Event                                        & event );
      void fillHltDiMuons (const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands ,
                           const edm::Event                                        & event   ,
                           const edm::EventSetup                                   & eventSetup);
      void fillHltTracks  (const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands ,
                           const edm::Event                                        & event    ); 
      void fillHltV0      (const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands ,
                           const edm::Event                                        & event   ,
                           const edm::EventSetup                                   & eventSetup); 
      void fillHltMuVtx   (const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                           const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                           const edm::Event                                        & event);
      void fillHltV0Vtx   (const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                           const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands,
                           const edm::Event                                        & event);
      
      
      // ----------member data ---------------------------
      // Services
      edm::Service<TFileService> outfile_;

      edm::InputTag triggerResultTag_;
      edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
      edm::InputTag triggerSummTag_;
      edm::EDGetTokenT<trigger::TriggerEvent> triggerSummToken_;

      edm::InputTag puTag_;
      edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;
      edm::InputTag offlinePVTag_;
      edm::EDGetTokenT<reco::VertexCollection> offlinePVToken_;
      edm::InputTag beamspotTag_;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
      edm::InputTag lumiScalerTag_;
      edm::EDGetTokenT<LumiScalersCollection> lumiScalerToken_; 
      
      edm::InputTag genTag_;
      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

      edm::InputTag l3candTag_;
      edm::EDGetTokenT<reco::RecoChargedCandidateCollection> l3candToken_;
      edm::InputTag tkcandTag_;
      edm::EDGetTokenT<reco::RecoChargedCandidateCollection> tkcandToken_;
      
      edm::InputTag mumuVtxTag_;
      edm::EDGetTokenT<reco::VertexCollection> mumuVtxToken_;
      edm::InputTag  v0VtxTag_;
      edm::EDGetTokenT<reco::VertexCollection> v0VtxToken_;
      //edm::EDGetTokenT<std::vector<reco::Track>> TrakCollection_;
      //
      
      edm::InputTag bcandTag_;
      edm::EDGetTokenT<pat::CompositeCandidateCollection> bcandToken_;
      
      ntupleEvent event_;
      std::map<std::string,TTree*> outTree_;
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
MuMuTrkHLT::MuMuTrkHLT(const edm::ParameterSet& iConfig):
triggerResultTag_       (iConfig.getUntrackedParameter<edm::InputTag>("triggerResult")), 
triggerResultToken_     (consumes<edm::TriggerResults>(triggerResultTag_)),
triggerSummTag_         (iConfig.getUntrackedParameter<edm::InputTag>("triggerSummary")), 
triggerSummToken_       (consumes<trigger::TriggerEvent>(triggerSummTag_)),
//puTag_                  (iConfig.getUntrackedParameter<edm::InputTag>("puInfoTag")),
//puToken_                (consumes<std::vector< PileupSummaryInfo>>(puTag_)), 
/*offlinePVTag_           (iConfig.getParameter<edm::InputTag>("offlineVtx")), 
offlinePVToken_         (consumes<reco::VertexCollection>(offlinePVTag_)), 
beamspotTag_            (iConfig.getParameter<edm::InputTag>("beamspot")), 
beamspotToken_          (consumes<reco::BeamSpot>(beamspotTag_)), 
//lumiScalerTag_          (iConfig.getUntrackedParameter<edm::InputTag>("lumiScalerTag")),
//lumiScalerToken_        (consumes<LumiScalersCollection>(lumiScalerTag_)), 
*/genTag_                 (iConfig.getUntrackedParameter<edm::InputTag>("genParticlesTag")),
genToken_               (consumes<reco::GenParticleCollection>(genTag_)), 
l3candTag_              (iConfig.getParameter<edm::InputTag>("L3CandidatesTag")),
l3candToken_            (consumes<reco::RecoChargedCandidateCollection>(l3candTag_)),
tkcandTag_              (iConfig.getParameter<edm::InputTag>("TkCandidatesTag")),
tkcandToken_            (consumes<reco::RecoChargedCandidateCollection>(tkcandTag_)),
mumuVtxTag_             (iConfig.getUntrackedParameter<edm::InputTag>("MuMuVtxTag")),
mumuVtxToken_           (consumes<reco::VertexCollection>(mumuVtxTag_)),
v0VtxTag_               (iConfig.getUntrackedParameter<edm::InputTag>("V0VtxTag")),
v0VtxToken_             (consumes<reco::VertexCollection>(v0VtxTag_))
/*bcandTag_               (iConfig.getUntrackedParameter<edm::InputTag>("bCandTag")),
bcandToken_             (consumes<pat::CompositeCandidateCollection>(bcandTag_))
*/{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


MuMuTrkHLT::~MuMuTrkHLT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTrkHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   beginEvent();
   // Fill general info
   event_.runNumber             = iEvent.id().run();
   event_.luminosityBlockNumber = iEvent.id().luminosityBlock();
   event_.eventNumber           = iEvent.id().event();

   // Fill trigger information for probe muon
  edm::Handle<edm::TriggerResults>   triggerResults;
  edm::Handle<trigger::TriggerEvent> triggerEvent;

  if (iEvent.getByToken(triggerResultToken_, triggerResults) &&
      iEvent.getByToken(triggerSummToken_  , triggerEvent)) {
      
    edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResults);
    fillHlt(triggerResults, triggerEvent, triggerNames_, iEvent, false);
  }
  else 
    edm::LogError("") << "Trigger collection for probe muon not found !!!";

  /*if (!iEvent.isRealData()) {
    edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
    if ( iEvent.getByToken(puToken_,puInfo)){
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) 
      {
        if(PVI->getBunchCrossing()==0){
          event_.trueNI   = PVI->getTrueNumInteractions();
          continue;
        }
      }
    } 
    else  
      edm::LogError("") << "PU collection not found !!!";
  }
  */
  // Fill MC GEN info
  /*if (!iEvent.isRealData()) {
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genToken_, genParticles);
      fillGen(genParticles, iEvent);
  }  
  

  // Handle the online muon collection and fill online muons
  edm::Handle<reco::RecoChargedCandidateCollection> l3cands;
  if (iEvent.getByToken(l3candToken_, l3cands)){
    fillHltMuons  (l3cands, iEvent);
    fillHltDiMuons(l3cands, iEvent, iSetup);
  }  
  else
    edm::LogWarning("") << "Online muon collection not found !!!";


  // Handle the online track collection and fill online tracks
  edm::Handle<reco::RecoChargedCandidateCollection> tkcands;
  if (iEvent.getByToken(tkcandToken_, tkcands)){
    fillHltTracks(tkcands, iEvent);
    fillHltV0(tkcands, iEvent, iSetup);
  }
  else
    edm::LogWarning("") << "Online track collection not found !!!";


  //Handle the online mumu vtx collection and fill container
  edm::Handle<reco::VertexCollection> hlt_muvtx;
  if (iEvent.getByToken(mumuVtxToken_, hlt_muvtx) && 
      iEvent.getByToken( l3candToken_, l3cands  ) )
    fillHltMuVtx(hlt_muvtx, l3cands, iEvent);
  else
    edm::LogWarning("") << "Online dimuon vertex collection not found !!!";

  // Handle the online v0 vtx collection and fill container
  edm::Handle<reco::VertexCollection> hlt_v0vtx;
  if (iEvent.getByToken(v0VtxToken_, hlt_v0vtx) && 
      iEvent.getByToken( tkcandToken_, tkcands  ) )
    fillHltV0Vtx(hlt_v0vtx, tkcands, iEvent);
  else
    edm::LogWarning("") << "Online dimuon vertex collection not found !!!";
    */

  /*
  //edm::Handle<std::vector<reco::Track> > trak;
  //edm::Handle<std::vector<pat::GenericParticle> > trak;
  //iEvent.getByToken(TrakCollection_,trak);
  //for (std::vector<reco::Track>::const_iterator trakCand = trak->begin(), trakend=trak->end(); trakCand!= trakend; ++trakCand){
  //for (std::vector<pat::GenericParticle>::const_iterator trakCand = trak->begin(), trakend=trak->end(); trakCand!= trakend; ++trakCand){
    //double trakpt = trakCand->pt();
    //if ( trakpt < 0.7 ) continue;
  //}
   */outTree_["ntupleTree"] -> Fill();
}

void MuMuTrkHLT::fillHlt( const edm::Handle<edm::TriggerResults>   & triggerResults, 
                          const edm::Handle<trigger::TriggerEvent> & triggerEvent  ,
                          const edm::TriggerNames                  & triggerNames  ,
                          const edm::Event                         & event         ,
                          bool                                       isTag         )
{    
   
  for (unsigned int itrig=0; itrig < triggerNames.size(); ++itrig) 
  {
    LogDebug ("triggers") << triggerNames.triggerName(itrig) ;
    if (triggerResults->accept(itrig)) 
    {
      std::string pathName = triggerNames.triggerName(itrig);
      if (isTag) event_.hltTag.triggers.push_back(pathName);
      else       event_.hlt   .triggers.push_back(pathName);
    }
  }
     
     
  const trigger::size_type nFilters(triggerEvent->sizeFilters());
  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
  {
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();

    if ( ( filterTag.find ("Double"        ) !=std::string::npos  ||
           filterTag.find ("Jpsi"          ) !=std::string::npos  ||
           filterTag.find ("LowMass"       ) !=std::string::npos  ||
           filterTag.find ("PsiPrime"      ) !=std::string::npos  ||
           filterTag.find ("Displaced"     ) !=std::string::npos 
//            filterTag.find ("DoubleMu"  ) !=std::string::npos ||
//            filterTag.find ("DiMuonGlb" ) !=std::string::npos
           ) 
            //&& 
           //filterTag.find ("Tau"       ) ==std::string::npos   &&
           //filterTag.find ("EG"        ) ==std::string::npos   &&
           //filterTag.find ("MultiFit"  ) ==std::string::npos  
       )   
    {
      trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
      const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
      
      for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
      {  
        trigger::size_type objKey = objectKeys.at(iKey);
        const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
        
        HLTObjCand hltObj;
        
        hltObj.filterTag = filterTag;
  
        hltObj.pt  = triggerObj.pt();
        hltObj.eta = triggerObj.eta();
        hltObj.phi = triggerObj.phi();
        
        if (isTag)       event_.hltTag.objects.push_back(hltObj);
        else             event_.hlt   .objects.push_back(hltObj);
      }    
    }       
  }

}

void MuMuTrkHLT::fillGen(const edm::Handle<reco::GenParticleCollection> & genParticles ,
                         const edm::Event                               & event )
{

  int B0Id     =    511;
  int JpsiId   =    443;
  int kshort   =    310;
  //int kplus    =    321;
  int piplus   =    211;
  int muId     =    13; 

  for ( size_t i=0; i< genParticles->size(); ++i) 
  {   
    const reco::GenParticle &p = (*genParticles)[i];
    // only save interesting p
    int id = fabs( p.pdgId() );
    if (id != B0Id) continue; 
    
    GenParticleCand theGen;
    theGen.pdgId  = p.pdgId();
    theGen.pt     = p.pt() ;
    theGen.eta    = p.eta();
    theGen.phi    = p.phi();
    theGen.energy = p.energy();
    theGen.status = p.status();
    
    unsigned int n_moms = p.numberOfMothers();
    if (n_moms == 0 ){
      theGen.pdgMother.push_back(0);
      theGen.pdgRealMother.push_back(0);
    }   
    else {
      for (unsigned int im=0; im < n_moms; ++im){
        theGen.pdgMother.push_back(p.motherRef(im)->pdgId());
      }   
    }   
    
    bool boolJpsi   = false;
    bool boolKshort = false;
    bool boolMuM    = false;
    bool boolMuP    = false;
    for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
    {   
      const reco::Candidate *des = p.daughter(ides);
      if( des->pdgId() == JpsiId ) 
      {   
        boolJpsi = true;
        if (des->numberOfDaughters()!=2)      continue;  
        for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
        {   
          const reco::Candidate *mumu = des->daughter(imumu);
          if ( mumu->pdgId() == -muId ) 
          {   
            theGen.mumPt  = mumu->pt() ;
            theGen.mumEta = mumu->eta();
            theGen.mumPhi = mumu->phi();
          }   
          else if ( mumu->pdgId() == muId )
          {   
            theGen.mupPt  = mumu->pt() ;
            theGen.mupEta = mumu->eta();
            theGen.mupPhi = mumu->phi();
          }   
          else continue;    
        }   
      }   
      else if( des->pdgId() == -muId ) 
      {   
        boolMuM = true;
        theGen.mumPt  = des->pt() ;
        theGen.mumEta = des->eta();
        theGen.mumPhi = des->phi();
      }   
      else if( des->pdgId() == muId ) 
      {   
        boolMuP = true;
        theGen.mupPt  = des->pt() ;
        theGen.mupEta = des->eta();
        theGen.mupPhi = des->phi();
      }   
      if( fabs(des->pdgId()) == kshort)  
      {   
        boolKshort = true;
        for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
        {   
          const reco::Candidate *mumu = des->daughter(imumu);
          if( mumu->pdgId() == piplus )   {   
            theGen.pipPt  = mumu->pt() ;
            theGen.pipEta = mumu->eta();
            theGen.pipPhi = mumu->phi();
          }   
          else if ( mumu->pdgId() == -piplus )   {   
            theGen.pimPt  = mumu->pt() ;
            theGen.pimEta = mumu->eta();
            theGen.pimPhi = mumu->phi();
          }   
          else continue;    
        }   
      }
    } // end for des
    if (!( (boolJpsi || (boolMuP && boolMuM) ) && boolKshort)) continue;
    event_.genParticles.push_back(theGen);
  }  // end for genParticles
}


void 
MuMuTrkHLT::fillHltMuons (const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                          const edm::Event                                        & event )
{
  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    HLTMuCand theL3Mu;

    reco::RecoChargedCandidateRef candref(l3cands, il3);
    theL3Mu.pt      = candref -> pt();
    theL3Mu.eta     = candref -> eta();
    theL3Mu.phi     = candref -> phi();
    theL3Mu.charge  = candref -> charge();

    event_.hlt_mu.push_back(theL3Mu);
  }
}

void
MuMuTrkHLT::fillHltDiMuons(const edm::Handle<reco::RecoChargedCandidateCollection>& l3cands, const edm::Event& event,
               const edm::EventSetup& eventSetup)
{
  if (l3cands->size() < 2) return;
  HLTDimuonCand theDimuon;
  
  for( unsigned int il3 = 0; il3 < l3cands->size(); ++il3) 
  {
    //     HLTMuCand oneL3Mu;
    reco::RecoChargedCandidateRef oneref(l3cands, il3);
    
    theDimuon.mu1pt  = oneref -> pt();
    theDimuon.mu1eta = oneref -> eta();
    theDimuon.mu1phi = oneref -> phi();
  
    theDimuon.DCA    = -20  ;

    for( unsigned int jl3 = il3+1; jl3 < l3cands->size(); ++jl3) 
    {
      reco::RecoChargedCandidateRef tworef(l3cands, jl3);
      //if ( tworef -> charge() * oneref -> charge() > 0) continue;
       
      theDimuon.mu2pt  = tworef -> pt();
      theDimuon.mu2eta = tworef -> eta();
      theDimuon.mu2phi = tworef -> phi();
    
      theDimuon.charge =  tworef -> charge() * oneref -> charge();
      event_.hlt_dimu.push_back(theDimuon);
    } 
  }
}

void
MuMuTrkHLT::fillHltTracks(const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands,
                          const edm::Event                                        & event   ) 
{
  for( unsigned int itk = 0; itk < trkcands->size(); ++itk) 
  {
    HLTTkCand theTk;

    reco::RecoChargedCandidateRef candref(trkcands, itk);
    theTk.pt      = candref -> pt();
    theTk.eta     = candref -> eta();
    theTk.phi     = candref -> phi();
    theTk.charge  = candref -> charge();
    
    event_.hlt_tk.push_back(theTk);
  }
}

void
MuMuTrkHLT::fillHltV0(const edm::Handle<reco::RecoChargedCandidateCollection>& trkcands, const edm::Event& event,
               const edm::EventSetup& eventSetup)
{
  if (trkcands->size() < 2) return;
  HLTV0Cand theDitrack;
  
  for( unsigned int itrk = 0; itrk < trkcands->size(); ++itrk) 
  {
    reco::RecoChargedCandidateRef oneref(trkcands, itrk);
    
    theDitrack.tk1pt  = oneref -> pt();
    theDitrack.tk1eta = oneref -> eta();
    theDitrack.tk1phi = oneref -> phi();
  
    theDitrack.DCA    = -20  ;

    for( unsigned int jtrk = itrk+1; jtrk < trkcands->size(); ++jtrk) 
    {
      reco::RecoChargedCandidateRef tworef(trkcands, jtrk);
      //if ( tworef -> charge() * oneref -> charge() > 0) continue;
       
      theDitrack.tk2pt  = tworef -> pt();
      theDitrack.tk2eta = tworef -> eta();
      theDitrack.tk2phi = tworef -> phi();
    
      theDitrack.charge =  tworef -> charge() * oneref -> charge();
      event_.hlt_v0.push_back(theDitrack);
    } 
  }
}

void 
MuMuTrkHLT::fillHltMuVtx(const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                               const edm::Handle<reco::RecoChargedCandidateCollection> & l3cands,
                               const edm::Event                                        & event)
{

  const reco::VertexCollection & hltVertices = *hltVertexHandle.product();

  HLTMuMuVtxCand theVtx;
  for (unsigned ivtx = 0 ;  ivtx < hltVertices.size(); ivtx++){
    theVtx.x  = hltVertices.at(ivtx).position().x();
    theVtx.y  = hltVertices.at(ivtx).position().y();
    theVtx.z  = hltVertices.at(ivtx).position().z();

    theVtx.ex = hltVertices.at(ivtx).xError();
    theVtx.ey = hltVertices.at(ivtx).yError();
    theVtx.ez = hltVertices.at(ivtx).zError();
    
    theVtx.CL = TMath::Prob(hltVertices.at(ivtx).chi2(), hltVertices.at(ivtx).ndof() );
    
    bool foundMu0 = false;
    bool foundMu1 = false;

    reco::Vertex::trackRef_iterator trki;
    for (trki  = hltVertices.at(ivtx).tracks_begin(); trki != hltVertices.at(ivtx).tracks_end(); ++trki) 
    {
     reco::RecoChargedCandidateRef tmp1Ref(l3cands, (*trki).key());
     if (! (tmp1Ref -> track().isNull()) && !foundMu0 )  {
       theVtx.mu1pt = tmp1Ref -> track() -> pt();
         foundMu0 = true;
       }
       else if (! (tmp1Ref -> track().isNull()) && !foundMu1 ){
       theVtx.mu2pt = tmp1Ref -> track() -> pt();
         foundMu1 = true;
       }
    }
    event_.hlt_muvtx.push_back(theVtx);
  }
    
}

void 
MuMuTrkHLT::fillHltV0Vtx(const edm::Handle<reco::VertexCollection>               & hltVertexHandle ,
                               const edm::Handle<reco::RecoChargedCandidateCollection> & trkcands,
                               const edm::Event                                        & event)
{

  const reco::VertexCollection & hltVertices = *hltVertexHandle.product();

  HLTV0VtxCand theVtx;
  for (unsigned ivtx = 0 ;  ivtx < hltVertices.size(); ivtx++){
    theVtx.x  = hltVertices.at(ivtx).position().x();
    theVtx.y  = hltVertices.at(ivtx).position().y();
    theVtx.z  = hltVertices.at(ivtx).position().z();

    theVtx.ex = hltVertices.at(ivtx).xError();
    theVtx.ey = hltVertices.at(ivtx).yError();
    theVtx.ez = hltVertices.at(ivtx).zError();
    
    theVtx.CL = TMath::Prob(hltVertices.at(ivtx).chi2(), hltVertices.at(ivtx).ndof() );
    
    bool foundTk0 = false;
    bool foundTk1 = false;

    reco::Vertex::trackRef_iterator trki;
    for (trki  = hltVertices.at(ivtx).tracks_begin(); trki != hltVertices.at(ivtx).tracks_end(); ++trki) 
    {
     reco::RecoChargedCandidateRef tmp1Ref(trkcands, (*trki).key());
     if (! (tmp1Ref -> track().isNull()) && !foundTk0 )  {
       theVtx.tk1pt = tmp1Ref -> track() -> pt();
         foundTk0 = true;
       }
       else if (! (tmp1Ref -> track().isNull()) && !foundTk1 ){
       theVtx.tk2pt = tmp1Ref -> track() -> pt();
         foundTk1 = true;
       }
    }
    event_.hlt_v0vtx.push_back(theVtx);
  }
    
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTrkHLT::beginJob()
{
  outTree_["ntupleTree"] = outfile_-> make<TTree>("ntupleTree","ntupleTree");
  outTree_["ntupleTree"] -> Branch("event" ,&event_, 64000,2);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTrkHLT::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTrkHLT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void 
MuMuTrkHLT::beginEvent()
{
  event_.hlt.triggers.clear();
  event_.hlt.objects.clear();

  event_.hltTag.triggers.clear();
  event_.hltTag.objects.clear();

  event_.genParticles.clear();
  event_.bcands.clear();
  
  //event_.L1muons.clear();
  event_.hlt_mu.clear();
  event_.hlt_tk.clear();
  event_.hlt_dimu.clear();
  event_.hlt_v0.clear();
  event_.hlt_muvtx.clear();
  event_.hlt_v0vtx.clear();

  event_.nVtx       = -1;
  event_.trueNI     = -1;
  
  //nGoodVtx = 0; 
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTrkHLT);
