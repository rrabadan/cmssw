// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class OniaV0Fitter : public edm::EDProducer {
  public:
      explicit OniaV0Fitter(const edm::ParameterSet&);
      ~OniaV0Fitter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

   // ----------member data ---------------------------
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> oniav0_cand_;
  double mass_;
  double massv0_;
  std::vector<double> OniaV0MassCuts_;
  std::vector<double> V0MassCuts_;
  std::vector<double> MassTraks_;
  std::string product_name_;

  template<typename T>
  struct GreaterByVProb {
     typedef T first_argument_type;
     typedef T second_argument_type;
     bool operator()( const T & t1, const T & t2 ) const {
        return t1.userFloat("vProb") > t2.userFloat("vProb");
     }   
  };  
};

OniaV0Fitter::OniaV0Fitter(const edm::ParameterSet& iConfig) {
  thebeamspot_    = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"));
  oniav0_cand_    = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaV0"));
  mass_           = iConfig.getParameter<double>("mass_constraint");
  massv0_         = iConfig.getParameter<double>("massv0_constraint");
  OniaV0MassCuts_ = iConfig.getParameter<std::vector<double>>("OniaV0MassCuts");
  V0MassCuts_     = iConfig.getParameter<std::vector<double>>("V0MassCuts");
  MassTraks_      = iConfig.getParameter<std::vector<double>>("MassTraks");
  product_name_   = iConfig.getParameter<std::string>("product_name");
  produces<pat::CompositeCandidateCollection>(product_name_);
}

OniaV0Fitter::~OniaV0Fitter() {
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void OniaV0Fitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  // Grab paramenters

  edm::Handle<pat::CompositeCandidateCollection> PsiV0CandHandle;
  iEvent.getByToken(oniav0_cand_, PsiV0CandHandle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > PsiV0CandRefitColl(new pat::CompositeCandidateCollection);

  // the beamspot 
  reco::Vertex theBeamSpotV;
  Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  int indexPsiV0=-1;
  for (pat::CompositeCandidateCollection::const_iterator oniav0=PsiV0CandHandle->begin(); oniav0!=PsiV0CandHandle->end(); ++oniav0) {
    indexPsiV0++;

    //offline V0 vtx fit
    TrackRef reftk1 = (dynamic_cast<const reco::RecoCandidate*>(oniav0->daughter("V0")->daughter("track0")))->track();
    TrackRef reftk2 = (dynamic_cast<const reco::RecoCandidate*>(oniav0->daughter("V0")->daughter("track1")))->track();

    vector<TransientTrack> V0TT;
    V0TT.push_back((*theB).build(&reftk1));
    V0TT.push_back((*theB).build(&reftk2));

    KinematicParticleFactoryFromTransientTrack pFactory;
     
    const ParticleMass muonMass(0.1056583);
    float muonSigma = muonMass*1E-6;
    const ParticleMass trackMass1(MassTraks_[0]);
    float trackSigma1 = trackMass1*1E-6;
    const ParticleMass trackMass2(MassTraks_[1]);
    float trackSigma2 = trackMass2*1E-6;
    
    std::vector<RefCountedKinematicParticle> allV0Daughters;
    allV0Daughters.push_back(pFactory.particle (V0TT[0], trackMass1, float(0), float(0), trackSigma1));
    allV0Daughters.push_back(pFactory.particle (V0TT[1], trackMass2, float(0), float(0), trackSigma2));

    KinematicParticleVertexFitter vertexfitter;
    RefCountedKinematicTree V0Tree = vertexfitter.fit(allV0Daughters);
    
    if (V0Tree->isEmpty()) continue;
    
    V0Tree->movePointerToTheTop();
    RefCountedKinematicParticle fitV0 = V0Tree->currentParticle();
    RefCountedKinematicVertex   V0DecayVertex = V0Tree->currentDecayVertex();
    double v0_ma_fit = 14000.;
    double v0_x2_fit = 10000.;
    double v0_vp_fit = -9999.;

    if (fitV0->currentState().isValid()) { 
      v0_ma_fit = fitV0->currentState().mass();
      v0_x2_fit = V0DecayVertex->chiSquared();
      v0_vp_fit = ChiSquaredProbability(v0_x2_fit, (double)(V0DecayVertex->degreesOfFreedom()));

      if ( v0_ma_fit < V0MassCuts_[0] || v0_ma_fit > V0MassCuts_[1] || v0_vp_fit < 0.01 ) continue;

      int    v0_ch_fit = (oniav0->daughter("V0"))->charge();
      double v0_px_fit = fitV0->currentState().kinematicParameters().momentum().x();
      double v0_py_fit = fitV0->currentState().kinematicParameters().momentum().y();
      double v0_pz_fit = fitV0->currentState().kinematicParameters().momentum().z();
      double v0_p_fit  = fitV0->currentState().globalMomentum().mag();
      double v0_en_fit = sqrt(v0_ma_fit*v0_ma_fit+v0_px_fit*v0_px_fit+v0_py_fit*v0_py_fit+v0_pz_fit*v0_pz_fit);

      double v0_vx_fit = V0DecayVertex->position().x();
      double v0_vy_fit = V0DecayVertex->position().y();
      double v0_vz_fit = V0DecayVertex->position().z();

      TVector3 vtx;
      TVector3 pvtx;
      VertexDistanceXY vdistXY;

      //-------------V0 info
      vtx.SetXYZ(v0_vx_fit,v0_vy_fit,0);
      TVector3 pperp(v0_px_fit, v0_py_fit, 0);
      AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
      pvtx.SetXYZ(theBeamSpotV .position().x(),theBeamSpotV .position().y(),0);
      TVector3 vdiff = vtx - pvtx;
      double cosAlphav0 = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
      Measurement1D distXY = vdistXY.distance(reco::Vertex(*V0DecayVertex), theBeamSpotV);
      double ctauv0 = distXY.value()*cosAlphav0 * v0_ma_fit/pperp.Perp();
      GlobalError v1e = (reco::Vertex(*V0DecayVertex)).error();
      GlobalError v2e = theBeamSpotV .error();
      AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
      double ctauErrv0 = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*v0_ma_fit/(pperp.Perp2());

      //Constrained vtx fit for V0
      KinematicConstrainedVertexFitter constVertexFitter;
      MultiTrackKinematicConstraint *v0_mc = new  TwoTrackMassKinematicConstraint(massv0_);
      RefCountedKinematicTree constrainedV0Tree = constVertexFitter.fit(allV0Daughters,v0_mc);
      
      if (constrainedV0Tree->isEmpty()) continue;
      
      constrainedV0Tree->movePointerToTheTop();
      RefCountedKinematicParticle constfitV0 = constrainedV0Tree->currentParticle();
      
      if (constfitV0->currentState().isValid()) {

        const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(oniav0->daughter("onia")->daughter("muon1") );
        const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(oniav0->daughter("onia")->daughter("muon2") );

        //TrackRef refmu1 = (dynamic_cast<const pat::Muon*>(oniav0->daughter("onia")->daughter("muon1")))->innerTrack();
        //TrackRef refmu2 = (dynamic_cast<const pat::Muon*>(oniav0->daughter("onia")->daughter("muon1")))->innerTrack();

        TrackRef refmu1 = pmu1->innerTrack();
        TrackRef refmu2 = pmu2->innerTrack();

        vector<TransientTrack> MuMu;
        MuMu.push_back((*theB).build(&refmu1));
        MuMu.push_back((*theB).build(&refmu2));

        std::vector<RefCountedKinematicParticle> allPsiV0Daughters;
        allPsiV0Daughters.push_back(pFactory.particle (MuMu[0], muonMass, float(0), float(0), muonSigma));
        allPsiV0Daughters.push_back(pFactory.particle (MuMu[1], muonMass, float(0), float(0), muonSigma));
        allPsiV0Daughters.push_back(constfitV0);

        MultiTrackKinematicConstraint *onia_mtc = new  TwoTrackMassKinematicConstraint(mass_);
        RefCountedKinematicTree PsiV0Tree = constVertexFitter.fit(allPsiV0Daughters,onia_mtc);

        if (PsiV0Tree->isEmpty()) continue;
          
        PsiV0Tree->movePointerToTheTop();
        RefCountedKinematicParticle fitPsiV0 = PsiV0Tree->currentParticle();
        RefCountedKinematicVertex PsiV0DecayVertex = PsiV0Tree->currentDecayVertex();
        
        // Get PsiV0 reffited
        double oniav0_ma_fit = 14000.;
        double oniav0_x2_fit = 10000.;
        double oniav0_vp_fit = -9999.;
        if ( fitPsiV0->currentState().isValid()){

          if ( oniav0_ma_fit > OniaV0MassCuts_[0] && oniav0_ma_fit < OniaV0MassCuts_[1] && oniav0_vp_fit > 0.01 ) {
            int    oniav0_ch_fit = oniav0->charge();         
            double oniav0_px_fit = fitPsiV0->currentState().kinematicParameters().momentum().x();
            double oniav0_py_fit = fitPsiV0->currentState().kinematicParameters().momentum().y();
            double oniav0_pz_fit = fitPsiV0->currentState().kinematicParameters().momentum().z();
            double oniav0_en_fit = sqrt(oniav0_ma_fit*oniav0_ma_fit+oniav0_px_fit*oniav0_px_fit+
                                        oniav0_py_fit*oniav0_py_fit+oniav0_pz_fit*oniav0_pz_fit);
            double oniav0_vx_fit = PsiV0DecayVertex->position().x();
            double oniav0_vy_fit = PsiV0DecayVertex->position().y();
            double oniav0_vz_fit = PsiV0DecayVertex->position().z();

            vtx.SetXYZ(oniav0_vx_fit,oniav0_vy_fit,0);
            pperp.SetX(v0_px_fit); pperp.SetY(v0_py_fit); 
            AlgebraicVector3 vpperpp(pperp.x(),pperp.y(),0);
            vdiff = vtx - pvtx;
            double cosAlphab = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            distXY = vdistXY.distance(reco::Vertex(*PsiV0DecayVertex), theBeamSpotV);
            double ctaub = distXY.value()*cosAlphab * oniav0_ma_fit/pperp.Perp();
            v1e = (reco::Vertex(*PsiV0DecayVertex)).error();
            vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrb = sqrt(ROOT::Math::Similarity(vpperpp,vXYe)) * oniav0_ma_fit/(pperp.Perp2());

            reco::CompositeCandidate recoPsiV0(oniav0_ch_fit,math::XYZTLorentzVector(oniav0_px_fit,oniav0_py_fit,
                                            oniav0_pz_fit,oniav0_en_fit),math::XYZPoint(oniav0_vx_fit,oniav0_vy_fit,oniav0_vz_fit),511);
            pat::CompositeCandidate patPsiV0(recoPsiV0);
            patPsiV0.addUserFloat("vProb",oniav0_vp_fit);
            patPsiV0.addUserFloat("vChi2",oniav0_x2_fit);
            patPsiV0.addUserFloat("cosAlpha",cosAlphab);
            patPsiV0.addUserFloat("ctauPV",ctaub);
            patPsiV0.addUserFloat("ctauErrPV",ctauErrb);
            patPsiV0.addUserInt("bIndex",indexPsiV0); 

// get first muon
            bool child = PsiV0Tree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMu1 = PsiV0Tree->currentParticle();
            if (!child) break;
            float m1_ma_fit = fitMu1->currentState().mass();
            int   m1_ch_fit = fitMu1->currentState().particleCharge();
            float m1_px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
            float m1_py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
            float m1_pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
            float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
            reco::CompositeCandidate recoMu1(m1_ch_fit,math::XYZTLorentzVector(m1_px_fit,m1_py_fit,m1_pz_fit,m1_en_fit),
                                             math::XYZPoint(oniav0_vx_fit,oniav0_vy_fit,oniav0_vz_fit),13);
            pat::CompositeCandidate patMu1(recoMu1);
// get second muon
            child = PsiV0Tree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitMu2 = PsiV0Tree->currentParticle();
            if (!child) break;
            float m2_ma_fit = fitMu2->currentState().mass();
            int   m2_ch_fit = fitMu2->currentState().particleCharge();
            float m2_px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
            float m2_py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
            float m2_pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
            float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
            reco::CompositeCandidate recoMu2(m2_ch_fit,math::XYZTLorentzVector(m2_px_fit,m2_py_fit,m2_pz_fit,m2_en_fit),
                                             math::XYZPoint(oniav0_vx_fit,oniav0_vy_fit,oniav0_vz_fit),13);
            pat::CompositeCandidate patMu2(recoMu2);

            //muon Quality
            int muontag1 = 0;
            int muontag2 = 0;
            if ( pmu1->isGlobalMuon() )            muontag1 += 1;
            if ( pmu2->isGlobalMuon() )            muontag2 += 1;
            if ( pmu1->isSoftMuon(theBeamSpotV) )  muontag1 += 2;
            if ( pmu2->isSoftMuon(theBeamSpotV) )  muontag2 += 2;
            if ( pmu1->isTightMuon(theBeamSpotV) ) muontag1 += 4;
            if ( pmu2->isTightMuon(theBeamSpotV) ) muontag2 += 4;

            patMu1.addUserInt("muontag",muontag1);
            patMu2.addUserInt("muontag",muontag2);
            
            //Define psi from two muons
            pat::CompositeCandidate patPsi;
            patPsi.addDaughter(patMu1,"muon1");
            patPsi.addDaughter(patMu2,"muon2");  
            patPsi.setP4(patMu1.p4()+patMu2.p4());

            const pat::CompositeCandidate *dimuonC = dynamic_cast<const pat::CompositeCandidate *>(oniav0->daughter("onia"));
            patPsi.addUserFloat("vProb",    dimuonC->userFloat("vProb"));
            patPsi.addUserFloat("cosAlpha", dimuonC->userFloat("cosAlpha"));
            patPsi.addUserFloat("ppdlPV",   dimuonC->userFloat("ppdlPV"));
            patPsi.addUserFloat("ppdlErrPV",dimuonC->userFloat("ppdlErrPV"));

            child = PsiV0Tree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitV0 = PsiV0Tree->currentParticle();
            if (!child) break;
            reco::CompositeCandidate recoV0(v0_ch_fit,math::XYZTLorentzVector(v0_px_fit,v0_py_fit,v0_pz_fit,v0_en_fit),
                                             math::XYZPoint(oniav0_vx_fit,oniav0_vy_fit,oniav0_vz_fit),13);
            pat::CompositeCandidate patV0(recoV0);
            patV0.addUserFloat("vProb",     v0_vp_fit);
            patV0.addUserFloat("cosAlpha",  cosAlphav0);
            patV0.addUserFloat("ctau",      ctauv0);
            patV0.addUserFloat("ctauErr",   ctauErrv0);
// get tk1
            child = V0Tree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitTrk1 = V0Tree->currentParticle();
            if (!child) break;
            float tk_ma_fit = fitTrk1->currentState().mass();
            int   tk_ch_fit = fitTrk1->currentState().particleCharge();
            float tk_px_fit = fitTrk1->currentState().kinematicParameters().momentum().x();
            float tk_py_fit = fitTrk1->currentState().kinematicParameters().momentum().y();
            float tk_pz_fit = fitTrk1->currentState().kinematicParameters().momentum().z();
            float tk_en_fit = sqrt(tk_ma_fit*tk_ma_fit+tk_px_fit*tk_px_fit+tk_py_fit*tk_py_fit+tk_pz_fit*tk_pz_fit);
            reco::CompositeCandidate recoTrk1(tk_ch_fit,math::XYZTLorentzVector(tk_px_fit,tk_py_fit,tk_pz_fit,tk_en_fit),
                                             math::XYZPoint(v0_vx_fit,v0_vy_fit,v0_vz_fit));
            pat::CompositeCandidate patTrk1(recoTrk1);

// get tk2
            child = V0Tree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitTrk2 = V0Tree->currentParticle();
            if (!child) break;
            float tk2_ma_fit = fitTrk2->currentState().mass();
            int   tk2_ch_fit = fitTrk2->currentState().particleCharge();
            float tk2_px_fit = fitTrk2->currentState().kinematicParameters().momentum().x();
            float tk2_py_fit = fitTrk2->currentState().kinematicParameters().momentum().y();
            float tk2_pz_fit = fitTrk2->currentState().kinematicParameters().momentum().z();
            float tk2_en_fit = sqrt(tk2_ma_fit*tk2_ma_fit+tk2_px_fit*tk2_px_fit+tk2_py_fit*tk2_py_fit+tk2_pz_fit*tk2_pz_fit);
            reco::CompositeCandidate recoTrk2(tk2_ch_fit,math::XYZTLorentzVector(tk2_px_fit,tk2_py_fit,tk2_pz_fit,tk2_en_fit),
                                             math::XYZPoint(v0_vx_fit,v0_vy_fit,v0_vz_fit));
            pat::CompositeCandidate patTrk2(recoTrk2);

            patV0.addDaughter(patTrk1, "track1");
            patV0.addDaughter(patTrk2, "track2");
            //patV0.setP4(patTrk1.p4()+patTrk2.p4());

            patPsiV0.addDaughter(patPsi,"onia");
            patPsiV0.addDaughter(patV0,"ditrak");
            
            PsiV0CandRefitColl->push_back(patPsiV0);
          }
        }
      }
    }
  }
  // now sort by vProb
  OniaV0Fitter::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(PsiV0CandRefitColl->begin(),PsiV0CandRefitColl->end(),vPComparator);
  iEvent.put(std::move(PsiV0CandRefitColl));
}

// ------------ method called once each job just before starting event loop  ------------
void OniaV0Fitter::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void OniaV0Fitter::endJob() {}

// ------------ method called when starting to processes a run  ------------
void OniaV0Fitter::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void OniaV0Fitter::endRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void OniaV0Fitter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void OniaV0Fitter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaV0Fitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaV0Fitter);
