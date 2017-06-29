#ifndef  ntupleTree_h
#define  ntupleTree_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>

class GenParticleCand {
public:
  Int_t   pdgId; 
  Int_t   status; 
  Float_t energy; 
  Float_t pt; 
  Float_t eta; 
  Float_t phi; 
  
  Float_t mumPt ;
  Float_t mumEta;
  Float_t mumPhi;
  Float_t mupPt ;
  Float_t mupEta;
  Float_t mupPhi;
  Float_t pipPt  ;
  Float_t pipEta ;
  Float_t pipPhi ;
  Float_t pimPt   ;
  Float_t pimEta  ;
  Float_t pimPhi  ;

  std::vector<Int_t>  pdgMother; 
  std::vector<Int_t>  pdgRealMother; 

  GenParticleCand(){};
  virtual ~GenParticleCand(){};
  
  ClassDef(GenParticleCand,1)
};

class BCand {
public:

  Float_t  Mu1Pt  ;  
  Float_t  Mu2Pt  ;  
  Float_t  Mu1Eta ;  
  Float_t  Mu1Phi ;  
  Float_t  Mu2Eta ;  
  Float_t  Mu2Phi ;  
  Int_t    Mu1Ch  ;  
  Int_t    Mu2Ch  ;  

  Float_t MuMuMass       ;  
  Float_t MuMuCL         ;  
  Float_t JpsiPosition_x ;  
  Float_t JpsiPosition_y ;  
  Float_t JpsiPosition_z ;  
  Float_t JpsiPt         ;  
  Float_t JpsiL          ;  
  Float_t JpsiSigma      ;  
  Float_t JpsiCosBS      ;  

  Float_t TrkPPt         ;  
  Float_t TrkPEta        ;  
  Float_t TrkPPhi        ;  
  Float_t TrkPd0Sign     ;  
  Float_t TrkMPt         ;  
  Float_t TrkMEta        ;  
  Float_t TrkMPhi        ;  
  Float_t TrkMd0Sign     ;  

  Float_t V0Mass        ;  
  Float_t barV0Mass     ;  
  Float_t V0Pt          ;  
  Float_t V0Eta         ;  
  Float_t V0Phi         ;  
  Float_t V0CL          ;  
  Float_t V0L           ;  
  Float_t V0Sigma       ;  
  Float_t V0CosBS       ;  
  Float_t BMass        ;  
  //Float_t barBMass     ;  
  Float_t BPt          ;  
  Float_t BEta         ;  
  Float_t BPhi         ;  

  Float_t JpsiV0Position_x;  
  Float_t JpsiV0Position_y;  
  Float_t JpsiV0Position_z;  
  Float_t JpsiV0CL        ;  
  Float_t JpsiV0L         ;  
  Float_t JpsiV0Sigma     ;  
  Float_t JpsiV0CosBS     ;  

  BCand(){};
  virtual ~BCand(){};

  ClassDef(BCand,1)
};

class L1MuonCand {
public:

  Float_t pt;           
  Float_t eta;          
  Float_t phi;          
  Int_t   charge;      
  Int_t   quality;      
  
  L1MuonCand(){};
  virtual ~L1MuonCand(){};

  ClassDef(L1MuonCand,1)

};


class HLTMuCand {
public:

  Float_t  pt     ;  
  Float_t  eta    ;  
  Float_t  phi    ;  
  Int_t    charge ;  

  HLTMuCand(){};
  virtual ~HLTMuCand(){};

  ClassDef(HLTMuCand,1)

};


class HLTTkCand {
public:

  Float_t  pt     ;  
  Float_t  eta    ;  
  Float_t  phi    ;  
  Int_t    charge ;  
  Float_t  d0Sign ;  

  Float_t  CL      ;  
  Float_t  LSigma  ;  
  Float_t  CosBS   ;  

  HLTTkCand(){};
  virtual ~HLTTkCand(){};

  ClassDef(HLTTkCand,1)

};


class HLTDimuonCand {
public:

  Float_t  DCA     ;  
  Float_t  mu1pt   ;  
  Float_t  mu2pt   ;  
  Float_t  mu1eta  ;  
  Float_t  mu2eta  ;  
  Float_t  mu1phi  ;  
  Float_t  mu2phi  ;  
  Float_t  charge  ;  

  HLTDimuonCand(){};
  virtual ~HLTDimuonCand(){};

  ClassDef(HLTDimuonCand,1)

};


class HLTMuMuVtxCand {
public:

  Float_t  mu1pt   ;  
  Float_t  mu2pt   ;  

  Float_t  CL      ;  
  Float_t  x       ;  
  Float_t  y       ;  
  Float_t  z       ;  
  Float_t  ex      ;  
  Float_t  ey      ;  
  Float_t  ez      ; 
   
  HLTMuMuVtxCand(){};
  virtual ~HLTMuMuVtxCand(){};

  ClassDef(HLTMuMuVtxCand,1)

};


class HLTV0Cand {
public:

  Float_t  DCA     ;  
  Float_t  tk1pt   ;  
  Float_t  tk2pt   ;  
  Float_t  tk1eta  ;  
  Float_t  tk2eta  ;  
  Float_t  tk1phi  ;  
  Float_t  tk2phi  ;  
  Float_t  charge  ;  

  HLTV0Cand(){};
  virtual ~HLTV0Cand(){};

  ClassDef(HLTV0Cand,1)

};


class HLTV0VtxCand {
public:

  Float_t  tk1pt   ;  
  Float_t  tk2pt   ;  

  Float_t  CL      ;  
  Float_t  x       ;  
  Float_t  y       ;  
  Float_t  z       ;  
  Float_t  ex      ;  
  Float_t  ey      ;  
  Float_t  ez      ; 
   
  HLTV0VtxCand(){};
  virtual ~HLTV0VtxCand(){};

  ClassDef(HLTV0VtxCand,1)

};


class HLTObjCand {
public:

  std::string filterTag; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  
  HLTObjCand(){};
  virtual ~HLTObjCand(){};

  ClassDef(HLTObjCand,1)

};

class HLTInfo {
public:
  std::vector<std::string>  triggers;  
  std::vector<HLTObjCand>   objects;   

  HLTInfo(){};
  virtual ~HLTInfo(){};
  bool match( const std::string & path ) {
    if (  std::find (  triggers.begin(), triggers.end(), path ) != triggers.end() )  return true;
    //     if (! iname.compare("HLT_Mu20_v1") == 0) continue;
    return false;
  }

  bool find( const std::string & path ) {
  for ( std::vector<std::string>::const_iterator it = triggers.begin(); it != triggers.end(); ++it ) {
//    std::cout << it << std::endl;
      if ( it-> compare(path) == 0) return true;
//       if ( it->find ( path ) != std::string::npos ) return true;
  }
  return false;
  }
  
  ClassDef(HLTInfo,1)

};


class ntupleEvent {
public:

  Int_t   runNumber;             
  Int_t   luminosityBlockNumber; 
  Int_t   eventNumber;           

  Int_t   nVtx;                    
  Int_t   nTrks;   
  Float_t trueNI;   

  //Float_t bxId;
  //Float_t instLumi; 
  std::vector <GenParticleCand>     genParticles; 
  std::vector <BCand>               bcands;         

  std::vector <HLTMuCand>           hlt_mu;      
  std::vector <HLTTkCand>           hlt_tk;      
  std::vector <HLTDimuonCand>       hlt_dimu;
  std::vector <HLTV0Cand>           hlt_v0;
  std::vector <HLTMuMuVtxCand>      hlt_muvtx;      
  std::vector <HLTV0VtxCand>        hlt_v0vtx;      

  //std::vector <L1MuonCand>          L1muons; 

  HLTInfo                           hlt;
  HLTInfo                           hltTag;  

  ntupleEvent(){};
  virtual ~ntupleEvent(){};

  ClassDef(ntupleEvent,1)
};

#endif
