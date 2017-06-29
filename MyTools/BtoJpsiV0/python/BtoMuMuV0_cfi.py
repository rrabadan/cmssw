import FWCore.ParameterSet.Config as cms

psiV0producer = cms.EDProducer('OniaV0Producer',
    Onia            = cms.InputTag('onia2MuMuPAT'),
    VCol            = cms.InputTag('oniaV0Tracks','Kshort'),
    thebeamspot_    = cms.InputTag('offlineBeamSpot'),
    OniaMassCuts    = cms.vdouble(2.946916,3.246916),      # J/psi mass window 3.096916 +/- 0.150
    V0MassCuts      = cms.vdouble(0.43,0.57),  # phi mass window 1.019461 +/- .015
    OniaV0MassCuts  = cms.vdouble(5.0,5.7),            # b-hadron mass window
    MassTraks       = cms.vdouble(0.1395706,0.1395706),         # traks masses
    hits            = cms.double(7),
    ptcut           = cms.double(0.7),
    chi2cut         = cms.double(10.),
    IPsigXYcut      = cms.double(2.0),
    OnlyBest        = cms.bool(True)    
)

psiV0Fitter = cms.EDProducer('OniaV0Fitter',
    OniaV0                = cms.InputTag('OniaV0Producer','OniaV0Candidates'),
    thebeamspot_          = cms.InputTag('offlineBeamSpot'),
    mass_constraint       = cms.double(3.096916),              # J/psi mass in GeV
    massv0_constraint     = cms.double(0.497611), 
    OniaV0MassCuts        = cms.vdouble(5.0,5.7),            # b-hadron mass window
    V0MassCuts            = cms.vdouble(0.43,0.57),
    MassTraks             = cms.vdouble(0.1395706,0.1395706),         # traks masses
    product_name          = cms.string('OniaV0Candidates')
)
