//
// $Id: PATTauProducer.cc,v 1.8 2008/01/25 16:17:02 gpetrucc Exp $
//

#include "PhysicsTools/PatAlgos/interface/PATTauProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"

#include <DataFormats/TauReco/interface/PFTau.h>
#include <DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h>
#include <DataFormats/TauReco/interface/CaloTau.h>
#include <DataFormats/TauReco/interface/CaloTauDiscriminatorByIsolation.h>

#include "PhysicsTools/PatUtils/interface/ObjectResolutionCalc.h"
#include "PhysicsTools/PatUtils/interface/LeptonLRCalc.h"

#include <vector>
#include <memory>


using namespace pat;


PATTauProducer::PATTauProducer(const edm::ParameterSet & iConfig) {
  // initialize the configurables
  tauSrc_         = iConfig.getParameter<edm::InputTag>( "tauSource" );
  addGenMatch_    = iConfig.getParameter<bool>         ( "addGenMatch" );
  addResolutions_ = iConfig.getParameter<bool>         ( "addResolutions" );
  useNNReso_      = iConfig.getParameter<bool>         ( "useNNResolutions" );
  addLRValues_    = iConfig.getParameter<bool>         ( "addLRValues" );
  genPartSrc_     = iConfig.getParameter<edm::InputTag>( "genParticleMatch" );
  tauResoFile_    = iConfig.getParameter<std::string>  ( "tauResoFile" );
  tauLRFile_      = iConfig.getParameter<std::string>  ( "tauLRFile" ); 

  // construct resolution calculator
  if (addResolutions_) {
    theResoCalc_ = new ObjectResolutionCalc(edm::FileInPath(tauResoFile_).fullPath(), useNNReso_);
  }

  // produces vector of taus
  produces<std::vector<Tau> >();
}


PATTauProducer::~PATTauProducer() {
  if (addResolutions_) delete theResoCalc_;
}


void PATTauProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {     
  std::auto_ptr<std::vector<Tau> > patTaus(new std::vector<Tau>()); 

  edm::Handle<View<TauType> > anyTaus;
  try {
    iEvent.getByLabel(tauSrc_, anyTaus);
  } catch (const edm::Exception &e) {
    edm::LogWarning("DataSource") << "WARNING! No Tau collection found. This missing input will not block the job. Instead, an empty tau collection is being be produced.";
    iEvent.put(patTaus);
    return;
  }
   
  edm::Handle<edm::Association<reco::CandidateCollection> > genMatch;
  if (addGenMatch_) iEvent.getByLabel(genPartSrc_, genMatch); 

  // prepare LR calculation if required
  if (addLRValues_) {
    theLeptonLRCalc_ = new LeptonLRCalc(iSetup, "", "", edm::FileInPath(tauLRFile_).fullPath());
  }

  for (size_t idx = 0, ntaus = anyTaus->size(); idx < ntaus; ++idx) {
    edm::RefToBase<TauType> tausRef = anyTaus->refAt(idx);
    const TauType * originalTau = tausRef.get();

    Tau aTau(tausRef);

    if (typeid(originalTau) == typeid(const reco::PFTau *)) {
      const reco::PFTau *thePFTau = dynamic_cast<const reco::PFTau*>(originalTau);
      const reco::PFJet *pfJet    = dynamic_cast<const reco::PFJet*>(thePFTau->pfTauTagInfoRef()->pfjetRef().get());
      if(pfJet) {
        aTau.setEmEnergyFraction(pfJet->chargedEmEnergyFraction()+pfJet->neutralEmEnergyFraction());
        aTau.setEOverP(thePFTau->energy()/thePFTau->leadTrack()->p());
      }
    } else if (typeid(originalTau) == typeid(const reco::CaloTau *)) {
      const reco::CaloTau *theCaloTau = dynamic_cast<const reco::CaloTau*>(originalTau);
      const reco::CaloJet *caloJet    = dynamic_cast<const reco::CaloJet*>(theCaloTau->caloTauTagInfoRef()->calojetRef().get());
      if(caloJet) {
        aTau.setEmEnergyFraction(caloJet->emEnergyFraction());
        aTau.setEOverP(caloJet->energy()/theCaloTau->leadTrack()->p());
      }
    }

    // add MC match if demanded
    if (addGenMatch_) {
      reco::CandidateRef genTau = (*genMatch)[tausRef];
      if (genTau.isNonnull() && genTau.isAvailable() ) {
        aTau.setGenLepton(*genTau);
      } else {
        aTau.setGenLepton(reco::Particle(0, reco::Particle::LorentzVector(0,0,0,0))); // TQAF way of setting "null"
      }
    }

    // add resolution info if demanded
    if (addResolutions_) {
      (*theResoCalc_)(aTau);
    }
    // add lepton LR info if requested
    if (addLRValues_) {
      theLeptonLRCalc_->calcLikelihood(aTau, iEvent);
    }

    patTaus->push_back(aTau);
  }

  // sort taus in pT
  std::sort(patTaus->begin(), patTaus->end(), pTTauComparator_);

  // put genEvt object in Event
  iEvent.put(patTaus);

  // destroy the lepton LR calculator
  if (addLRValues_) delete theLeptonLRCalc_;

}

