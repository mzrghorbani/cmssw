//
// $Id: LeptonJetIsolationAngle.h,v 1.1 2008/01/15 13:20:26 lowette Exp $
//

#ifndef PhysicsTools_PatUtils_LeptonJetIsolationAngle_h
#define PhysicsTools_PatUtils_LeptonJetIsolationAngle_h

/**
  \class    LeptonJetIsolationAngle LeptonJetIsolationAngle.h "PhysicsTools/PatUtils/interface/LeptonJetIsolationAngle.h"
  \brief    Calculates a lepton's jet isolation angle

   LeptonJetIsolationAngle calculates an isolation angle w.r.t. a list of
   given jets as the minimal angle to a jet in Euclidean space, as defined in
   CMS Note 2006/024

  \author   Steven Lowette
  \version  $Id: LeptonJetIsolationAngle.h,v 1.1 2008/01/15 13:20:26 lowette Exp $
*/


#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/PatUtils/interface/TrackerIsolationPt.h"


namespace pat {


  class LeptonJetIsolationAngle {

    public:

      LeptonJetIsolationAngle();
      ~LeptonJetIsolationAngle();

      float calculate(const Electron & anElectron, const edm::Handle<edm::View<reco::Track> > & trackHandle, const edm::Event & iEvent);
      float calculate(const Muon & aMuon, const edm::Handle<edm::View<reco::Track> > & trackHandle, const edm::Event & iEvent);

    private:

      float calculate(const HepLorentzVector & aLepton, const edm::Handle<edm::View<reco::Track> > & trackHandle, const edm::Event & iEvent);
      float spaceAngle(const HepLorentzVector & aLepton, const reco::CaloJet & aJet);

    private:

      TrackerIsolationPt trkIsolator_;

  };


}

#endif

