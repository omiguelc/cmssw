/* 
 * class PFRecoTauDiscriminationAgainstMuon
 * created : May 07 2008,
 * revised : always,
 * Authors : Sho Maruyama,M.Bachtis
 */

#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"


#include <string>

using namespace reco;

class PFRecoTauDiscriminationAgainstMuon : public PFTauDiscriminationProducerBase 
{
 public:
  explicit PFRecoTauDiscriminationAgainstMuon(const edm::ParameterSet& iConfig)
    : PFTauDiscriminationProducerBase(iConfig) 
  {   
    discriminatorOption_  = iConfig.getParameter<std::string>("discriminatorOption");  
    hop_  = iConfig.getParameter<double>("HoPMin");  
    a  = iConfig.getParameter<double>("a");  
    b  = iConfig.getParameter<double>("b");  
    c  = iConfig.getParameter<double>("c");  	 
    maxNumberOfMatches_ = iConfig.exists("maxNumberOfMatches") ? iConfig.getParameter<int>("maxNumberOfMatches") : 0;
    checkNumMatches_ = iConfig.exists("checkNumMatches") ? iConfig.getParameter<bool>("checkNumMatches") : false;
  }

  ~PFRecoTauDiscriminationAgainstMuon() override {} 

  double discriminate(const PFTauRef& pfTau) const override;

 private:  
  std::string discriminatorOption_;
  double hop_;
  double a;
  double b;
  double c;
  int maxNumberOfMatches_;
  bool checkNumMatches_;
};

double PFRecoTauDiscriminationAgainstMuon::discriminate(const PFTauRef& thePFTauRef) const
{
  bool decision = true;

  if ( thePFTauRef->hasMuonReference() ) {
    const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(thePFTauRef->leadPFChargedHadrCand().get());
    if (pfcand != nullptr) {
      MuonRef muonref = pfcand->muonRef();
      if ( discriminatorOption_ == "noSegMatch" ) {
        if ( muonref ->numberOfMatches() > maxNumberOfMatches_ ) decision = false;
      } else if (discriminatorOption_ == "twoDCut") {
        double seg = muon::segmentCompatibility(*muonref);
        double calo= muonref->caloCompatibility(); 
        double border = calo * a + seg * b +c;
        if ( border > 0 ) decision = false; 
      } else if ( discriminatorOption_ == "merePresence" ) {
        decision = false;
      } else if (discriminatorOption_ == "combined" ) { // testing purpose only
        unsigned int muType = 0;
        if      ( muonref->isGlobalMuon()  ) muType = 1;
        else if ( muonref->isCaloMuon()    ) muType = 2;
        else if ( muonref->isTrackerMuon() ) muType = 3;
        float muonEnergyFraction = 0.;
        const reco::PFJet* pfJetPtr = dynamic_cast<const reco::PFJet*>(thePFTauRef->pfTauTagInfoRef()->pfjetRef().get());
        if (pfJetPtr) {
          muonEnergyFraction = pfJetPtr->chargedMuEnergyFraction();
        } else throw cms::Exception("Type Mismatch") << "The PFTau was not made from PFJets, and this outdated algorithm was not updated to cope with PFTaus made from other Jets.\n";

        bool eta_veto = false;
        bool phi_veto = false;
        if ( fabs(muonref->eta()) > 2.3 || (fabs(muonref->eta()) > 1.4 && fabs(muonref->eta()) < 1.6)) eta_veto = true;
        if ( muonref->phi() < 0.1 && muonref->phi() > -0.1) phi_veto = true;
        if ( muType != 1 || muonref ->numberOfMatches() > 0 || eta_veto || phi_veto || muonEnergyFraction > 0.9 ) decision = false; // as place holder
      } else if ( discriminatorOption_ == "noAllArbitrated" || discriminatorOption_ == "noAllArbitratedWithHOP" ) {
        if(checkNumMatches_ && muonref ->numberOfMatches() > maxNumberOfMatches_) decision = false;
        if ( muon::isGoodMuon(*muonref, muon::AllArbitrated) ) decision = false;      
      } else if ( discriminatorOption_ == "HOP" ) { 
        decision = true; // only calo. muon cut requested: keep all tau candidates, regardless of signals in muon system
      } else {
        throw edm::Exception(edm::errors::UnimplementedFeature) 
	  << " Invalid Discriminator option = " << discriminatorOption_ << " --> please check cfi file !!\n";
      }
    } else throw cms::Exception("Type Mismatch") << "The PFTau was not made from PFCandidates, and this outdated algorithm was not updated to cope with PFTaus made from other Candidates.\n";
  } // valid muon ref

  // Additional calo. muon cut: veto one prongs compatible with MIP signature
  if ( discriminatorOption_ == "HOP" || discriminatorOption_ == "noAllArbitratedWithHOP" ) {
    if ( thePFTauRef->leadPFChargedHadrCand().isNonnull() ) {
      const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(thePFTauRef->leadPFChargedHadrCand().get());
      if (pfcand != nullptr) {
        double muonCaloEn = pfcand->hcalEnergy() + pfcand->ecalEnergy();
        if ( thePFTauRef->decayMode() == 0 && muonCaloEn < (hop_*pfcand->p()) ) decision = false;
      }
    } else throw cms::Exception("Type Mismatch") << "The PFTau was not made from PFCandidates, and this outdated algorithm was not updated to cope with PFTaus made from other Candidates.\n";
  }

  return (decision ? 1. : 0.);
} 

DEFINE_FWK_MODULE(PFRecoTauDiscriminationAgainstMuon);
