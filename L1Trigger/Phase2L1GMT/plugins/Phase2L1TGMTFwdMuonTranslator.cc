// -*- C++ -*-

#ifndef PHASE2GMT_SAFWDMUONTRANSLATOR
#define PHASE2GMT_SAFWDMUONTRANSLATOR

// system include files
#include <memory>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuonPhase2/interface/SAMuon.h"
#include "DataFormats/L1TMuonPhase2/interface/EMTFTrack.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/TPUtils.h"

//
// class declaration
//
using namespace Phase2L1GMT;
using namespace l1t;

class Phase2L1TGMTFwdMuonTranslator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1TGMTFwdMuonTranslator(const edm::ParameterSet&);
  ~Phase2L1TGMTFwdMuonTranslator() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  l1t::SAMuon Convertl1tMuon(const l1t::Muon& mu, const int bx_);
  l1t::MuonStubRefVector selectLayerBX(const l1t::MuonStubRefVector& all, int bx, uint layer);
  void associateStubs(l1t::SAMuon&, const l1t::MuonStubRefVector&);

  l1t::SAMuon ConvertEMTFTrack(const l1t::phase2::EMTFTrack& track, const int bx_);

  // ----------member data ---------------------------
  edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;
  edm::EDGetTokenT<l1t::MuonStubCollection> stubToken_;

  edm::EDGetTokenT<l1t::phase2::EMTFTrackCollection> emtfTrackToken_;
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
Phase2L1TGMTFwdMuonTranslator::Phase2L1TGMTFwdMuonTranslator(const edm::ParameterSet& iConfig)
    : muonToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
      stubToken_(consumes<l1t::MuonStubCollection>(iConfig.getParameter<edm::InputTag>("stubs"))), 
      emtfTrackToken_(consumes<l1t::phase2::EMTFTrackCollection>(iConfig.getParameter<edm::InputTag>("emtfTracks"))) 
{
  produces<std::vector<l1t::SAMuon> >("prompt").setBranchAlias("prompt");
  produces<std::vector<l1t::SAMuon> >("displaced").setBranchAlias("displaced");
}

Phase2L1TGMTFwdMuonTranslator::~Phase2L1TGMTFwdMuonTranslator() {}

// ------------ method called to produce the data  ------------
void Phase2L1TGMTFwdMuonTranslator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<l1t::MuonBxCollection> muon;
  iEvent.getByToken(muonToken_, muon);
  
  edm::Handle<l1t::MuonStubCollection> stubHandle;
  iEvent.getByToken(stubToken_, stubHandle);
  
  edm::Handle<l1t::phase2::EMTFTrackCollection> emtf_tracks;
  iEvent.getByToken(emtfTrackToken_, emtf_tracks);
  
  // Process Stubs
  l1t::MuonStubRefVector stubs;
  
  for (uint i = 0; i < stubHandle->size(); ++i) {
    l1t::MuonStubRef stub(stubHandle, i);
    
    if (stub->bxNum() == 0)
      stubs.push_back(stub);
  }

  // Collect Prompt Muons
  std::vector<SAMuon> prompt;
  std::vector<SAMuon> displaced;
  //  std::vector<l1t::MuonStub> hybridStubs;
  //  edm::RefProd stubRefProd = iEvent.getRefBeforePut<std::vector<l1t::MuonStub> >("hybridStubs");
  //  l1t::MuonStubRef::key_type idxStub =0;

  for (unsigned int i = 0; i < muon->size(0); ++i) {
    const l1t::Muon& mu = muon->at(0, i);
    l1t::SAMuon samuon = Convertl1tMuon(mu, 0);
    
    if (samuon.tfType() == l1t::tftype::bmtf)
      continue;
    
    //now associate the stubs
    associateStubs(samuon, stubs);
    prompt.push_back(samuon);
  }

  for (unsigned int track_id = 0; track_id < emtf_tracks->size(); ++track_id) {
    auto& track = emtf_tracks->at(track_id);

    if (!track.valid()) {
      continue;
    }

    std::cout << "Found track " << track.unconstrained() << std::endl;

    auto samuon = ConvertEMTFTrack(track, 0);

    //now associate the stubs
    associateStubs(samuon, stubs);

    // Add To Collections
    if (track.unconstrained()) {
        displaced.push_back(samuon);
    } else {
        prompt.push_back(samuon);
    }
  }
  
  // Output Prompt Muon Collection
  std::unique_ptr<std::vector<l1t::SAMuon> > prompt_ptr = std::make_unique<std::vector<l1t::SAMuon> >(prompt);
  std::unique_ptr<std::vector<l1t::SAMuon> > displaced_ptr = std::make_unique<std::vector<l1t::SAMuon> >(displaced);
  iEvent.put(std::move(prompt_ptr), "prompt");
  iEvent.put(std::move(displaced_ptr), "displaced");
}

// ===  FUNCTION  ============================================================
//         Name:  Phase2L1TGMTFwdMuonTranslator::Convertl1tMuon
//  Description:
// ===========================================================================
SAMuon Phase2L1TGMTFwdMuonTranslator::Convertl1tMuon(const l1t::Muon& mu, const int bx_) {
  ap_uint<BITSSAQUAL> qual = mu.hwQual();
  int charge = mu.charge() > 0 ? 0 : 1;

  ap_uint<BITSPT> pt = round(mu.pt() / LSBpt);
  ap_int<BITSPHI> phi = round(mu.phi() / LSBphi);
  ap_int<BITSETA> eta = round(mu.eta() / LSBeta);
  // FIXME: Below are not well defined in phase1 GMT
  // Using the version from Correlator for now
  ap_int<BITSSAZ0> z0 = 0;  // No tracks info in Phase 1
  // Use 2 bits with LSB = 30cm for BMTF and 25cm for EMTF currently, but subjet to change
  ap_int<BITSSAD0> d0 = mu.hwDXY();

  //Here do not use the word format to GT but use the word format expected by GMT
  int bstart = 0;
  wordtype word(0);
  bstart = wordconcat<wordtype>(word, bstart, 1, 1);
  bstart = wordconcat<wordtype>(word, bstart, charge, 1);
  bstart = wordconcat<wordtype>(word, bstart, pt, BITSPT);
  bstart = wordconcat<wordtype>(word, bstart, phi, BITSPHI);
  bstart = wordconcat<wordtype>(word, bstart, eta, BITSETA);
  bstart = wordconcat<wordtype>(word, bstart, 0, BITSSAD0);
  bstart = wordconcat<wordtype>(word, bstart, qual, 8);

  SAMuon samuon(mu.p4(), charge, pt.to_uint(), eta.to_int(), phi.to_int(), z0.to_int(), d0.to_int(), qual.to_uint());
  if (mu.tfMuonIndex() >= 0 && mu.tfMuonIndex() <= 17)
    samuon.setTF(tftype::emtf_pos);
  else if (mu.tfMuonIndex() >= 18 && mu.tfMuonIndex() <= 35)
    samuon.setTF(tftype::omtf_pos);
  else if (mu.tfMuonIndex() >= 36 && mu.tfMuonIndex() <= 71)
    samuon.setTF(tftype::bmtf);
  else if (mu.tfMuonIndex() >= 72 && mu.tfMuonIndex() <= 89)
    samuon.setTF(tftype::omtf_neg);
  else
    samuon.setTF(tftype::emtf_neg);
  samuon.setWord(word);

  return samuon;
}  // -----  end of function Phase2L1TGMTFwdMuonTranslator::Convertl1tMuon  -----

l1t::MuonStubRefVector Phase2L1TGMTFwdMuonTranslator::selectLayerBX(const l1t::MuonStubRefVector& all,
                                                                    int bx,
                                                                    uint layer) {
  l1t::MuonStubRefVector out;
  for (const auto& stub : all) {
    if (stub->bxNum() == bx && stub->tfLayer() == layer)
      out.push_back(stub);
  }
  return out;
}

void Phase2L1TGMTFwdMuonTranslator::associateStubs(l1t::SAMuon& mu, const l1t::MuonStubRefVector& stubs) {
  for (unsigned int layer = 0; layer <= 4; ++layer) {
    l1t::MuonStubRefVector selectedStubs = selectLayerBX(stubs, 0, layer);
    int bestStubINT = -1;
    float dPhi = 1000.0;
    for (uint i = 0; i < selectedStubs.size(); ++i) {
      const l1t::MuonStubRef& stub = selectedStubs[i];
      float deltaPhi =
          (stub->quality() & 0x1) ? stub->offline_coord1() - mu.p4().phi() : stub->offline_coord2() - mu.p4().phi();
      if (deltaPhi > M_PI)
        deltaPhi = deltaPhi - 2 * M_PI;
      if (deltaPhi < -M_PI)
        deltaPhi = deltaPhi + 2 * M_PI;
      deltaPhi = fabs(deltaPhi);
      float deltaEta = (stub->etaQuality() == 0 || (stub->etaQuality() & 0x1))
                           ? fabs(stub->offline_eta1() - mu.p4().eta())
                           : fabs(stub->offline_eta2() - mu.p4().eta());
      if (deltaPhi < 0.3 && deltaEta < 0.3 && deltaPhi < dPhi) {
        dPhi = deltaPhi;
        bestStubINT = i;
      }
    }
    if (bestStubINT >= 0) {
      mu.addStub(selectedStubs[bestStubINT]);
    }
  }
}

SAMuon Phase2L1TGMTFwdMuonTranslator::ConvertEMTFTrack(const l1t::phase2::EMTFTrack& track, const int bx_) {
  // Convert EMTF Phi and Theta to Global Phi and Eta
  float track_phi = emtf::phase2::tp::calc_phi_glob_rad_from_loc(
        track.sector(),
        emtf::phase2::tp::calc_phi_loc_deg_from_int(track.modelPhi())
  );
  float track_theta = emtf::phase2::tp::calc_theta_rad_from_int(track.modelEta());
  float track_eta = - std::log(std::tan(track_theta / 2));

  // Calculate Lorentz Vector
  // Muon mass taken from L1Trigger/L1TMuon/plugins/L1TMuonProducer.cc
  math::PtEtaPhiMLorentzVector p4(
      track.emtfPt() * LSBpt, 
      track_eta, track_phi, 
      0.0
  );

  // Quantize Values
  ap_uint<BITSSAQUAL> qual   = 15;  // Not sure how this should be handled
  int                 charge = (track.emtfQ() == 1) ? 1 : 0; 
  ap_uint<BITSPT>     pt     = track.emtfPt();  // Already quantized by EMTF in the right units
  ap_int<BITSPHI>     phi    = round(track_phi / LSBphi);
  ap_int<BITSETA>     eta    = round(track_theta / LSBeta);
  ap_int<BITSSAZ0>    z0     = 0;  // Not implemented by EMTF yet
  ap_int<BITSSAD0>    d0     = track.emtfD0();  // Already quantized by EMTF in the right units

  //Here do not use the word format to GT but use the word format expected by GMT
  int bstart = 0;
  wordtype word(0);
  bstart = wordconcat<wordtype>(word, bstart, 1, 1);
  bstart = wordconcat<wordtype>(word, bstart, charge, 1);
  bstart = wordconcat<wordtype>(word, bstart, pt, BITSPT);
  bstart = wordconcat<wordtype>(word, bstart, phi, BITSPHI);
  bstart = wordconcat<wordtype>(word, bstart, eta, BITSETA);
  bstart = wordconcat<wordtype>(word, bstart, 0, BITSSAD0);
  bstart = wordconcat<wordtype>(word, bstart, qual, 8);

  SAMuon samuon(
          p4, 
          charge, pt.to_uint(), 
          eta.to_int(), phi.to_int(), 
          z0.to_int(), d0.to_int(), 
          qual.to_uint()
  );

  if (track.endcap() == 1)
    samuon.setTF(tftype::emtf_pos);
  else
    samuon.setTF(tftype::emtf_neg);
  
  samuon.setWord(word);

  return samuon;
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void Phase2L1TGMTFwdMuonTranslator::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void Phase2L1TGMTFwdMuonTranslator::endStream() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1TGMTFwdMuonTranslator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TGMTFwdMuonTranslator);
#endif
