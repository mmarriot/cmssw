// Input for test_catch2_SimpleJetConstituentTableProducer.cc: PF candidates, and jets and taus whose constituents are
// Ptrs to them (which needs the ProductID of the candidates, so they are made in one module).
//   ""        std::vector<reco::PFCandidate>  pt 10, 20, 30
//   "others"  std::vector<reco::PFCandidate>  one candidate, pt 40
//   "jets"    std::vector<reco::PFJet>        pt 40: candidates 2, 0; pt 5: none; pt 20: candidate 1, others 0
//   "taus"    std::vector<reco::PFTau>        pt 20: candidate 1; pt 40: candidates 0, 2

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <memory>
#include <vector>

class ConstituentTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit ConstituentTestInputProducer(edm::ParameterSet const&)
      : candidatesToken_(produces()),
        othersToken_(produces("others")),
        jetsToken_(produces("jets")),
        tausToken_(produces("taus")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const candidates = event.emplace(candidatesToken_, makeCandidates({10, 20, 30}));
    auto const others = event.emplace(othersToken_, makeCandidates({40}));
    auto ptr = [](auto const& handle, unsigned int key) { return reco::CandidatePtr(handle, key); };

    std::vector<reco::PFJet> jets;
    jets.push_back(makeJet(40, {ptr(candidates, 2), ptr(candidates, 0)}));
    jets.push_back(makeJet(5, {}));
    jets.push_back(makeJet(20, {ptr(candidates, 1), ptr(others, 0)}));
    event.emplace(jetsToken_, std::move(jets));

    std::vector<reco::PFTau> taus(2);
    taus[0].setP4(p4(20));
    taus[0].setsignalCands({ptr(candidates, 1)});
    taus[1].setP4(p4(40));
    taus[1].setsignalCands({ptr(candidates, 0), ptr(candidates, 2)});
    event.emplace(tausToken_, std::move(taus));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  static reco::Particle::LorentzVector p4(double pt) { return {pt, 0, 0, pt}; }

  static std::vector<reco::PFCandidate> makeCandidates(std::vector<double> const& pts) {
    std::vector<reco::PFCandidate> out;
    for (double pt : pts)
      out.emplace_back(0, p4(pt), reco::PFCandidate::gamma);
    return out;
  }

  static reco::PFJet makeJet(double pt, reco::Jet::Constituents const& constituents) {
    return reco::PFJet(p4(pt), reco::Particle::Point(), reco::PFJet::Specific(), constituents);
  }

  const edm::EDPutTokenT<std::vector<reco::PFCandidate>> candidatesToken_, othersToken_;
  const edm::EDPutTokenT<std::vector<reco::PFJet>> jetsToken_;
  const edm::EDPutTokenT<std::vector<reco::PFTau>> tausToken_;
};

DEFINE_FWK_MODULE(ConstituentTestInputProducer);
