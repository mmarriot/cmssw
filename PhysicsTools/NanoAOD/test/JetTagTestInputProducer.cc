// Input for test_catch2_JetTagTableProducer.cc: jets and jet tags keyed on them (which needs the ProductID of the
// jets, so they are made in one module).
//   "jets" 3 jets of pt 10, 20, 30; "otherJets" 3 jets
//   "a": 0.1, 0.2, 0.3 and "b": -1, 0.5, 2 on "jets"; "other": 1, 1, 1 on "otherJets"; "empty": no jets

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <vector>

class JetTagTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit JetTagTestInputProducer(edm::ParameterSet const&)
      : jetsToken_(produces("jets")),
        otherJetsToken_(produces("otherJets")),
        aToken_(produces("a")),
        bToken_(produces("b")),
        otherToken_(produces("other")),
        emptyToken_(produces("empty")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    std::vector<reco::BasicJet> jetList;
    for (double pt : {10., 20., 30.})
      jetList.emplace_back(reco::Particle::LorentzVector(pt, 0, 0, pt), reco::Particle::Point());
    auto const jets = event.emplace(jetsToken_, jetList);
    auto const otherJets = event.emplace(otherJetsToken_, jetList);
    event.emplace(aToken_, tags(jets, {0.1f, 0.2f, 0.3f}));
    event.emplace(bToken_, tags(jets, {-1.f, 0.5f, 2.f}));
    event.emplace(otherToken_, tags(otherJets, {1.f, 1.f, 1.f}));
    event.emplace(emptyToken_);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  static reco::JetTagCollection tags(edm::OrphanHandle<std::vector<reco::BasicJet>> const& jets,
                                     std::vector<float> const& values) {
    reco::JetTagCollection out{edm::RefToBaseProd<reco::Jet>(jets)};
    for (std::size_t i = 0; i < values.size(); ++i)
      out.setValue(i, values[i]);
    return out;
  }

  const edm::EDPutTokenT<std::vector<reco::BasicJet>> jetsToken_, otherJetsToken_;
  const edm::EDPutTokenT<reco::JetTagCollection> aToken_, bToken_, otherToken_, emptyToken_;
};

DEFINE_FWK_MODULE(JetTagTestInputProducer);
