#ifndef RecoTICL_Linking_TracksterLinkingbyCornettoDevice_h
#define RecoTICL_Linking_TracksterLinkingbyCornettoDevice_h

// Cornetto linking with the first pass (steps 1 to 3 of the grouping) from CornettoGroupsProducer, which runs it on the
// device: the group label of each input trackster, its root or itself. The second pass, the emission rule and the
// merged tracksters are the ones of TracksterLinkingbyCornetto. The groups producer reads the same trackster
// collections and has the same pair test parameters.

#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TracksterLinkingbyCornetto.h"

namespace ticl {

  class TracksterLinkingbyCornettoDevice : public TracksterLinkingbyCornetto {
  public:
    TracksterLinkingbyCornettoDevice(const edm::ParameterSet& conf,
                                     edm::ConsumesCollector iC,
                                     cms::Ort::ONNXRuntime const* onnxRuntime = nullptr)
        : TracksterLinkingbyCornetto(conf, iC, onnxRuntime),
          groupsToken_(iC.consumes<std::vector<unsigned int>>(conf.getParameter<edm::InputTag>("groups"))) {}

    void setEvent(edm::Event& iEvent, edm::EventSetup const&) override { groups_ = &iEvent.get(groupsToken_); }

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc) {
      TracksterLinkingbyCornetto::fillPSetDescription(iDesc);
      iDesc.add<edm::InputTag>("groups", edm::InputTag("ticlCornettoGroups"))
          ->setComment("The group label of each input trackster, from CornettoGroupsProducer.");
    }

  protected:
    // The groups of the labels: groups with the trackster itself as label are kept, so the groups are the ones of
    // the first pass of TracksterLinkingbyCornetto.
    std::vector<std::vector<unsigned int>> firstPassGroups(const Features& f) const override {
      if (groups_->size() != f.energy.size())
        throw edm::Exception(edm::errors::LogicError)
            << "TracksterLinkingbyCornettoDevice: " << groups_->size() << " group labels for " << f.energy.size()
            << " tracksters. The groups producer must read the trackster collections of the linking.";
      return cornetto::groups(*groups_, [](unsigned int, unsigned int) { return true; });
    }

  private:
    const edm::EDGetTokenT<std::vector<unsigned int>> groupsToken_;
    const std::vector<unsigned int>* groups_ = nullptr;
  };

}  // namespace ticl

#endif
