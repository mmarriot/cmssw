#include <cstdint>
#include <optional>
#include <vector>

#include "DataFormats/Common/interface/MultiSpan.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "RecoTICL/Linking/interface/CornettoFeatures.h"
#include "RecoTICL/Linking/interface/CornettoPairTest.h"

#include "CornettoGroupsAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // The first pass of the Cornetto linking on the device (CornettoGroupsAlgo), for TracksterLinkingbyCornettoDevice:
  // the group label of each trackster of the linking, in the order of TracksterLinksProducer. The features are the
  // ones of TracksterLinkingbyCornetto, computed on the host from the trackster collections.
  class CornettoGroupsProducer : public stream::SynchronizingEDProducer<> {
  public:
    CornettoGroupsProducer(const edm::ParameterSet& config)
        : SynchronizingEDProducer(config), pairTest_(ticl::cornetto::makePairTest(config)), groupsToken_{produces()} {
      for (auto const& tag : config.getParameter<std::vector<edm::InputTag>>("tracksters_collections"))
        trackstersTokens_.emplace_back(consumes<std::vector<ticl::Trackster>>(tag));
    }

    void acquire(const device::Event& iEvent, const device::EventSetup&) override {
      // The input of TracksterLinksProducer: the collections in their order, the empty ones skipped by MultiSpan.
      edm::MultiSpan<ticl::Trackster> tracksters;
      for (auto const& token : trackstersTokens_)
        tracksters.add(iEvent.get(token));
      const uint32_t n = tracksters.size();
      labels_ = cms::alpakatools::make_host_buffer<uint32_t[]>(iEvent.queue(), n);
      if (n == 0)
        return;

      const auto f = ticl::cornetto::computeFeatures(tracksters);
      auto features = cms::alpakatools::make_host_buffer<ticl::cornetto::Feature[]>(iEvent.queue(), n);
      for (uint32_t i = 0; i < n; ++i)
        features[i] = ticl::cornetto::feature(f, i);
      auto featuresOnDevice = cms::alpakatools::make_device_buffer<ticl::cornetto::Feature[]>(iEvent.queue(), n);
      alpaka::memcpy(iEvent.queue(), featuresOnDevice, features);

      auto labelsOnDevice = cms::alpakatools::make_device_buffer<uint32_t[]>(iEvent.queue(), n);
      algo_.run(iEvent.queue(), pairTest_, featuresOnDevice.data(), n, labelsOnDevice.data());
      alpaka::memcpy(iEvent.queue(), *labels_, labelsOnDevice);
    }

    void produce(device::Event& iEvent, const device::EventSetup&) override {
      const auto n = alpaka::getExtentProduct(*labels_);
      iEvent.emplace(groupsToken_, labels_->data(), labels_->data() + n);
      labels_.reset();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::vector<edm::InputTag>>("tracksters_collections", {edm::InputTag("ticlTrackstersCLUEstering")})
          ->setComment("The trackster collections of the linking, in its order.");
      ticl::cornetto::fillPairTestDescription(desc);
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    std::vector<edm::EDGetTokenT<std::vector<ticl::Trackster>>> trackstersTokens_;
    const ticl::cornetto::PairTest pairTest_;
    const edm::EDPutTokenT<std::vector<unsigned int>> groupsToken_;
    const CornettoGroupsAlgo algo_{};
    // The labels, copied from the device in acquire and put in produce. Pinned, so that the copy does not block.
    std::optional<cms::alpakatools::host_buffer<uint32_t[]>> labels_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(CornettoGroupsProducer);
