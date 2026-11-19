// Input for test_catch2_MuonTrackLinkTableProducer.cc: muons with references to tracks (which needs the ProductIDs of
// the tracks, so everything is made in one module).
//   track collections "general" (3), "standAloneAtVertex" (2), "standAlone" (2), "global" (1), "other" (1)
//   "muons": 0 inner general 1, outer standAloneAtVertex 1, global global 0;
//            1 inner general 0, outer standAlone 0;
//            2 no tracks
//   "otherMuons": one muon with its outer track in "other"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <vector>

class MuonTrackLinkTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit MuonTrackLinkTestInputProducer(edm::ParameterSet const&)
      : generalToken_(produces("general")),
        atVertexToken_(produces("standAloneAtVertex")),
        standAloneToken_(produces("standAlone")),
        globalToken_(produces("global")),
        otherToken_(produces("other")),
        muonsToken_(produces("muons")),
        otherMuonsToken_(produces("otherMuons")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const general = event.emplace(generalToken_, 3);
    auto const atVertex = event.emplace(atVertexToken_, 2);
    auto const standAlone = event.emplace(standAloneToken_, 2);
    auto const global = event.emplace(globalToken_, 1);
    auto const other = event.emplace(otherToken_, 1);

    std::vector<reco::Muon> muons(3);
    muons[0].setInnerTrack(reco::TrackRef(general, 1));
    muons[0].setOuterTrack(reco::TrackRef(atVertex, 1));
    muons[0].setGlobalTrack(reco::TrackRef(global, 0));
    muons[1].setInnerTrack(reco::TrackRef(general, 0));
    muons[1].setOuterTrack(reco::TrackRef(standAlone, 0));
    event.emplace(muonsToken_, std::move(muons));

    std::vector<reco::Muon> otherMuons(1);
    otherMuons[0].setOuterTrack(reco::TrackRef(other, 0));
    event.emplace(otherMuonsToken_, std::move(otherMuons));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  const edm::EDPutTokenT<reco::TrackCollection> generalToken_, atVertexToken_, standAloneToken_, globalToken_,
      otherToken_;
  const edm::EDPutTokenT<reco::MuonCollection> muonsToken_, otherMuonsToken_;
};

DEFINE_FWK_MODULE(MuonTrackLinkTestInputProducer);
