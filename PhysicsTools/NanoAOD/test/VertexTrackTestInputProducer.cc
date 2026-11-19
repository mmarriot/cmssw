// Input for test_catch2_VertexTrackTableProducer.cc: vertices with references to tracks (which needs the ProductID
// of the tracks, so everything is made in one module).
//   "tracks": 3 tracks; "otherTracks": 1 track;
//   "" (the vertices), with their (track, weight):
//     0  z = -1: (0, 1), (2, 0.5)
//     1  z = 1:  (1, 0.25)
//     2  z = 2:  (2, 1), (0, 0.75)
//   "otherVertices": 1 vertex with the track of "otherTracks".

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <utility>
#include <vector>

class VertexTrackTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit VertexTrackTestInputProducer(edm::ParameterSet const&)
      : tracksToken_(produces("tracks")),
        otherTracksToken_(produces("otherTracks")),
        verticesToken_(produces()),
        otherVerticesToken_(produces("otherVertices")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const tracks = event.emplace(tracksToken_, 3);
    auto const otherTracks = event.emplace(otherTracksToken_, 1);
    auto vertex = [](double z, std::vector<std::pair<reco::TrackRef, float>> const& links) {
      reco::Vertex vertex(reco::Vertex::Point(0, 0, z), reco::Vertex::Error(), 1., 1., links.size());
      for (auto const& [track, weight] : links)
        vertex.add(reco::TrackBaseRef(track), weight);
      return vertex;
    };
    reco::VertexCollection vertices;
    vertices.push_back(vertex(-1., {{reco::TrackRef(tracks, 0), 1.f}, {reco::TrackRef(tracks, 2), 0.5f}}));
    vertices.push_back(vertex(1., {{reco::TrackRef(tracks, 1), 0.25f}}));
    vertices.push_back(vertex(2., {{reco::TrackRef(tracks, 2), 1.f}, {reco::TrackRef(tracks, 0), 0.75f}}));
    event.emplace(verticesToken_, std::move(vertices));
    event.emplace(otherVerticesToken_, 1, vertex(0., {{reco::TrackRef(otherTracks, 0), 1.f}}));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  const edm::EDPutTokenT<reco::TrackCollection> tracksToken_, otherTracksToken_;
  const edm::EDPutTokenT<reco::VertexCollection> verticesToken_, otherVerticesToken_;
};

DEFINE_FWK_MODULE(VertexTrackTestInputProducer);
