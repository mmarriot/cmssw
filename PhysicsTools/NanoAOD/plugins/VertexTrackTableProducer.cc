// The tracks of the vertices of a reco::Vertex collection, as a table with one row per (vertex, track) pair, ordered
// by vertex and, within a vertex, in the order of Vertex::tracks():
// - vertexIdx: the row of the vertex in the vertex table, made with the same src, cut and maxLen
//   (e.g. SimpleVertexFlatTableProducer);
// - trackIdx: the index of the track in the collection `tracks`, which is its row in a table of that collection that
//   keeps its order and has no selection (e.g. SimpleTrackFlatTableProducer);
// - weight: the weight of the track in the vertex fit (Vertex::trackWeight, which is stored in steps of 1/255).
// A track can belong to more than one vertex (the primary and a secondary vertex, or several primary vertices with
// small weights), so this is a link table, in the layout of the jet constituent tables, rather than a column of the
// track table.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

class VertexTrackTableProducer final : public edm::global::EDProducer<> {
public:
  explicit VertexTrackTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        precision_(cfg.getParameter<int>("precision")),
        maxLen_(cfg.existsAs<unsigned int>("maxLen") ? cfg.getParameter<unsigned int>("maxLen")
                                                     : std::numeric_limits<unsigned int>::max()),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        cut_(cfg.getParameter<std::string>("cut")),
        verticesToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        tracksToken_(consumes(cfg.getParameter<edm::InputTag>("tracks"))),
        token_(produces()) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    // the same selection as SimpleFlatTableProducer, so that vertexIdx is a row of the vertex table
    std::vector<reco::Vertex const*> selected;
    std::size_t nLinks = 0;
    auto const vertices = event.getHandle(verticesToken_);
    if (vertices.isValid() || !skipNonExistingSrc_) {  // dereferencing an invalid handle throws
      selected.reserve(vertices->size());
      for (auto const& vertex : *vertices) {
        if (selected.size() >= maxLen_)
          break;
        if (cut_(vertex)) {
          selected.push_back(&vertex);
          nLinks += vertex.tracksSize();
        }
      }
    }

    std::vector<int32_t> vertexIdx, trackIdx;
    std::vector<float> weights;
    vertexIdx.reserve(nLinks);
    trackIdx.reserve(nLinks);
    weights.reserve(nLinks);
    if (!selected.empty()) {
      auto const tracks = event.getHandle(tracksToken_);
      const edm::ProductID tracksId = tracks.id();  // throws if the tracks are missing
      for (std::size_t row = 0; row < selected.size(); ++row) {
        for (auto const& track : selected[row]->tracks()) {
          // a track of another collection would give an index into the wrong table
          if (track.id() != tracksId)
            throw cms::Exception("Configuration") << "a track of vertex " << row << " is in the product " << track.id()
                                                  << ", not in 'tracks' (" << tracksId << ")";
          vertexIdx.push_back(static_cast<int32_t>(row));
          trackIdx.push_back(static_cast<int32_t>(track.key()));
          // trackWeight looks the track up again: linear in the tracks of the vertex, so quadratic per vertex
          weights.push_back(selected[row]->trackWeight(track));
        }
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(vertexIdx.size(), name_, false);
    table->setDoc(doc_);
    table->addColumn<int32_t>("vertexIdx", vertexIdx, "row of the vertex in the vertex table");
    table->addColumn<int32_t>("trackIdx", trackIdx, "index of the track in its collection (row in the track table)");
    table->addColumn<float>("weight", weights, "weight of the track in the vertex fit (steps of 1/255)", precision_);
    event.put(token_, std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input vertex collection, a std::vector<reco::Vertex>");
    desc.add<std::string>("cut", "")->setComment("selection of the vertices: must be that of the vertex table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected vertices: must be that of the vertex table");
    desc.add<edm::InputTag>("tracks")->setComment(
        "the collection of the vertices' tracks (reco::Track or a derived type), whose table trackIdx indexes");
    desc.add<std::string>("name")->setComment("name of the table of the (vertex, track) pairs");
    desc.add<std::string>("doc", "")->setComment("documentation of the table");
    desc.add<int>("precision", -1)->setComment("mantissa bits kept in the weight (zero or negative: all)");
    desc.add<bool>("skipNonExistingSrc", false)->setComment("an empty table for missing vertices instead of failing");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const std::string name_, doc_;
  const int precision_;
  const unsigned int maxLen_;
  const bool skipNonExistingSrc_;
  const StringCutObjectSelector<reco::Vertex> cut_;
  const edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  const edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> token_;
};

DEFINE_FWK_MODULE(VertexTrackTableProducer);
