// The tracks of each muon (reco::Muon: inner, outer, global, TPFMS, picky, DYT track) as rows of track tables, in an
// extension of the muon table: one column per (track type, track collection), the row of the muon's track of that
// type in that collection's table (kept in collection order), -1 when the muon has none or it is in another of the
// configured collections. A track type can have several collections: the outer track of a muon is the standalone
// track updated at the vertex (standAloneMuons:UpdatedAtVtx) when the update succeeded and the standalone track
// otherwise (standAloneMuons), and the two collections are not in the same order, so a string expression on the key of
// the reference cannot tell the rows apart. A track that is in none of the collections configured for its type is -1
// in every column of that type, with a warning; with strict, it throws.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {
  const std::map<std::string, reco::Muon::MuonTrackType> trackTypes{{"inner", reco::Muon::InnerTrack},
                                                                    {"outer", reco::Muon::OuterTrack},
                                                                    {"global", reco::Muon::CombinedTrack},
                                                                    {"tpfms", reco::Muon::TPFMS},
                                                                    {"picky", reco::Muon::Picky},
                                                                    {"dyt", reco::Muon::DYT}};
}  // namespace

class MuonTrackLinkTableProducer final : public edm::global::EDProducer<> {
public:
  explicit MuonTrackLinkTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        strict_(cfg.getParameter<bool>("strict")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        maxLen_(cfg.existsAs<unsigned int>("maxLen") ? cfg.getParameter<unsigned int>("maxLen")
                                                     : std::numeric_limits<unsigned int>::max()),
        cut_(cfg.getParameter<std::string>("cut")),
        muonsToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        tableToken_(produces()) {
    for (auto const& pset : cfg.getParameter<std::vector<edm::ParameterSet>>("tracks"))
      columns_.push_back({pset.getParameter<std::string>("name"),
                          pset.getParameter<std::string>("doc"),
                          trackTypes.at(pset.getParameter<std::string>("type")),
                          consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("src"))});
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const muons = event.getHandle(muonsToken_);
    std::vector<reco::Muon const*> selected;  // the rows of the muon table, as SimpleFlatTableProducer selects them
    if (muons.isValid() || !skipNonExistingSrc_) {  // dereferencing an invalid handle throws
      for (std::size_t i = 0; i < muons->size() && selected.size() < maxLen_; ++i)
        if (cut_((*muons)[i]))
          selected.push_back(&(*muons)[i]);
    }

    std::vector<edm::ProductID> products;
    products.reserve(columns_.size());
    for (auto const& column : columns_)
      products.push_back(event.getHandle(column.tracksToken).id());

    auto table = std::make_unique<nanoaod::FlatTable>(selected.size(), name_, false, true);
    std::vector<int32_t> rows(selected.size());
    for (std::size_t k = 0; k < columns_.size(); ++k) {
      for (std::size_t row = 0; row < selected.size(); ++row) {
        auto const track = selected[row]->muonTrack(columns_[k].type);
        rows[row] = track.isNonnull() && track.id() == products[k] ? static_cast<int32_t>(track.key()) : -1;
      }
      table->addColumn<int32_t>(columns_[k].name, rows, columns_[k].doc);
    }
    // every track of a configured type should be in one of the collections configured for that type: those that are
    // not are -1 in every column of the type, with a warning (strict: an error), one message per event and type
    for (std::size_t k = 0; k < columns_.size(); ++k) {
      auto const type = columns_[k].type;
      if (std::any_of(columns_.begin(), columns_.begin() + k, [type](Column const& c) { return c.type == type; }))
        continue;  // checked with the first column of the type
      std::size_t unknown = 0;
      for (reco::Muon const* muon : selected) {
        auto const track = muon->muonTrack(type);
        bool known = track.isNull();
        for (std::size_t j = k; j < columns_.size() && !known; ++j)
          known = columns_[j].type == type && products[j] == track.id();
        unknown += known ? 0 : 1;
      }
      if (unknown == 0)
        continue;
      if (strict_)
        throw cms::Exception("Configuration") << unknown << " muon tracks of the type of column " << columns_[k].name
                                              << " are in no collection configured for that type";
      edm::LogWarning("MuonTrackLinkTableProducer")
          << unknown << " muon tracks of the type of column " << columns_[k].name
          << " are in no collection configured for that type: written as -1";
    }
    event.put(tableToken_, std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("muons"))->setComment("the muons, std::vector<reco::Muon>");
    desc.add<std::string>("cut", "")->setComment("selection of the muons: must be that of the muon table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected muons: must be that of the muon table");
    desc.add<std::string>("name")->setComment("name of the muon table, extended with one column per track collection");
    desc.add<bool>("skipNonExistingSrc", false)->setComment("an empty table for missing muons instead of failing");
    desc.add<bool>("strict", false)
        ->setComment(
            "throw for a muon track in none of the collections configured for its type, instead of writing -1 with a "
            "warning");
    edm::ParameterSetDescription track;
    track.add<std::string>("name")->setComment("name of the column");
    track.add<std::string>("doc", "")->setComment("documentation of the column");
    track.ifValue(edm::ParameterDescription<std::string>("type", "inner", true, edm::Comment("the muon track")),
                  edm::allowedValues<std::string>("inner", "outer", "global", "tpfms", "picky", "dyt"));
    track.add<edm::InputTag>("src")->setComment("a track collection the muon tracks of that type are in");
    desc.addVPSet("tracks", track, {})
        ->setComment("one column per (track type, track collection); list every collection a type can be in");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  struct Column {
    std::string name, doc;
    reco::Muon::MuonTrackType type;
    edm::EDGetTokenT<reco::TrackCollection> tracksToken;
  };

  const std::string name_;
  const bool strict_, skipNonExistingSrc_;
  const unsigned int maxLen_;
  const StringCutObjectSelector<reco::Muon> cut_;
  const edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> tableToken_;
  std::vector<Column> columns_;
};

DEFINE_FWK_MODULE(MuonTrackLinkTableProducer);
