// What each PF candidate (reco::PFCandidate) was made from, as rows of the tables of those objects:
// - an extension of the PF candidate table with the row of the candidate's track (trackRef), GSF track (gsfTrackRef),
//   muon (muonRef) and TICL candidate, -1 when it has none;
// - for each configured PF cluster collection, a table with one row per (candidate, cluster) pair of the PF block
//   elements the candidate was made from (elementsInBlocks), in candidate order.
// The rows are those of tables that keep the order of their collections: the candidate table with the same cut and
// maxLen, and the track, muon and cluster tables without either. A reference into another collection than the
// configured one throws, rather than giving a row of the wrong table. Reading the block elements needs the PF blocks,
// which RECO keeps and AOD does not: with no cluster collections configured, the blocks are not read. They are
// consumed, so the framework reads them before the module runs rather than when the first element is dereferenced.
// TICL candidates: the Phase-2 reconstruction makes one PF candidate per TICL candidate (pfTICL, in the order of the
// TICL candidates) and appends them to the others (particleFlowTmp merges particleFlowTmpBarrel and pfTICL, then
// PFLinker copies the list in order into particleFlow). The PF candidates themselves keep no link to them, so the
// last candidates of src must be copies of ticlCandidates, in order: this is checked on every event.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "DataFormats/TrackReco/interface/Track.h"
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

namespace {
  // The row of the object a reference points to, in a table of the collection `product` kept in order; -1 for a
  // null reference.
  template <typename R>
  int32_t rowOf(R const& ref, edm::ProductID const& product, char const* what) {
    if (ref.isNull())
      return -1;
    if (ref.id() != product)
      throw cms::Exception("Configuration") << "the reference of a PF candidate for " << what << " points into product "
                                            << ref.id() << ", not the configured " << product;
    return static_cast<int32_t>(ref.key());
  }
}  // namespace

class PFCandidateLinkTableProducer final : public edm::global::EDProducer<> {
public:
  explicit PFCandidateLinkTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        maxLen_(cfg.existsAs<unsigned int>("maxLen") ? cfg.getParameter<unsigned int>("maxLen")
                                                     : std::numeric_limits<unsigned int>::max()),
        cut_(cfg.getParameter<std::string>("cut")),
        candidatesToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        tracksToken_(optionalToken<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
        gsfTracksToken_(optionalToken<reco::GsfTrackCollection>(cfg.getParameter<edm::InputTag>("gsfTracks"))),
        muonsToken_(optionalToken<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
        ticlToken_(optionalToken<reco::PFCandidateCollection>(cfg.getParameter<edm::InputTag>("ticlCandidates"))),
        tableToken_(produces()) {
    auto const& clusters = cfg.getParameter<std::vector<edm::ParameterSet>>("clusters");
    if (!clusters.empty())
      blocksToken_ = consumes(cfg.getParameter<edm::InputTag>("blocks"));
    for (auto const& pset : clusters) {
      auto const& name = pset.getParameter<std::string>("name");
      if (name == name_)
        throw cms::Exception("Configuration")
            << "a cluster link table cannot be named as the candidate table, " << name;
      clusters_.push_back({consumes<reco::PFClusterCollection>(pset.getParameter<edm::InputTag>("src")),
                           produces<nanoaod::FlatTable>(name),
                           name,
                           pset.getParameter<std::string>("doc")});
    }
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const& candidates = event.get(candidatesToken_);
    std::vector<std::size_t> selected;  // the rows of the candidate table, as SimpleFlatTableProducer selects them
    for (std::size_t i = 0; i < candidates.size() && selected.size() < maxLen_; ++i)
      if (cut_(candidates[i]))
        selected.push_back(i);

    auto table = std::make_unique<nanoaod::FlatTable>(selected.size(), name_, false, true);
    auto addRows = [&](auto const& token, auto ref, char const* column, char const* doc) {
      if (token.isUninitialized())
        return;
      const auto product = event.getHandle(token).id();
      std::vector<int32_t> rows;
      rows.reserve(selected.size());
      for (std::size_t i : selected)
        rows.push_back(rowOf((candidates[i].*ref)(), product, column));
      table->addColumn<int32_t>(column, rows, doc);
    };
    addRows(tracksToken_, &reco::PFCandidate::trackRef, "trackIdx", "row of the track (trackRef), -1 if none");
    addRows(gsfTracksToken_,
            &reco::PFCandidate::gsfTrackRef,
            "gsfTrackIdx",
            "row of the GSF track (gsfTrackRef), -1 if none");
    addRows(muonsToken_, &reco::PFCandidate::muonRef, "muonIdx", "row of the muon (muonRef), -1 if none");
    if (!ticlToken_.isUninitialized())
      table->addColumn<int32_t>("ticlCandidateIdx",
                                ticlRows(candidates, event.get(ticlToken_), selected),
                                "row of the TICL candidate it was made from, -1 if none");
    event.put(tableToken_, std::move(table));

    if (!clusters_.empty())
      putClusterLinks(event, candidates, selected);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("particleFlow"))->setComment("PF candidates");
    desc.add<std::string>("cut", "")->setComment("selection of the candidates: must be that of the candidate table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected candidates: must be that of the candidate table");
    desc.add<std::string>("name")->setComment("name of the candidate table, extended with the columns of the rows");
    desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))
        ->setComment("the tracks of trackRef, for a trackIdx column (empty: no column)");
    desc.add<edm::InputTag>("gsfTracks", edm::InputTag("electronGsfTracks"))
        ->setComment("the GSF tracks of gsfTrackRef, for a gsfTrackIdx column (empty: no column)");
    desc.add<edm::InputTag>("muons", edm::InputTag("muons"))
        ->setComment("the muons of muonRef, for a muonIdx column (empty: no column)");
    desc.add<edm::InputTag>("ticlCandidates", edm::InputTag())
        ->setComment(
            "the PF candidates made from TICL candidates, one per TICL candidate in order, which must be the last "
            "candidates of src (pfTICL in Phase 2), for a ticlCandidateIdx column (empty: no column)");
    desc.add<edm::InputTag>("blocks", edm::InputTag("particleFlowBlock"))
        ->setComment("the PF blocks of the candidates' elements, read only when there are cluster collections");
    edm::ParameterSetDescription cluster;
    cluster.add<edm::InputTag>("src")->setComment("a PF cluster collection");
    cluster.add<std::string>("name")->setComment("name of the table of the (candidate, cluster) pairs");
    cluster.add<std::string>("doc", "")->setComment("documentation of that table");
    desc.addVPSet("clusters", cluster, {})
        ->setComment(
            "one table per PF cluster collection of the (candidate, cluster) pairs of the PF block elements; elements "
            "of other collections are left out");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  // a token for an optional input: none (uninitialized) for an empty InputTag
  template <typename T>
  edm::EDGetTokenT<T> optionalToken(edm::InputTag const& tag) {
    return tag.label().empty() ? edm::EDGetTokenT<T>() : consumes<T>(tag);
  }

  struct ClusterLinks {
    edm::EDGetTokenT<reco::PFClusterCollection> clustersToken;
    edm::EDPutTokenT<nanoaod::FlatTable> tableToken;
    std::string name, doc;
  };

  std::vector<int32_t> ticlRows(reco::PFCandidateCollection const& candidates,
                                reco::PFCandidateCollection const& ticl,
                                std::vector<std::size_t> const& selected) const {
    if (ticl.size() > candidates.size())
      throw cms::Exception("Configuration")
          << ticl.size() << " TICL PF candidates but only " << candidates.size() << " PF candidates";
    const std::size_t first = candidates.size() - ticl.size();
    for (std::size_t i = first; i < candidates.size(); ++i) {
      auto const &candidate = candidates[i], &copy = ticl[i - first];
      if (candidate.pdgId() != copy.pdgId() || candidate.charge() != copy.charge() || candidate.p4() != copy.p4())
        throw cms::Exception("Configuration") << "PF candidate " << i << " is not TICL PF candidate " << i - first
                                              << ": the PF candidates do not end with the TICL PF candidates";
    }
    std::vector<int32_t> rows;
    rows.reserve(selected.size());
    for (std::size_t i : selected)
      rows.push_back(i < first ? -1 : static_cast<int32_t>(i - first));
    return rows;
  }

  void putClusterLinks(edm::Event& event,
                       reco::PFCandidateCollection const& candidates,
                       std::vector<std::size_t> const& selected) const {
    const auto blocks = event.getHandle(blocksToken_).id();
    std::vector<edm::ProductID> products;
    products.reserve(clusters_.size());
    for (auto const& links : clusters_)
      products.push_back(event.getHandle(links.clustersToken).id());
    std::vector<std::vector<int32_t>> candidateRows(clusters_.size()), clusterRows(clusters_.size());
    for (int32_t row = 0; std::size_t i : selected) {
      for (auto const& [block, element] : candidates[i].elementsInBlocks()) {
        if (block.id() != blocks)
          throw cms::Exception("Configuration") << "a PF candidate is made of elements of PF blocks " << block.id()
                                                << ", not of the configured blocks " << blocks;
        auto const& cluster = block->elements()[element].clusterRef();
        if (cluster.isNull())
          continue;
        for (std::size_t k = 0; k < products.size(); ++k) {
          if (cluster.id() == products[k]) {
            candidateRows[k].push_back(row);
            clusterRows[k].push_back(static_cast<int32_t>(cluster.key()));
            break;
          }
        }
      }
      ++row;
    }
    for (std::size_t k = 0; k < clusters_.size(); ++k) {
      auto table = std::make_unique<nanoaod::FlatTable>(candidateRows[k].size(), clusters_[k].name, false);
      table->setDoc(clusters_[k].doc);
      table->addColumn<int32_t>("pfCandIdx", candidateRows[k], "row of the PF candidate in " + name_);
      table->addColumn<int32_t>("clusterIdx", clusterRows[k], "row of the PF cluster");
      event.put(clusters_[k].tableToken, std::move(table));
    }
  }

  const std::string name_;
  const unsigned int maxLen_;
  const StringCutObjectSelector<reco::PFCandidate> cut_;
  const edm::EDGetTokenT<reco::PFCandidateCollection> candidatesToken_;
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  const edm::EDGetTokenT<reco::PFCandidateCollection> ticlToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> tableToken_;
  edm::EDGetTokenT<reco::PFBlockCollection> blocksToken_;
  std::vector<ClusterLinks> clusters_;
};

DEFINE_FWK_MODULE(PFCandidateLinkTableProducer);
