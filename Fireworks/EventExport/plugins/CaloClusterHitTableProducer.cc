// Cluster -> hit links as a flat table: one row per (cluster, hit) pair, rows ordered by cluster.
// Hit rows are resolved through the DetId column of an existing hit table, so the indices always
// match the exported hit rows, whichever producer made that table.

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class CaloClusterHitTableProducer : public edm::global::EDProducer<> {
public:
  explicit CaloClusterHitTableProducer(const edm::ParameterSet& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        hitIdColumn_(cfg.getParameter<std::string>("hitIdColumn")),
        clustersToken_(consumes(cfg.getParameter<edm::InputTag>("clusters"))),
        hitTableToken_(consumes(cfg.getParameter<edm::InputTag>("hitTable"))) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& hits = event.get(hitTableToken_);
    const int col = hits.columnIndex(hitIdColumn_);
    if (col < 0 || hits.columnType(col) != nanoaod::FlatTable::ColumnType::UInt32)
      throw cms::Exception("Configuration")
          << "hit table '" << hits.name() << "' has no uint32 column '" << hitIdColumn_ << "'";

    // DetId -> row as a sorted vector: a single allocation and binary-search lookups.
    std::vector<std::pair<uint32_t, int32_t>> rowOf;
    rowOf.reserve(hits.size());
    int32_t row = 0;
    for (uint32_t id : hits.columnData<uint32_t>(col))
      rowOf.emplace_back(id, row++);
    std::sort(rowOf.begin(), rowOf.end());
    auto findRow = [&rowOf](uint32_t id) -> int32_t {
      auto it = std::lower_bound(rowOf.begin(), rowOf.end(), std::make_pair(id, std::numeric_limits<int32_t>::min()));
      return (it != rowOf.end() && it->first == id) ? it->second : -1;
    };

    const auto& clusters = event.get(clustersToken_);
    size_t n = 0;
    for (const auto& cluster : clusters)
      n += cluster.hitsAndFractions().size();

    std::vector<int32_t> clusterIdx, hitIdx;
    std::vector<float> fraction;
    clusterIdx.reserve(n), hitIdx.reserve(n), fraction.reserve(n);
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (const auto& [id, frac] : clusters[i].hitsAndFractions()) {
        clusterIdx.push_back(static_cast<int32_t>(i));
        hitIdx.push_back(findRow(id.rawId()));
        fraction.push_back(frac);
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(n, name_, false);
    table->setDoc(doc_);
    table->addColumn<int32_t>("clusterIdx", clusterIdx, "index of the cluster in the input collection");
    table->addColumn<int32_t>("hitIdx", hitIdx, "row in the hit table, -1 if the hit is not in it");
    table->addColumn<float>("fraction", fraction, "fraction of the hit energy assigned to the cluster");
    event.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("name")->setComment("name of the output table");
    desc.add<std::string>("doc", "")->setComment("table documentation");
    desc.add<edm::InputTag>("clusters")->setComment("any collection viewable as reco::CaloCluster");
    desc.add<edm::InputTag>("hitTable")->setComment("FlatTable holding the hits");
    desc.add<std::string>("hitIdColumn", "detId")->setComment("uint32 DetId column of the hit table");
    descriptions.addDefault(desc);
  }

private:
  const std::string name_, doc_, hitIdColumn_;
  const edm::EDGetTokenT<edm::View<reco::CaloCluster>> clustersToken_;
  const edm::EDGetTokenT<nanoaod::FlatTable> hitTableToken_;
};

DEFINE_FWK_MODULE(CaloClusterHitTableProducer);
