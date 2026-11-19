// One row per rechit of an MTD rechit collection (BTL or ETL), in collection order: the raw DetId, the row and column
// of the cell in its module, the energy, time and time error, and the global position of the hit. The DetId of an ETL
// hit is that of its sensor, so a hit is identified by its DetId, row and column; the DetId of a BTL hit is that of
// its crystal.
// The position is the one the MTD clusterizer starts from (MTDThresholdClusterizer): the centre of the ETL pad, or the
// centre of the BTL crystal moved along the bar by the position measured from the time difference of its two ends
// (FTLRecHit::position, not limited to the length of the bar).

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomUtil.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class MTDRecHitTableProducer final : public edm::global::EDProducer<> {
public:
  explicit MTDRecHitTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        precision_(cfg.getParameter<int>("precision")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        cut_(cfg.getParameter<std::string>("cut")),
        hitsToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        geometryToken_(esConsumes()),
        topologyToken_(esConsumes()),
        tableToken_(produces()) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const& setup) const override {
    std::vector<uint32_t> detIds;
    std::vector<uint8_t> rows, columns;
    std::vector<float> energies, times, timeErrors, xs, ys, zs;
    auto const hits = event.getHandle(hitsToken_);
    if (hits.isValid() || !skipNonExistingSrc_) {  // dereferencing an invalid handle throws
      mtd::MTDGeomUtil geometry;
      geometry.setGeometry(&setup.getData(geometryToken_));
      geometry.setTopology(&setup.getData(topologyToken_));
      // exact without a cut, an upper bound with one
      const std::size_t n = hits->size();
      detIds.reserve(n);
      rows.reserve(n);
      columns.reserve(n);
      for (auto* column : {&energies, &times, &timeErrors, &xs, &ys, &zs})
        column->reserve(n);
      for (FTLRecHit const& hit : *hits) {
        if (!cut_(hit))
          continue;
        const DetId id = hit.id();
        const GlobalPoint position = geometry.isBTL(id)
                                         ? geometry.globalPosition(id, LocalPoint(hit.position(), 0.f, 0.f))
                                         : geometry.position(id, hit.row(), hit.column()).second;
        detIds.push_back(id.rawId());
        rows.push_back(hit.row());
        columns.push_back(hit.column());
        energies.push_back(hit.energy());
        times.push_back(hit.time());
        timeErrors.push_back(hit.timeError());
        xs.push_back(position.x());
        ys.push_back(position.y());
        zs.push_back(position.z());
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(detIds.size(), name_, false);
    table->setDoc(doc_);
    table->addColumn<uint32_t>("detId", detIds, "raw DetId of the hit: of the crystal in BTL, of the sensor in ETL");
    table->addColumn<uint8_t>("row", rows, "row of the cell in its module");
    table->addColumn<uint8_t>("column", columns, "column of the cell in its module");
    table->addColumn<float>("energy", energies, "energy [MeV]", precision_);
    table->addColumn<float>("time", times, "time [ns]", precision_);
    table->addColumn<float>("timeError", timeErrors, "uncertainty of the time [ns]", precision_);
    table->addColumn<float>("x", xs, "global x of the hit [cm]", precision_);
    table->addColumn<float>("y", ys, "global y of the hit [cm]", precision_);
    table->addColumn<float>("z", zs, "global z of the hit [cm]", precision_);
    event.put(tableToken_, std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("input MTD rechit collection (BTL or ETL)");
    desc.add<std::string>("cut", "")->setComment("selection of the hits (string expression on FTLRecHit)");
    desc.add<std::string>("name")->setComment("name of the table");
    desc.add<std::string>("doc", "")->setComment("documentation of the table");
    desc.add<int>("precision", -1)->setComment("mantissa bits kept in the float columns (zero or negative: all)");
    desc.add<bool>("skipNonExistingSrc", false)->setComment("empty table for a missing collection instead of failing");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const std::string name_, doc_;
  const int precision_;
  const bool skipNonExistingSrc_;
  const StringCutObjectSelector<FTLRecHit> cut_;
  const edm::EDGetTokenT<FTLRecHitCollection> hitsToken_;
  const edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> geometryToken_;
  const edm::ESGetToken<MTDTopology, MTDTopologyRcd> topologyToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> tableToken_;
};

DEFINE_FWK_MODULE(MTDRecHitTableProducer);
