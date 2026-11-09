// Flat table of calorimeter rechits with their geometry: one row per hit, rows in input-collection order.
// Works for any rechit collection whose elements provide detid(), energy() and time(): ECAL, HCAL, HGCal.

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include <memory>
#include <string>
#include <vector>

namespace {
  // HGCal: global layer (CE-E first); HCAL: depth; otherwise 0.
  int16_t layerOf(DetId id, const hgcal::RecHitTools& tools) {
    switch (id.det()) {
      case DetId::HGCalEE:
      case DetId::HGCalHSi:
      case DetId::HGCalHSc:
        return static_cast<int16_t>(tools.getLayerWithOffset(id));
      case DetId::Hcal:
        return static_cast<int16_t>(HcalDetId(id).depth());
      default:
        return 0;
    }
  }

  // Approximate half-width of the cell front face [cm]: HGCal silicon apothem, HGCal scintillator
  // half r*dphi, otherwise half an edge of the cached cell's front face (ECAL crystal, HCAL tower).
  float halfWidthOf(DetId id, const GlobalPoint& pos, const hgcal::RecHitTools& tools, const CaloGeometry& geom) {
    if (tools.isSilicon(id))
      return tools.getRadiusToSide(id);
    if (tools.isScintillator(id))
      return 0.5f * pos.perp() * tools.getScintDEtaDPhi(id).second;
    const auto cell = geom.getGeometry(id);
    if (!cell)
      return 0.f;
    const auto& corners = cell->getCorners();
    return 0.5f * (corners[1] - corners[0]).mag();
  }
}  // namespace

template <typename C>
class CaloRecHitTableProducer : public edm::global::EDProducer<> {
public:
  explicit CaloRecHitTableProducer(const edm::ParameterSet& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        skipMissing_(cfg.getParameter<bool>("skipMissing")),
        srcTokens_(edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("src"),
                                         [this](const edm::InputTag& tag) { return consumes<C>(tag); })),
        geomToken_(esConsumes()) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup& setup) const override {
    const auto& geom = setup.getData(geomToken_);
    hgcal::RecHitTools tools;  // per-event setup is a few lookups; keeps the module stateless
    tools.setGeometry(geom);

    std::vector<edm::Handle<C>> inputs;
    size_t n = 0;
    for (const auto& token : srcTokens_) {
      auto handle = event.getHandle(token);
      if (!handle.isValid() && skipMissing_)
        continue;
      n += handle->size();  // dereferencing an invalid handle throws, reporting the missing product
      inputs.push_back(handle);
    }

    std::vector<uint32_t> detId;
    std::vector<uint8_t> det, subdet;
    std::vector<int16_t> layer;
    std::vector<float> energy, time, x, y, z, eta, phi, halfWidth;
    for (auto* v : {&energy, &time, &x, &y, &z, &eta, &phi, &halfWidth})
      v->reserve(n);
    detId.reserve(n), det.reserve(n), subdet.reserve(n), layer.reserve(n);

    for (const auto& input : inputs) {
      for (const auto& hit : *input) {
        const DetId id = hit.detid();
        const GlobalPoint pos = tools.getPosition(id);
        detId.push_back(id.rawId());
        det.push_back(id.det());
        subdet.push_back(id.subdetId());
        layer.push_back(layerOf(id, tools));
        energy.push_back(hit.energy());
        time.push_back(hit.time());
        x.push_back(pos.x());
        y.push_back(pos.y());
        z.push_back(pos.z());
        eta.push_back(pos.eta());
        phi.push_back(pos.phi());
        halfWidth.push_back(halfWidthOf(id, pos, tools, geom));
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(n, name_, false);
    table->setDoc(doc_);
    table->addColumn<uint32_t>("detId", detId, "raw DetId");
    table->addColumn<uint8_t>("det", det, "DetId::Detector (3 Ecal, 4 Hcal, 8 HGCalEE, 9 HGCalHSi, 10 HGCalHSc)");
    table->addColumn<uint8_t>("subdet", subdet, "DetId subdetector");
    table->addColumn<int16_t>("layer", layer, "HGCal global layer (CE-E first), HCAL depth, 0 otherwise");
    table->addColumn<float>("energy", energy, "energy [GeV]");
    table->addColumn<float>("time", time, "time [ns]");
    table->addColumn<float>("x", x, "cell centre x [cm]");
    table->addColumn<float>("y", y, "cell centre y [cm]");
    table->addColumn<float>("z", z, "cell centre z [cm]");
    table->addColumn<float>("eta", eta, "cell centre pseudorapidity");
    table->addColumn<float>("phi", phi, "cell centre azimuth [rad]");
    table->addColumn<float>("halfWidth", halfWidth, "approximate half-width of the cell front face [cm]");
    event.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("name")->setComment("name of the output table");
    desc.add<std::string>("doc", "")->setComment("table documentation");
    desc.add<std::vector<edm::InputTag>>("src")->setComment("rechit collections, concatenated in this order");
    desc.add<bool>("skipMissing", false)->setComment("ignore missing input collections");
    descriptions.addDefault(desc);
  }

private:
  const std::string name_, doc_;
  const bool skipMissing_;
  const std::vector<edm::EDGetTokenT<C>> srcTokens_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
};

using EcalRecHitTableProducer = CaloRecHitTableProducer<EcalRecHitCollection>;
using HBHERecHitTableProducer = CaloRecHitTableProducer<HBHERecHitCollection>;
using HFRecHitTableProducer = CaloRecHitTableProducer<HFRecHitCollection>;
using HORecHitTableProducer = CaloRecHitTableProducer<HORecHitCollection>;
using HGCRecHitTableProducer = CaloRecHitTableProducer<HGCRecHitCollection>;
DEFINE_FWK_MODULE(EcalRecHitTableProducer);
DEFINE_FWK_MODULE(HBHERecHitTableProducer);
DEFINE_FWK_MODULE(HFRecHitTableProducer);
DEFINE_FWK_MODULE(HORecHitTableProducer);
DEFINE_FWK_MODULE(HGCRecHitTableProducer);
