/** \class MuChamberTableProducer MuChamberTableProducer.cc DPGAnalysis/MuonTools/plugins/MuChamberTableProducer.cc
 *
 * EDProducer : the muon chambers that hold a rec-hit or a segment of the configured collections in the event, one
 * row per chamber, ordered by DetId, with the geometry of the chamber: its centre, its local x and y axes in global
 * coordinates, and the half sizes of its (rectangular or trapezoidal) bounds. With it, the chambers can be used as
 * volumes without the geometry, e.g. for occupancy or efficiency maps.
 * The chamber of a DT DetId is its DTChamberId, of a CSC DetId its chamber (layer 0), of an RPC DetId its chamber
 * (roll 0) and of a GEM DetId its superchamber (layer 0, eta partition 0).
 */

#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DPGAnalysis/MuonTools/interface/MuNtupleUtils.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  DetId chamberOf(DetId id) {
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        return DTChamberId(id.rawId());
      case MuonSubdetId::CSC:
        return CSCDetId(id).chamberId();
      case MuonSubdetId::RPC:
        return RPCDetId(id).chamberId();
      case MuonSubdetId::GEM: {
        const GEMDetId gem(id);
        return GEMDetId(gem.region(), gem.ring(), gem.station(), 0, gem.chamber(), 0);
      }
      default:
        throw cms::Exception("MuChamberTableProducer") << "DetId " << id.rawId() << " is not a DT, CSC, RPC or GEM one";
    }
  }

  // Appends the chambers of the DetIds of an IdToHitRange collection (whose DetIds are those of its objects).
  template <typename C>
  void addChambers(nano_mu::EDTokenHandle<C> const& token, edm::Event const& event, std::vector<uint32_t>& chambers) {
    auto const collection = token.conditionalGet(event);
    if (collection.isValid())
      for (auto id = collection->id_begin(); id != collection->id_end(); ++id)
        chambers.push_back(chamberOf(*id).rawId());
  }
}  // namespace

class MuChamberTableProducer final : public edm::global::EDProducer<> {
public:
  explicit MuChamberTableProducer(edm::ParameterSet const& config)
      : m_name{config.getParameter<std::string>("name")},
        m_precision{config.getParameter<int>("precision")},
        m_dtSegments{config, consumesCollector(), "dtSegmentSrc"},
        m_cscSegments{config, consumesCollector(), "cscSegmentSrc"},
        m_gemSegments{config, consumesCollector(), "gemSegmentSrc"},
        m_dtRecHits{config, consumesCollector(), "dtRecHitSrc"},
        m_cscRecHits{config, consumesCollector(), "cscRecHitSrc"},
        m_rpcRecHits{config, consumesCollector(), "rpcRecHitSrc"},
        m_gemRecHits{config, consumesCollector(), "gemRecHitSrc"},
        m_geometryToken{esConsumes()},
        m_token{produces()} {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const& environment) const override {
    std::vector<uint32_t> chambers;
    addChambers(m_dtSegments, event, chambers);
    addChambers(m_cscSegments, event, chambers);
    addChambers(m_gemSegments, event, chambers);
    addChambers(m_dtRecHits, event, chambers);
    addChambers(m_cscRecHits, event, chambers);
    addChambers(m_rpcRecHits, event, chambers);
    addChambers(m_gemRecHits, event, chambers);
    std::sort(chambers.begin(), chambers.end());
    chambers.erase(std::unique(chambers.begin(), chambers.end()), chambers.end());

    auto const& geometry = environment.getData(m_geometryToken);
    const auto n = chambers.size();
    std::vector<float> x(n), y(n), z(n), ux(n), uy(n), uz(n), vx(n), vy(n), vz(n);
    std::vector<float> halfWidthBottom(n), halfWidthTop(n), halfLength(n), halfThickness(n);
    for (std::size_t i = 0; i < n; ++i) {
      auto const* chamber = geometry.idToDet(DetId(chambers[i]));
      if (chamber == nullptr)
        throw cms::Exception("MuChamberTableProducer") << "no geometry for the chamber " << chambers[i];
      auto const& surface = chamber->surface();
      const auto centre = surface.position();
      const auto u = surface.toGlobal(LocalVector(1, 0, 0));
      const auto v = surface.toGlobal(LocalVector(0, 1, 0));
      x[i] = centre.x(), y[i] = centre.y(), z[i] = centre.z();
      ux[i] = u.x(), uy[i] = u.y(), uz[i] = u.z();
      vx[i] = v.x(), vy[i] = v.y(), vz[i] = v.z();
      auto const& bounds = surface.bounds();
      if (auto const* trapezoid = dynamic_cast<TrapezoidalPlaneBounds const*>(&bounds)) {
        // half bottom edge, half top edge, half thickness, half apothem
        const auto parameters = trapezoid->parameters();
        halfWidthBottom[i] = parameters[0], halfWidthTop[i] = parameters[1];
        halfThickness[i] = parameters[2], halfLength[i] = parameters[3];
      } else {
        halfWidthBottom[i] = halfWidthTop[i] = 0.5f * bounds.width();
        halfLength[i] = 0.5f * bounds.length(), halfThickness[i] = 0.5f * bounds.thickness();
      }
    }

    auto table = std::make_unique<nanoaod::FlatTable>(n, m_name, false);
    table->setDoc("muon chambers with a rec-hit or a segment, ordered by DetId (GEM: superchambers)");
    table->addColumn<uint32_t>("detId", chambers, "raw DetId of the chamber, subdetector = (detId >> 25) & 7");
    auto addFloat = [&](std::string const& name, std::vector<float> const& values, std::string const& doc) {
      table->addColumn<float>(name, values, doc, m_precision);
    };
    addFloat("x", x, "centre x in global coordinates - cm");
    addFloat("y", y, "centre y in global coordinates - cm");
    addFloat("z", z, "centre z in global coordinates - cm");
    addFloat("ux", ux, "local x axis (width) in global coordinates, x component");
    addFloat("uy", uy, "local x axis (width) in global coordinates, y component");
    addFloat("uz", uz, "local x axis (width) in global coordinates, z component");
    addFloat("vx", vx, "local y axis (length) in global coordinates, x component");
    addFloat("vy", vy, "local y axis (length) in global coordinates, y component");
    addFloat("vz", vz, "local y axis (length) in global coordinates, z component");
    addFloat("halfWidthBottom", halfWidthBottom, "half width along the local x axis at local y = -halfLength - cm");
    addFloat("halfWidthTop", halfWidthTop, "half width along the local x axis at local y = +halfLength - cm");
    addFloat("halfLength", halfLength, "half length along the local y axis - cm");
    addFloat("halfThickness", halfThickness, "half thickness along the local z axis - cm");
    event.put(m_token, std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("name", "muChamber")->setComment("name of the table");
    desc.add<int>("precision", -1)->setComment("mantissa bits kept in the floats (zero or negative: all)");
    const std::string comment{"collection whose chambers are in the table, \"none\" for none"};
    desc.add<edm::InputTag>("dtSegmentSrc", edm::InputTag{"dt4DSegments"})->setComment(comment);
    desc.add<edm::InputTag>("cscSegmentSrc", edm::InputTag{"cscSegments"})->setComment(comment);
    desc.add<edm::InputTag>("gemSegmentSrc", edm::InputTag{"gemSegments"})->setComment(comment);
    desc.add<edm::InputTag>("dtRecHitSrc", edm::InputTag{"dt1DRecHits"})->setComment(comment);
    desc.add<edm::InputTag>("cscRecHitSrc", edm::InputTag{"csc2DRecHits"})->setComment(comment);
    desc.add<edm::InputTag>("rpcRecHitSrc", edm::InputTag{"rpcRecHits"})->setComment(comment);
    desc.add<edm::InputTag>("gemRecHitSrc", edm::InputTag{"gemRecHits"})->setComment(comment);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const std::string m_name;
  const int m_precision;
  const nano_mu::EDTokenHandle<DTRecSegment4DCollection> m_dtSegments;
  const nano_mu::EDTokenHandle<CSCSegmentCollection> m_cscSegments;
  const nano_mu::EDTokenHandle<GEMSegmentCollection> m_gemSegments;
  const nano_mu::EDTokenHandle<DTRecHitCollection> m_dtRecHits;
  const nano_mu::EDTokenHandle<CSCRecHit2DCollection> m_cscRecHits;
  const nano_mu::EDTokenHandle<RPCRecHitCollection> m_rpcRecHits;
  const nano_mu::EDTokenHandle<GEMRecHitCollection> m_gemRecHits;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> m_geometryToken;
  const edm::EDPutTokenT<nanoaod::FlatTable> m_token;
};

DEFINE_FWK_MODULE(MuChamberTableProducer);
