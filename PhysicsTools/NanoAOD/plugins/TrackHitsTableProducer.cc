// The valid hits of the tracks of a reco::Track collection (or of a derived type, such as reco::GsfTrack), in the
// layout of the collection tables of SimpleCollectionFlatTableProducer:
// - an extension of the track table with the number of valid hits of each track, n<hitsName>;
// - a table <hitsName> with one row per valid hit, ordered by track and, within a track, in the order of
//   Track::recHits(): the raw DetId, the global position and the layer.
// Tracker hits are stored without their local position, which is computed again from the cluster with the cluster
// parameter estimators of a TransientTrackingRecHitBuilder, without the track direction: the position the local
// reconstruction gives the hit, not the one of the final fit. A strip hit measures one coordinate only: along the
// strip its position is the strip centre. Muon hits (DT and CSC segments, RPC and GEM hits) keep their position, which
// for a segment is its position in the chamber. Hits of other types (MTD hits, and matched or projected strip hits,
// which are not on the tracks of the standard reconstruction) have NaN positions.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/trackerHitRTTI.h"
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
#include "FWCore/Utilities/interface/ESInputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonTopologies/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {
  // The local position of a tracker hit from its cluster, without the track direction; none for the hit types that
  // are not a single cluster (matched and projected strip hits) and for clusters without an estimator.
  std::optional<LocalPoint> clusterPosition(TrackingRecHit const& hit,
                                            GeomDetUnit const& det,
                                            TkTransientTrackingRecHitBuilder const& builder) {
    if (!trackerHitRTTI::isSingle(hit))
      return std::nullopt;
    auto const& cluster = static_cast<TrackerSingleRecHit const&>(hit).omniClusterRef();
    if (cluster.isPhase2() && builder.phase2TrackerClusterParameterEstimator() != nullptr)
      return builder.phase2TrackerClusterParameterEstimator()->localParameters(cluster.phase2OTCluster(), det).first;
    if (cluster.isStrip() && builder.stripClusterParameterEstimator() != nullptr)
      return builder.stripClusterParameterEstimator()->localParameters(cluster.stripCluster(), det).first;
    if (cluster.isPixel() && builder.pixelClusterParameterEstimator() != nullptr)
      return std::get<0>(builder.pixelClusterParameterEstimator()->getParameters(cluster.pixelCluster(), det));
    return std::nullopt;
  }

  // Tracker: the layer (barrel) or disk (endcap); muon system: the station; 0 otherwise.
  uint8_t layerOf(DetId id, TrackerTopology const& topology) {
    if (id.det() == DetId::Tracker)
      return topology.layer(id);
    if (id.det() != DetId::Muon)
      return 0;
    switch (id.subdetId()) {
      case MuonSubdetId::DT:
        return DTChamberId(id).station();
      case MuonSubdetId::CSC:
        return CSCDetId(id).station();
      case MuonSubdetId::RPC:
        return RPCDetId(id).station();
      case MuonSubdetId::GEM:
        return GEMDetId(id).station();
      case MuonSubdetId::ME0:
        return ME0DetId(id).station();
      default:
        return 0;
    }
  }
}  // namespace

class TrackHitsTableProducer final : public edm::global::EDProducer<> {
public:
  explicit TrackHitsTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        hitsName_(cfg.getParameter<std::string>("hitsName")),
        hitsDoc_(cfg.getParameter<std::string>("hitsDoc")),
        precision_(cfg.getParameter<int>("precision")),
        maxLen_(cfg.existsAs<unsigned int>("maxLen") ? cfg.getParameter<unsigned int>("maxLen")
                                                     : std::numeric_limits<unsigned int>::max()),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        cut_(cfg.getParameter<std::string>("cut")),
        tracksToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        // the track extras and the hits are reached through references from the tracks; consuming them (they have the
        // label of the tracks in RECO) lets the framework read them before produce() rather than on first use
        extrasToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        hitsCollectionToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        geometryToken_(esConsumes()),
        topologyToken_(esConsumes()),
        builderToken_(esConsumes(cfg.getParameter<edm::ESInputTag>("hitBuilder"))),
        countsToken_(produces()),
        hitsToken_(produces(hitsName_ + "Table")) {
    if (hitsName_ == name_)
      throw cms::Exception("Configuration") << "hitsName must differ from the name of the track table, " << name_;
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const& setup) const override {
    // the same selection as SimpleFlatTableProducer, so that the extension has the rows of the track table
    std::vector<reco::Track const*> selected;
    std::size_t nHits = 0;
    auto const tracks = event.getHandle(tracksToken_);
    if (tracks.isValid() || !skipNonExistingSrc_) {  // dereferencing an invalid handle throws
      selected.reserve(tracks->size());
      for (auto const& track : *tracks) {
        if (selected.size() >= maxLen_)
          break;
        if (cut_(track)) {
          // AOD keeps the tracks but not their extras and hits: say so rather than fail on the first dereference
          if (selected.empty() && !track.recHitsOk())
            throw cms::Exception("ProductNotFound")
                << "the hits of the tracks of src are not in the event: this table needs RECO or FEVT (track extras "
                   "and hits), not AOD";
          selected.push_back(&track);
          nHits += track.recHitsSize();  // an upper bound: the invalid hits are left out
        }
      }
    }

    auto const& geometry = setup.getData(geometryToken_);
    auto const& topology = setup.getData(topologyToken_);
    auto const* builder = dynamic_cast<TkTransientTrackingRecHitBuilder const*>(&setup.getData(builderToken_));
    if (builder == nullptr)
      throw cms::Exception("Configuration") << "hitBuilder is not a TkTransientTrackingRecHitBuilder";

    std::vector<uint16_t> counts;
    std::vector<uint32_t> detIds;
    std::vector<float> x, y, z;
    std::vector<uint8_t> layers;
    counts.reserve(selected.size());
    for (auto* column : {&x, &y, &z})
      column->reserve(nHits);
    detIds.reserve(nHits);
    layers.reserve(nHits);
    for (reco::Track const* track : selected) {
      uint16_t count = 0;  // uint16, as the counts of the collection tables of SimpleCollectionFlatTableProducer
      for (TrackingRecHit const* hit : track->recHits()) {
        if (!hit->isValid())
          continue;
        const DetId id = hit->geographicalId();
        const auto position = globalPosition(*hit, geometry, *builder);
        detIds.push_back(id.rawId());
        x.push_back(position.x());
        y.push_back(position.y());
        z.push_back(position.z());
        layers.push_back(layerOf(id, topology));
        ++count;
      }
      counts.push_back(count);
    }

    auto countsTable = std::make_unique<nanoaod::FlatTable>(counts.size(), name_, false, true);
    countsTable->addColumn<uint16_t>("n" + hitsName_, counts, "number of valid hits (rows in " + hitsName_ + ")");
    auto hitsTable = std::make_unique<nanoaod::FlatTable>(detIds.size(), hitsName_, false);
    hitsTable->setDoc(hitsDoc_);
    hitsTable->addColumn<uint32_t>(
        "detId", detIds, "raw DetId of the hit (of the chamber for a muon segment), det = detId >> 28");
    hitsTable->addColumn<float>("x", x, "global x [cm], NaN if not known", precision_);
    hitsTable->addColumn<float>("y", y, "global y [cm], NaN if not known", precision_);
    hitsTable->addColumn<float>("z", z, "global z [cm], NaN if not known", precision_);
    hitsTable->addColumn<uint8_t>(
        "layer", layers, "tracker: layer (barrel) or disk (endcap) as TrackerTopology::layer; muon: station");
    event.put(countsToken_, std::move(countsTable));
    event.put(hitsToken_, std::move(hitsTable));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment(
        "input track collection, of reco::Track or a derived type, with its track extras and hits (RECO or FEVT, not "
        "AOD)");
    desc.add<std::string>("cut", "")->setComment("selection of the tracks: must be that of the track table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected tracks: must be that of the track table");
    desc.add<std::string>("name")->setComment("name of the track table, extended with the n<hitsName> column");
    desc.add<std::string>("hitsName")
        ->setComment(
            "name of the table of the hits, one row per valid hit; n<hitsName> must not be a column of the track table "
            "already");
    desc.add<std::string>("hitsDoc", "")->setComment("documentation of the table of the hits");
    desc.add<int>("precision", -1)->setComment("mantissa bits kept in the positions (zero or negative: all)");
    desc.add<edm::ESInputTag>("hitBuilder", edm::ESInputTag("", "WithTrackAngle"))
        ->setComment(
            "TransientTrackingRecHitBuilder whose cluster parameter estimators give the tracker hit positions");
    desc.add<bool>("skipNonExistingSrc", false)->setComment("empty tables for missing tracks instead of failing");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  static GlobalPoint globalPosition(TrackingRecHit const& hit,
                                    GlobalTrackingGeometry const& geometry,
                                    TkTransientTrackingRecHitBuilder const& builder) {
    constexpr float nan = std::numeric_limits<float>::quiet_NaN();
    const DetId id = hit.geographicalId();
    if (id.det() == DetId::Tracker) {
      auto const* det = geometry.idToDetUnit(id);
      if (det != nullptr)
        if (const auto local = clusterPosition(hit, static_cast<GeomDetUnit const&>(*det), builder))
          return det->toGlobal(*local);
    } else if (id.det() == DetId::Muon) {
      if (auto const* det = geometry.idToDet(id))
        return det->toGlobal(hit.localPosition());
    }
    return GlobalPoint(nan, nan, nan);
  }

  const std::string name_, hitsName_, hitsDoc_;
  const int precision_;
  const unsigned int maxLen_;
  const bool skipNonExistingSrc_;
  const StringCutObjectSelector<reco::Track> cut_;
  const edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  const edm::EDGetTokenT<reco::TrackExtraCollection> extrasToken_;
  const edm::EDGetTokenT<TrackingRecHitCollection> hitsCollectionToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topologyToken_;
  const edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> builderToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> countsToken_, hitsToken_;
};

DEFINE_FWK_MODULE(TrackHitsTableProducer);
