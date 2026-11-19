// Input for testTrackHitsTable.sh: tracks with their track extras and hits, and the expected rows of the hit table.
//   "" (tracks, extras, hits), pt 10, 20, 30:
//     0  a pixel hit on a two-pixel cluster of the first barrel pixel module, a missing hit, an RPC hit at the centre
//        of the first RPC roll
//     1  no hits
//     2  an RPC hit at local (10, -5, 0) cm of the first RPC roll
//   "aod": one track whose track extra is not in the event (as in AOD);
//   "Expected": a NanoAOD table of the valid hits in order, with the DetId, the layer and the position from the
//   geometry: the cluster barycentre on the module surface (the CPE adds corrections of tens of microns) and the
//   position of the RPC hits.

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <cstdint>
#include <memory>
#include <vector>

class TrackHitsTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit TrackHitsTestInputProducer(edm::ParameterSet const&)
      : clustersToken_(produces()),
        hitsToken_(produces()),
        extrasToken_(produces()),
        tracksToken_(produces()),
        missingExtrasToken_(produces("aod")),
        aodTracksToken_(produces("aod")),
        expectedToken_(produces("Expected")),
        trackerToken_(esConsumes()),
        rpcToken_(esConsumes()),
        topologyToken_(esConsumes()) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const& setup) const override {
    auto const& tracker = setup.getData(trackerToken_);
    auto const& topology = setup.getData(topologyToken_);
    auto const& pixel = dynamic_cast<PixelGeomDetUnit const&>(*tracker.detsPXB().front());
    auto const& roll = *setup.getData(rpcToken_).rolls().front();

    // a cluster of two pixels along the local x axis
    auto clusters = std::make_unique<SiPixelClusterCollectionNew>();
    {
      SiPixelClusterCollectionNew::FastFiller filler(*clusters, pixel.geographicalId());
      const uint16_t adcs[] = {20000, 20000}, rows[] = {40, 41}, columns[] = {26, 26};
      filler.push_back(SiPixelCluster(2, adcs, rows, columns, 40, 26));
    }
    auto const clusterHandle = event.put(clustersToken_, std::move(clusters));
    auto const& cluster = clusterHandle->data().front();
    const auto clusterPosition =
        pixel.surface().toGlobal(pixel.specificTopology().localPosition(MeasurementPoint(cluster.x(), cluster.y())));

    const LocalPoint rpcLocal(10, -5, 0);
    auto hits = std::make_unique<TrackingRecHitCollection>();
    hits->push_back(new SiPixelRecHit(LocalPoint(),
                                      LocalError(1e-6, 0, 1e-6),
                                      0,
                                      pixel,
                                      edmNew::makeRefTo(clusterHandle, &clusterHandle->data().front())));
    hits->push_back(new InvalidTrackingRecHit(pixel, TrackingRecHit::missing));
    hits->push_back(new RPCRecHit(roll.id(), 0, LocalPoint(0, 0, 0)));
    hits->push_back(new RPCRecHit(roll.id(), 0, rpcLocal));

    auto const hitsRef = event.getRefBeforePut(hitsToken_);
    auto const extrasRef = event.getRefBeforePut(extrasToken_);
    auto extras = std::make_unique<reco::TrackExtraCollection>(3);
    (*extras)[0].setHits(hitsRef, 0, 3);
    (*extras)[1].setHits(hitsRef, 3, 0);
    (*extras)[2].setHits(hitsRef, 3, 1);
    auto tracks = std::make_unique<reco::TrackCollection>();
    for (unsigned int i = 0; i < 3; ++i) {
      tracks->push_back(track(10. * (i + 1)));
      tracks->back().setExtra(reco::TrackExtraRef(extrasRef, i));
    }
    event.put(hitsToken_, std::move(hits));
    event.put(extrasToken_, std::move(extras));
    event.put(tracksToken_, std::move(tracks));

    // the track extras of "aod" are declared but never put
    auto aodTracks = std::make_unique<reco::TrackCollection>(1, track(10.));
    aodTracks->front().setExtra(reco::TrackExtraRef(event.getRefBeforePut(missingExtrasToken_), 0));
    event.put(aodTracksToken_, std::move(aodTracks));

    const RPCDetId rpcId = roll.id();
    const auto rpcCentre = roll.toGlobal(LocalPoint(0, 0, 0)), rpcHit = roll.toGlobal(rpcLocal);
    auto expected = std::make_unique<nanoaod::FlatTable>(3, "Expected", false);
    expected->addColumn<uint32_t>(
        "detId", std::vector<uint32_t>{pixel.geographicalId().rawId(), rpcId.rawId(), rpcId.rawId()}, "");
    expected->addColumn<uint8_t>("layer",
                                 std::vector<uint8_t>{static_cast<uint8_t>(topology.layer(pixel.geographicalId())),
                                                      static_cast<uint8_t>(rpcId.station()),
                                                      static_cast<uint8_t>(rpcId.station())},
                                 "");
    expected->addColumn<float>("x", std::vector<float>{clusterPosition.x(), rpcCentre.x(), rpcHit.x()}, "");
    expected->addColumn<float>("y", std::vector<float>{clusterPosition.y(), rpcCentre.y(), rpcHit.y()}, "");
    expected->addColumn<float>("z", std::vector<float>{clusterPosition.z(), rpcCentre.z(), rpcHit.z()}, "");
    event.put(expectedToken_, std::move(expected));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  static reco::Track track(double pt) {
    return reco::Track(1., 1., reco::Track::Point(), reco::Track::Vector(pt, 0, 0), 1, reco::Track::CovarianceMatrix());
  }

  const edm::EDPutTokenT<SiPixelClusterCollectionNew> clustersToken_;
  const edm::EDPutTokenT<TrackingRecHitCollection> hitsToken_;
  const edm::EDPutTokenT<reco::TrackExtraCollection> extrasToken_;
  const edm::EDPutTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDPutTokenT<reco::TrackExtraCollection> missingExtrasToken_;
  const edm::EDPutTokenT<reco::TrackCollection> aodTracksToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> expectedToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerToken_;
  const edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topologyToken_;
};

DEFINE_FWK_MODULE(TrackHitsTestInputProducer);
