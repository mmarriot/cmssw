#include "catch2/catch_all.hpp"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  // the Run 3 ideal geometry, from the XML description, without alignment
  constexpr char const* configuration = R"(from FWCore.TestProcessor.TestProcess import *
from DPGAnalysis.MuonTools.muChamberTableProducer_cfi import muChamberTableProducer
process = TestProcess()
process.load('Configuration.Geometry.GeometryExtended2025Reco_cff')
for producer in process.es_producers_().values():
    if hasattr(producer, 'applyAlignment'):
        producer.applyAlignment = False
process.toTest = muChamberTableProducer.clone(dtSegmentSrc = 'dtSegments', cscSegmentSrc = 'none',
    gemSegmentSrc = 'none', dtRecHitSrc = 'none', cscRecHitSrc = 'cscRecHits', rpcRecHitSrc = 'rpcRecHits',
    gemRecHitSrc = 'gemRecHits')
process.moduleToTest(process.toTest)
)";

  std::vector<float> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<float>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }
}  // namespace

TEST_CASE("MuChamberTableProducer", "[MuChamberTableProducer]") {
  edm::test::TestProcessor::Config config{configuration};
  auto const dtToken = config.produces<DTRecSegment4DCollection>("dtSegments");
  auto const cscToken = config.produces<CSCRecHit2DCollection>("cscRecHits");
  auto const rpcToken = config.produces<RPCRecHitCollection>("rpcRecHits");
  auto const gemToken = config.produces<GEMRecHitCollection>("gemRecHits");
  edm::test::TestProcessor tester(config);

  // two DT chambers (one with two segments), two layers of one CSC chamber, an RPC roll, two GEM eta partitions in
  // the two layers of one superchamber
  const DTChamberId dt1(0, 1, 1), dt2(1, 2, 3);
  const CSCDetId cscLayer3(1, 2, 1, 5, 3), cscLayer4(1, 2, 1, 5, 4);
  const RPCDetId rpcRoll(0, 0, 1, 1, 1, 1, 1);
  const GEMDetId gemLayer1(1, 1, 1, 1, 10, 3), gemLayer2(1, 1, 1, 2, 10, 5);
  auto dtSegments = std::make_unique<DTRecSegment4DCollection>();
  const std::vector<DTRecSegment4D> two(2), one(1);
  dtSegments->put(dt1, two.begin(), two.end());
  dtSegments->put(dt2, one.begin(), one.end());
  auto cscHits = std::make_unique<CSCRecHit2DCollection>();
  const std::vector<CSCRecHit2D> cscHit(1);
  cscHits->put(cscLayer3, cscHit.begin(), cscHit.end());
  cscHits->put(cscLayer4, cscHit.begin(), cscHit.end());
  auto rpcHits = std::make_unique<RPCRecHitCollection>();
  const std::vector<RPCRecHit> rpcHit{RPCRecHit(rpcRoll, 0)};
  rpcHits->put(rpcRoll, rpcHit.begin(), rpcHit.end());
  auto gemHits = std::make_unique<GEMRecHitCollection>();
  const std::vector<GEMRecHit> gemHit1{GEMRecHit(gemLayer1, 0)}, gemHit2{GEMRecHit(gemLayer2, 0)};
  gemHits->put(gemLayer1, gemHit1.begin(), gemHit1.end());
  gemHits->put(gemLayer2, gemHit2.begin(), gemHit2.end());

  auto const event = tester.test(std::make_pair(dtToken, std::move(dtSegments)),
                                 std::make_pair(cscToken, std::move(cscHits)),
                                 std::make_pair(rpcToken, std::move(rpcHits)),
                                 std::make_pair(gemToken, std::move(gemHits)));
  auto const& table = *event.get<nanoaod::FlatTable>();
  REQUIRE(table.name() == "muChamber");

  // one row per chamber, ordered by DetId
  std::vector<uint32_t> expected{dt1.rawId(),
                                 dt2.rawId(),
                                 cscLayer3.chamberId().rawId(),
                                 rpcRoll.chamberId().rawId(),
                                 GEMDetId(1, 1, 1, 0, 10, 0).rawId()};
  std::sort(expected.begin(), expected.end());
  auto const detIds = table.columnData<uint32_t>(table.columnIndex("detId"));
  REQUIRE(std::vector<uint32_t>(detIds.begin(), detIds.end()) == expected);

  auto const x = column(table, "x"), y = column(table, "y"), z = column(table, "z");
  auto const ux = column(table, "ux"), uy = column(table, "uy"), uz = column(table, "uz");
  auto const vx = column(table, "vx"), vy = column(table, "vy"), vz = column(table, "vz");
  auto const bottom = column(table, "halfWidthBottom"), top = column(table, "halfWidthTop");
  auto const length = column(table, "halfLength"), thickness = column(table, "halfThickness");
  for (std::size_t i = 0; i < expected.size(); ++i) {
    INFO("chamber " << detIds[i]);
    // orthonormal local axes, positive sizes
    REQUIRE(ux[i] * ux[i] + uy[i] * uy[i] + uz[i] * uz[i] == Catch::Approx(1.f).margin(1e-5));
    REQUIRE(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i] == Catch::Approx(1.f).margin(1e-5));
    REQUIRE(ux[i] * vx[i] + uy[i] * vy[i] + uz[i] * vz[i] == Catch::Approx(0.f).margin(1e-5));
    REQUIRE(bottom[i] > 0);
    REQUIRE(top[i] > 0);
    REQUIRE(length[i] > 0);
    REQUIRE(thickness[i] > 0);
    // the centre in the right part of the detector, in cm
    const float r = std::hypot(x[i], y[i]);
    const uint32_t id = detIds[i];
    if (id == dt1.rawId() || id == rpcRoll.chamberId().rawId()) {  // MB1 and RB1 of wheel 0
      REQUIRE((r > 400 && r < 500));
      REQUIRE(std::abs(z[i]) < 130);
    } else if (id == dt2.rawId()) {  // MB2 of wheel +1
      REQUIRE((r > 480 && r < 560));
      REQUIRE((z[i] > 130 && z[i] < 400));
    } else if (id == cscLayer3.chamberId().rawId()) {  // ME+2/1
      REQUIRE((z[i] > 800 && z[i] < 900));
      REQUIRE(top[i] != Catch::Approx(bottom[i]));  // a trapezoid
    } else {                                        // GE+1/1
      REQUIRE((z[i] > 540 && z[i] < 600));
    }
    // rectangles: the DT chambers and the barrel RPC chambers
    if (id == dt1.rawId() || id == dt2.rawId() || id == rpcRoll.chamberId().rawId())
      REQUIRE(top[i] == bottom[i]);
  }
}
