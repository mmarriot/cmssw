#include "catch2/catch_all.hpp"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace {
  // the Phase-2 geometry D104 (BTL and ETL), from the XML description
  std::string configuration(std::string const& parameters = "") {
    return R"cfg(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.load('Configuration.Geometry.GeometryExtendedRun4D104Reco_cff')
process.toTest = cms.EDProducer('MTDRecHitTableProducer', src = cms.InputTag('hits'), name = cms.string('Hits'),
    doc = cms.string('the hits'))cfg" +
           parameters + R"cfg()
process.moduleToTest(process.toTest)
)cfg";
  }

  template <typename T>
  std::vector<T> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<T>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }

  nanoaod::FlatTable table(std::unique_ptr<FTLRecHitCollection> hits, std::string const& parameters = "") {
    edm::test::TestProcessor::Config config{configuration(parameters)};
    auto const token = config.produces<FTLRecHitCollection>("hits");
    edm::test::TestProcessor tester(config);
    auto const event = tester.test(std::make_pair(token, std::move(hits)));
    return *event.get<nanoaod::FlatTable>();
  }

  // DetIds of hits of a simulated event in this geometry: a BTL crystal (row 0, column 3 of its module), an ETL sensor
  const DetId btlCrystal(0x62810d43), etlSensor(0x630100b5);

  // FTLRecHit(id, row, column, energy, time, timeError, position, positionError)
  std::unique_ptr<FTLRecHitCollection> hits() {
    auto out = std::make_unique<FTLRecHitCollection>();
    out->push_back(FTLRecHit(btlCrystal, 0, 3, 4.f, 5.f, 0.25f, 0.f, 1.f));
    out->push_back(FTLRecHit(btlCrystal, 0, 3, 6.f, 7.f, 0.5f, 2.f, 1.f));  // 2 cm along the bar
    out->push_back(FTLRecHit(etlSensor, 5, 10, 0.125f, 11.f, 0.0625f, 0.f, 0.f));
    out->push_back(FTLRecHit(etlSensor, 6, 10, 0.25f, 12.f, 0.0625f, 0.f, 0.f));  // the next row
    out->push_back(FTLRecHit(etlSensor, 5, 11, 0.5f, 13.f, 0.0625f, 0.f, 0.f));   // the next column
    return out;
  }

  float distance(std::vector<float> const& x, std::vector<float> const& y, std::vector<float> const& z, int i, int j) {
    return std::hypot(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
  }
}  // namespace

TEST_CASE("MTDRecHitTableProducer: one row per hit with its position", "[MTDRecHitTableProducer]") {
  auto const out = table(hits());
  REQUIRE(out.name() == "Hits");
  REQUIRE(out.doc() == "the hits");
  REQUIRE_FALSE(out.extension());
  REQUIRE(out.nRows() == 5);
  REQUIRE(out.nColumns() == 9);
  const uint32_t btl = btlCrystal.rawId(), etl = etlSensor.rawId();
  REQUIRE(column<uint32_t>(out, "detId") == std::vector<uint32_t>{btl, btl, etl, etl, etl});
  REQUIRE(column<uint8_t>(out, "row") == std::vector<uint8_t>{0, 0, 5, 6, 5});
  REQUIRE(column<uint8_t>(out, "column") == std::vector<uint8_t>{3, 3, 10, 10, 11});
  REQUIRE(column<float>(out, "energy") == std::vector<float>{4.f, 6.f, 0.125f, 0.25f, 0.5f});
  REQUIRE(column<float>(out, "time") == std::vector<float>{5.f, 7.f, 11.f, 12.f, 13.f});
  REQUIRE(column<float>(out, "timeError") == std::vector<float>{0.25f, 0.5f, 0.0625f, 0.0625f, 0.0625f});

  auto const x = column<float>(out, "x"), y = column<float>(out, "y"), z = column<float>(out, "z");
  // BTL: the crystal centre is on the barrel cylinder; the second hit is 2 cm along the bar from it
  REQUIRE_THAT(std::hypot(x[0], y[0]), Catch::Matchers::WithinAbs(115.8, 0.3));
  REQUIRE_THAT(distance(x, y, z, 0, 1), Catch::Matchers::WithinAbs(2., 1e-4));
  // ETL: pads on a disc; the next row and the next column are one pad pitch (1.325 mm) away
  for (int i = 2; i < 5; ++i)
    REQUIRE_THAT(std::abs(z[i]), Catch::Matchers::WithinAbs(301., 2.));
  REQUIRE_THAT(distance(x, y, z, 2, 3), Catch::Matchers::WithinAbs(0.1325, 1e-4));
  REQUIRE_THAT(distance(x, y, z, 2, 4), Catch::Matchers::WithinAbs(0.1325, 1e-4));
}

TEST_CASE("MTDRecHitTableProducer: cut and precision", "[MTDRecHitTableProducer]") {
  auto const out = table(hits(), ", cut = cms.string('energy > 0.2 && time < 12.5'), precision = cms.int32(8)");
  REQUIRE(column<uint32_t>(out, "detId") ==
          std::vector<uint32_t>{btlCrystal.rawId(), btlCrystal.rawId(), etlSensor.rawId()});
  REQUIRE(column<uint8_t>(out, "row") == std::vector<uint8_t>{0, 0, 6});
  // 8 of the 23 mantissa bits are kept
  for (auto const& name : {"energy", "time", "timeError", "x", "y", "z"})
    for (float value : column<float>(out, name))
      REQUIRE((std::bit_cast<uint32_t>(value) & 0x7fff) == 0);
}

TEST_CASE("MTDRecHitTableProducer: missing input", "[MTDRecHitTableProducer]") {
  SECTION("fails by default") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration()});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("empty table with skipNonExistingSrc") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration(", skipNonExistingSrc = cms.bool(True)")});
    auto const event = tester.test();
    auto const& out = *event.get<nanoaod::FlatTable>();
    REQUIRE(out.nRows() == 0);
    REQUIRE(out.nColumns() == 9);
  }
}
