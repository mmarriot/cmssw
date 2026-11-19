#include "catch2/catch_all.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <cmath>
#include <string>
#include <vector>

// The input is made by JetTagTestInputProducer (this directory): jets of pt 10, 20, 30 with tags a = 0.1, 0.2, 0.3
// and b = -1, 0.5, 2; tags "other" of other jets, and an empty JetTagCollection "empty".
namespace {
  std::string configuration(std::string const& parameters = "") {
    return R"cfg(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.input = cms.EDProducer('JetTagTestInputProducer')
process.toTest = cms.EDProducer('JetTagTableProducer', src = cms.InputTag('input', 'jets'), name = cms.string('Jet'),
    variables = cms.PSet(
        tagA = cms.PSet(src = cms.VInputTag('input:a'), doc = cms.string('tag a')),
        tagB = cms.PSet(src = cms.VInputTag('input:b')),
        tagAB = cms.PSet(src = cms.VInputTag('input:a', 'input:b'), precision = cms.int32(6))))
)cfg" + parameters +
           R"cfg(
process.moduleToTest(process.toTest, cms.Task(process.input))
)cfg";
  }

  std::vector<float> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<float>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }
}  // namespace

TEST_CASE("JetTagTableProducer: tag values and sums as columns", "[JetTagTableProducer]") {
  edm::test::TestProcessor tester(edm::test::TestProcessor::Config{configuration()});
  auto const event = tester.test();
  auto const& table = *event.get<nanoaod::FlatTable>();
  REQUIRE(table.name() == "Jet");
  REQUIRE(table.extension());
  REQUIRE(table.nRows() == 3);
  REQUIRE(table.nColumns() == 3);
  REQUIRE(table.columnDoc(table.columnIndex("tagA")) == "tag a");
  REQUIRE(column(table, "tagA") == std::vector<float>{0.1f, 0.2f, 0.3f});
  REQUIRE(column(table, "tagB") == std::vector<float>{-1.f, 0.5f, 2.f});  // negative values as they are
  auto const sum = column(table, "tagAB");
  REQUIRE(sum.size() == 3);
  REQUIRE_THAT(sum[0], Catch::Matchers::WithinRel(-0.9f, 0.02f));  // 6 mantissa bits
  REQUIRE_THAT(sum[1], Catch::Matchers::WithinRel(0.7f, 0.02f));
  REQUIRE_THAT(sum[2], Catch::Matchers::WithinRel(2.3f, 0.02f));
}

TEST_CASE("JetTagTableProducer: the selection of the jet table", "[JetTagTableProducer]") {
  SECTION("cut") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration("process.toTest.cut = cms.string('pt > 15')")});
    auto const event = tester.test();
    REQUIRE(column(*event.get<nanoaod::FlatTable>(), "tagA") == std::vector<float>{0.2f, 0.3f});
  }
  SECTION("maxLen") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration("process.toTest.maxLen = cms.uint32(1)")});
    auto const event = tester.test();
    REQUIRE(column(*event.get<nanoaod::FlatTable>(), "tagB") == std::vector<float>{-1.f});
  }
}

TEST_CASE("JetTagTableProducer: tags of other jets, and no tags", "[JetTagTableProducer]") {
  SECTION("tags of another jet collection throw") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration("process.toTest.variables.tagA.src = ['input:other']")});
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
  }
  SECTION("an empty JetTagCollection gives NaN") {
    edm::test::TestProcessor tester(edm::test::TestProcessor::Config{
        configuration("process.toTest.variables.tagAB.src = ['input:a', 'input:empty']")});
    auto const event = tester.test();
    for (float value : column(*event.get<nanoaod::FlatTable>(), "tagAB"))
      REQUIRE(std::isnan(value));
  }
}
