#include "catch2/catch_all.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"

#include <string>
#include <vector>

// The input is made by ConstituentTestInputProducer (this directory): candidates of pt 10, 20, 30; jets of pt 40
// (candidates 2, 0), 5 (none) and 20 (candidate 1 and a candidate of another collection); taus of pt 20 (candidate 1)
// and 40 (candidates 0, 2).
namespace {
  std::string configuration(std::string const& plugin, std::string const& jets, std::string const& parameters) {
    return R"cfg(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.input = cms.EDProducer('ConstituentTestInputProducer')
process.toTest = cms.EDProducer(')cfg" +
           plugin + R"cfg(', name = cms.string('Constituents'), candIdxName = cms.string('pfCandIdx'),
    candIdxDoc = cms.string('row in the candidate table'), candidates = cms.InputTag('input'),
    jets = cms.InputTag('input', ')cfg" +
           jets + R"cfg('), jetCut = cms.string(''), jetConstCut = cms.string(''))
)cfg" + parameters +
           R"cfg(
process.moduleToTest(process.toTest, cms.Task(process.input))
)cfg";
  }

  std::vector<int> column(nanoaod::FlatTable const& table, std::string const& name) {
    REQUIRE(table.columnIndex(name) >= 0);
    auto const data = table.columnData<int>(table.columnIndex(name));
    return {data.begin(), data.end()};
  }

  // the (parent, candidate) columns of the constituent table
  std::pair<std::vector<int>, std::vector<int>> constituents(std::string const& plugin,
                                                             std::string const& parent,
                                                             std::string const& parameters = "") {
    edm::test::TestProcessor tester(
        edm::test::TestProcessor::Config{configuration(plugin, parent == "tauIdx" ? "taus" : "jets", parameters)});
    auto const event = tester.test();
    auto const& table = *event.get<nanoaod::FlatTable>("Constituents");
    REQUIRE(table.name() == "Constituents");
    REQUIRE(table.nColumns() == 2);
    return {column(table, parent), column(table, "pfCandIdx")};
  }
}  // namespace

TEST_CASE("SimpleJetConstituentTableProducer: reco::PFJet", "[SimpleJetConstituentTableProducer]") {
  auto const plugin = "SimplePFJetConstituentTableProducer";
  SECTION("all jets: one row per constituent in the candidates, in jet and daughter order") {
    auto const [jets, candidates] = constituents(plugin, "jetIdx");
    REQUIRE(jets == std::vector<int>{0, 0, 2});
    REQUIRE(candidates == std::vector<int>{2, 0, 1});
  }
  SECTION("jet cut: the index is the row among the selected jets") {
    auto const [jets, candidates] = constituents(plugin, "jetIdx", "process.toTest.jetCut = 'pt > 10'");
    REQUIRE(jets == std::vector<int>{0, 0, 1});
    REQUIRE(candidates == std::vector<int>{2, 0, 1});
  }
  SECTION("constituent cut") {
    auto const [jets, candidates] = constituents(plugin, "jetIdx", "process.toTest.jetConstCut = 'pt > 15'");
    REQUIRE(jets == std::vector<int>{0, 2});
    REQUIRE(candidates == std::vector<int>{2, 1});
  }
}

TEST_CASE("SimpleJetConstituentTableProducer: reco::PFTau signal candidates", "[SimpleJetConstituentTableProducer]") {
  auto const plugin = "SimplePFTauConstituentTableProducer";
  SECTION("all taus") {
    auto const [taus, candidates] = constituents(plugin, "tauIdx");
    REQUIRE(taus == std::vector<int>{0, 1, 1});
    REQUIRE(candidates == std::vector<int>{1, 0, 2});
  }
  SECTION("tau cut") {
    auto const [taus, candidates] = constituents(plugin, "tauIdx", "process.toTest.jetCut = 'pt > 30'");
    REQUIRE(taus == std::vector<int>{0, 0});
    REQUIRE(candidates == std::vector<int>{0, 2});
  }
}
