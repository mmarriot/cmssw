// Input for test_catch2_MuSegmentMatchTableProducer.cc: segments and muons whose matches refer to them (which needs
// the ProductIDs of the segments, so everything is made in one module).
//   "dtSegments": 2 in the chamber DT (0, 1, 1), 1 in DT (1, 2, 3); "otherDtSegments": 1 in DT (0, 1, 1);
//   "cscSegments": 2 in the chamber ME+2/1/5; "gemSegments": 1 in the superchamber GE+1/1/10;
//   "" (the muons), pt 10, 20, 30:
//     0  DT segments 1 (mask 1 << 10) and 0 (mask 0), CSC segment 1 (mask 2), GEM segment 0 (mask 4)
//     1  no matches
//     2  DT segment 2 (mask 1 << 18), CSC segment 0 (mask 1), and a segment match without any segment
//   "otherMuons": 1 muon matched to the segment of "otherDtSegments".

#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"

#include <cstdint>
#include <vector>

class MuSegmentMatchTestInputProducer : public edm::global::EDProducer<> {
public:
  explicit MuSegmentMatchTestInputProducer(edm::ParameterSet const&)
      : dtToken_(produces("dtSegments")),
        otherDtToken_(produces("otherDtSegments")),
        cscToken_(produces("cscSegments")),
        gemToken_(produces("gemSegments")),
        muonsToken_(produces()),
        otherMuonsToken_(produces("otherMuons")) {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    const DTChamberId dt1(0, 1, 1), dt2(1, 2, 3);
    const CSCDetId csc(1, 2, 1, 5, 0);
    const GEMDetId gem(1, 1, 1, 0, 10, 0);

    auto dtSegments = std::make_unique<DTRecSegment4DCollection>();
    const std::vector<DTRecSegment4D> two(2), one(1);
    dtSegments->put(dt1, two.begin(), two.end());
    dtSegments->put(dt2, one.begin(), one.end());
    auto const dt = event.put(dtToken_, std::move(dtSegments));
    auto otherDtSegments = std::make_unique<DTRecSegment4DCollection>();
    otherDtSegments->put(dt1, one.begin(), one.end());
    auto const otherDt = event.put(otherDtToken_, std::move(otherDtSegments));
    auto cscSegments = std::make_unique<CSCSegmentCollection>();
    const std::vector<CSCSegment> cscTwo(2);
    cscSegments->put(csc, cscTwo.begin(), cscTwo.end());
    auto const cscs = event.put(cscToken_, std::move(cscSegments));
    auto gemSegments = std::make_unique<GEMSegmentCollection>();
    const std::vector<GEMSegment> gemOne(1);
    gemSegments->put(gem, gemOne.begin(), gemOne.end());
    auto const gems = event.put(gemToken_, std::move(gemSegments));

    auto dtMatch = [&](DTRecSegment4DRef const& segment, uint32_t mask) {
      reco::MuonSegmentMatch match;
      match.dtSegmentRef = segment;
      match.mask = mask;
      return match;
    };
    auto cscMatch = [&](unsigned int key, uint32_t mask) {
      reco::MuonSegmentMatch match;
      match.cscSegmentRef = CSCSegmentRef(cscs, key);
      match.mask = mask;
      return match;
    };
    auto chamber = [](DetId id, std::vector<reco::MuonSegmentMatch> const& matches) {
      reco::MuonChamberMatch chamber;
      chamber.id = id;
      chamber.segmentMatches = matches;
      return chamber;
    };

    std::vector<reco::Muon> muons{muon(10), muon(20), muon(30)};
    reco::MuonSegmentMatch gemMatch;
    gemMatch.gemSegmentRef = GEMSegmentRef(gems, 0);
    gemMatch.mask = 4;
    auto gemChamber = chamber(gem, {});
    gemChamber.gemMatches.push_back(gemMatch);
    muons[0].setMatches(
        {chamber(dt1, {dtMatch(DTRecSegment4DRef(dt, 1), 1 << 10), dtMatch(DTRecSegment4DRef(dt, 0), 0)}),
         chamber(csc, {cscMatch(1, 2)}),
         gemChamber});
    muons[2].setMatches({chamber(dt2, {dtMatch(DTRecSegment4DRef(dt, 2), 1 << 18)}),
                         chamber(csc, {cscMatch(0, 1), reco::MuonSegmentMatch()})});
    event.emplace(muonsToken_, std::move(muons));

    std::vector<reco::Muon> otherMuons{muon(10)};
    otherMuons[0].setMatches({chamber(dt1, {dtMatch(DTRecSegment4DRef(otherDt, 0), 0)})});
    event.emplace(otherMuonsToken_, std::move(otherMuons));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    descriptions.addWithDefaultLabel(edm::ParameterSetDescription());
  }

private:
  static reco::Muon muon(double pt) { return reco::Muon(1, reco::Muon::LorentzVector(pt, 0, 0, pt)); }

  const edm::EDPutTokenT<DTRecSegment4DCollection> dtToken_, otherDtToken_;
  const edm::EDPutTokenT<CSCSegmentCollection> cscToken_;
  const edm::EDPutTokenT<GEMSegmentCollection> gemToken_;
  const edm::EDPutTokenT<reco::MuonCollection> muonsToken_, otherMuonsToken_;
};

DEFINE_FWK_MODULE(MuSegmentMatchTestInputProducer);
