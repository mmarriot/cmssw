/** \class MuSegmentMatchTableProducer MuSegmentMatchTableProducers.cc DPGAnalysis/MuonTools/plugins/MuSegmentMatchTableProducers.cc
 *
 * EDProducers : the muon <-> segment link tables, one per segment type (DT, CSC, GEM), from the segment matches
 * that the muon reconstruction stores in each muon (reco::Muon::matches()). One row per (muon, matched segment):
 * - muonIdx: the row of the muon in the muon table, made with the same src, cut and maxLen;
 * - segmentIdx: the index of the segment in its collection, which is its row in the segment table of that collection
 *   (DTSegmentFlatTableProducer, CSCSegmentFlatTableProducer, GEMSegmentFlatTableProducer keep the collection order);
 * - mask: the arbitration bits of the match (reco::MuonSegmentMatch::mask).
 * A segment can be matched to several muons, and a muon to several segments in a chamber, so this is a link table.
 */

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {
  // Where the matches of each segment type are in a reco::MuonChamberMatch, and their reference.
  struct DTMatches {
    using Collection = DTRecSegment4DCollection;
    static constexpr char const* kind = "DT";
    static constexpr char const* segments = "dt4DSegments";
    static auto const& matches(reco::MuonChamberMatch const& chamber) { return chamber.segmentMatches; }
    static auto const& ref(reco::MuonSegmentMatch const& match) { return match.dtSegmentRef; }
  };
  struct CSCMatches {
    using Collection = CSCSegmentCollection;
    static constexpr char const* kind = "CSC";
    static constexpr char const* segments = "cscSegments";
    static auto const& matches(reco::MuonChamberMatch const& chamber) { return chamber.segmentMatches; }
    static auto const& ref(reco::MuonSegmentMatch const& match) { return match.cscSegmentRef; }
  };
  struct GEMMatches {
    using Collection = GEMSegmentCollection;
    static constexpr char const* kind = "GEM";
    static constexpr char const* segments = "gemSegments";
    static auto const& matches(reco::MuonChamberMatch const& chamber) { return chamber.gemMatches; }
    static auto const& ref(reco::MuonSegmentMatch const& match) { return match.gemSegmentRef; }
  };
}  // namespace

template <typename M>
class MuSegmentMatchTableProducer final : public edm::global::EDProducer<> {
public:
  explicit MuSegmentMatchTableProducer(edm::ParameterSet const& config)
      : m_name{config.getParameter<std::string>("name")},
        m_doc{config.getParameter<std::string>("doc")},
        m_maxLen{config.existsAs<unsigned int>("maxLen") ? config.getParameter<unsigned int>("maxLen")
                                                         : std::numeric_limits<unsigned int>::max()},
        m_skipNonExistingSrc{config.getParameter<bool>("skipNonExistingSrc")},
        m_cut{config.getParameter<std::string>("cut")},
        m_muonToken{consumes(config.getParameter<edm::InputTag>("src"))},
        m_segmentToken{consumes(config.getParameter<edm::InputTag>("segmentSrc"))},
        m_token{produces()} {}

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    std::vector<int32_t> muonIdx, segmentIdx;
    std::vector<uint32_t> masks;
    auto const muons = event.getHandle(m_muonToken);
    auto const segments = event.getHandle(m_segmentToken);
    if (m_skipNonExistingSrc && !(muons.isValid() && segments.isValid())) {
      event.put(m_token, table(muonIdx, segmentIdx, masks));
      return;
    }

    const edm::ProductID segmentsId = segments.id();  // throws if the segments are missing
    int32_t row = 0;  // the same selection as SimpleFlatTableProducer, so that muonIdx is a row of the muon table
    for (auto const& muon : *muons) {
      if (static_cast<unsigned int>(row) >= m_maxLen)
        break;
      if (!m_cut(muon))
        continue;
      for (auto const& chamber : muon.matches())
        for (auto const& match : M::matches(chamber)) {
          auto const& segment = M::ref(match);
          if (segment.isNull())  // a match of another segment type
            continue;
          // a segment of another collection would give an index into the wrong table
          if (segment.id() != segmentsId)
            throw cms::Exception("Configuration") << "a " << M::kind << " segment of muon " << row << " is in product "
                                                  << segment.id() << ", not in 'segmentSrc' (" << segmentsId << ")";
          muonIdx.push_back(row);
          segmentIdx.push_back(static_cast<int32_t>(segment.key()));
          masks.push_back(match.mask);
        }
      ++row;
    }
    event.put(m_token, table(muonIdx, segmentIdx, masks));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag{"muons"})->setComment("input muon collection (reco::Muon or derived)");
    desc.add<std::string>("cut", "")->setComment("selection of the muons: must be that of the muon table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected muons: must be that of the muon table");
    desc.add<edm::InputTag>("segmentSrc", edm::InputTag{M::segments})
        ->setComment("the segment collection of the matches, whose table segmentIdx indexes");
    desc.add<std::string>("name", std::string("muon") + M::kind + "SegmentMatch")->setComment("name of the table");
    desc.add<std::string>("doc", std::string("muon <-> ") + M::kind + " segment matches")
        ->setComment("documentation of the table");
    desc.add<bool>("skipNonExistingSrc", false)
        ->setComment("an empty table for missing muons or segments instead of failing");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  std::unique_ptr<nanoaod::FlatTable> table(std::vector<int32_t> const& muonIdx,
                                            std::vector<int32_t> const& segmentIdx,
                                            std::vector<uint32_t> const& masks) const {
    auto out = std::make_unique<nanoaod::FlatTable>(muonIdx.size(), m_name, false);
    out->setDoc(m_doc);
    out->addColumn<int32_t>("muonIdx", muonIdx, "row of the muon in the muon table");
    out->addColumn<int32_t>(
        "segmentIdx",
        segmentIdx,
        std::string("index of the segment in its collection (row in the ") + M::kind + " segment table)");
    out->addColumn<uint32_t>("mask", masks, "arbitration bits of the match, reco::MuonSegmentMatch::mask");
    return out;
  }

  const std::string m_name, m_doc;
  const unsigned int m_maxLen;
  const bool m_skipNonExistingSrc;
  const StringCutObjectSelector<reco::Muon> m_cut;
  const edm::EDGetTokenT<edm::View<reco::Muon>> m_muonToken;
  const edm::EDGetTokenT<typename M::Collection> m_segmentToken;
  const edm::EDPutTokenT<nanoaod::FlatTable> m_token;
};

using MuDTSegmentMatchTableProducer = MuSegmentMatchTableProducer<DTMatches>;
using MuCSCSegmentMatchTableProducer = MuSegmentMatchTableProducer<CSCMatches>;
using MuGEMSegmentMatchTableProducer = MuSegmentMatchTableProducer<GEMMatches>;

DEFINE_FWK_MODULE(MuDTSegmentMatchTableProducer);
DEFINE_FWK_MODULE(MuCSCSegmentMatchTableProducer);
DEFINE_FWK_MODULE(MuGEMSegmentMatchTableProducer);
