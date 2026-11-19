#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <map>
#include <type_traits>

namespace {
  // The row of each candidate of a view, by the product and key of its Ptr. A view of a collection holds the keys
  // 0..n-1 of one product, and a view of Ptrs (a merged selection of packed candidates) keys of a few products, so
  // one key -> row vector per product turns the search for a constituent into two lookups (a linear search costs
  // tens of ms per event for the thousands of constituents of RECO jets among the PF candidates at 200 pileup).
  class CandidateRows {
  public:
    explicit CandidateRows(std::vector<reco::CandidatePtr> const &candidates) {
      for (int row = 0; auto const &ptr : candidates) {
        if (ptr.isNonnull()) {
          auto &rows = rows_[ptr.id()];
          if (ptr.key() >= rows.size())
            rows.resize(ptr.key() + 1, -1);
          if (rows[ptr.key()] < 0)  // the first row holding it, as a linear search would find
            rows[ptr.key()] = row;
        }
        ++row;
      }
    }

    // -1 if the candidate is not in the view
    int row(reco::CandidatePtr const &ptr) const {
      auto const found = rows_.find(ptr.id());
      return (found == rows_.end() || ptr.key() >= found->second.size()) ? -1 : found->second[ptr.key()];
    }

  private:
    std::map<edm::ProductID, std::vector<int>> rows_;
  };
}  // namespace

// The constituents of jets (their daughters) or of taus (their signal candidates) as a table of rows of a candidate
// table: for each constituent, the index of its candidate in `candidates` and of its jet or tau among the selected ones.
template <typename T, typename C = std::vector<typename T::ConstituentTypeFwdPtr>>
class SimpleJetConstituentTableProducer : public edm::stream::EDProducer<> {
public:
  using ConstituentsOutput = C;
  using ConstituentValueType = typename C::value_type;

  explicit SimpleJetConstituentTableProducer(const edm::ParameterSet &);
  ~SimpleJetConstituentTableProducer() override;

  ConstituentValueType const initptr(edm::Ptr<reco::Candidate> const &) const;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  static std::vector<reco::CandidatePtr> const &constituents(T const &object) {
    if constexpr (std::is_same_v<T, reco::PFTau>)
      return object.signalCands();
    else
      return object.daughterPtrVector();
  }

  void produce(edm::Event &, const edm::EventSetup &) override;

  const std::string name_;
  const std::string candIdxName_;
  const std::string candIdxDoc_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;

  const StringCutObjectSelector<T> jetCut_;
  const StringCutObjectSelector<ConstituentValueType> jetConstCut_;
};

//
// constructors and destructor
//
template <typename T, typename C>
SimpleJetConstituentTableProducer<T, C>::SimpleJetConstituentTableProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      candIdxName_(iConfig.getParameter<std::string>("candIdxName")),
      candIdxDoc_(iConfig.getParameter<std::string>("candIdxDoc")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      jetCut_(iConfig.getParameter<std::string>("jetCut")),
      jetConstCut_(iConfig.getParameter<std::string>("jetConstCut")) {
  produces<nanoaod::FlatTable>(name_);
  produces<ConstituentsOutput>();
}

template <typename T, typename C>
SimpleJetConstituentTableProducer<T, C>::~SimpleJetConstituentTableProducer() {}

template <typename T, typename C>
void SimpleJetConstituentTableProducer<T, C>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  auto outCands = std::make_unique<ConstituentsOutput>();

  auto jets = iEvent.getHandle(jet_token_);

  edm::Handle<reco::CandidateView> cands_;
  iEvent.getByToken(cand_token_, cands_);
  CandidateRows const candRows(cands_->ptrs());

  // First, select jets
  std::vector<T> jetsPassCut;
  for (unsigned jetIdx = 0; jetIdx < jets->size(); ++jetIdx) {
    const auto &jet = jets->at(jetIdx);
    if (!jetCut_(jet))
      continue;
    jetsPassCut.push_back(jets->at(jetIdx));
  }

  // Then loop over selected jets
  std::vector<int> parentJetIdx;
  std::vector<int> candIdx;
  for (unsigned jetIdx = 0; jetIdx < jetsPassCut.size(); ++jetIdx) {
    const auto &jet = jetsPassCut.at(jetIdx);

    // Loop over jet constituents
    std::vector<reco::CandidatePtr> const &daughters = constituents(jet);
    for (const auto &dauPtr : daughters) {
      // Apply cut on jet constituent
      typename C::value_type cand = initptr(dauPtr);
      if (!jetConstCut_(cand))
        continue;

      // Find jet constituent in candidate collection
      const int row = candRows.row(dauPtr);
      if (row < 0) {
        continue;
      }

      outCands->push_back(cand);
      parentJetIdx.push_back(jetIdx);
      candIdx.push_back(row);
    }
  }  // end jet loop

  auto candTable = std::make_unique<nanoaod::FlatTable>(outCands->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  candTable->template addColumn<int>(candIdxName_, candIdx, candIdxDoc_);

  std::string parentJetIdxName("jetIdx");
  std::string parentJetIdxDoc("Index of the parent jet");
  if constexpr (std::is_same<T, reco::GenJet>::value) {
    parentJetIdxName = "genJetIdx";
    parentJetIdxDoc = "Index of the parent gen jet";
  } else if constexpr (std::is_same<T, reco::PFTau>::value) {
    parentJetIdxName = "tauIdx";
    parentJetIdxDoc = "Index of the parent tau";
  }
  candTable->template addColumn<int>(parentJetIdxName, parentJetIdx, parentJetIdxDoc);

  iEvent.put(std::move(candTable), name_);
  iEvent.put(std::move(outCands));
}

template <typename T, typename C>
typename SimpleJetConstituentTableProducer<T, C>::ConstituentValueType const
SimpleJetConstituentTableProducer<T, C>::initptr(edm::Ptr<reco::Candidate> const &dau) const {
  return ConstituentValueType(dau);
}

template <typename T, typename C>
void SimpleJetConstituentTableProducer<T, C>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "FatJetPFCand");
  desc.add<std::string>("candIdxName", "PFCandIdx");
  desc.add<std::string>("candIdxDoc", "Index in PFCand table");
  desc.add<edm::InputTag>("jets", edm::InputTag("finalJetsAK8"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<std::string>("jetCut", "");
  desc.add<std::string>("jetConstCut", "");
  descriptions.addWithDefaultLabel(desc);
}

typedef SimpleJetConstituentTableProducer<pat::Jet, std::vector<edm::Ptr<pat::PackedCandidate>>>
    SimplePatJetConstituentTableProducer;
typedef SimpleJetConstituentTableProducer<reco::GenJet, std::vector<edm::Ptr<pat::PackedGenParticle>>>
    SimpleGenJetConstituentTableProducer;

// RECO jets and taus: the constituents are Ptrs to the PF candidates (particleFlow), read through reco::Candidate
typedef SimpleJetConstituentTableProducer<reco::PFJet, std::vector<reco::CandidatePtr>>
    SimplePFJetConstituentTableProducer;
typedef SimpleJetConstituentTableProducer<reco::PFTau, std::vector<reco::CandidatePtr>>
    SimplePFTauConstituentTableProducer;

DEFINE_FWK_MODULE(SimplePatJetConstituentTableProducer);
DEFINE_FWK_MODULE(SimpleGenJetConstituentTableProducer);
DEFINE_FWK_MODULE(SimplePFJetConstituentTableProducer);
DEFINE_FWK_MODULE(SimplePFTauConstituentTableProducer);
