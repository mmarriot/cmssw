// The values of jet tags (reco::JetTagCollection: b-tag and other jet tagger discriminators) as columns of an
// extension of a jet table, for jets that are not PAT jets (which carry them as bDiscriminator). Each column is the
// sum of one or more JetTag collections, such as probb + probbb, written as the taggers give them (including their
// own default values, such as -1 for a jet without tracks). Tagger inputs are written by JetTaggerTableProducer.
// A JetTagCollection holds one value per jet of the collection it is keyed on, in its order; the producer checks that
// this is the jet collection of the table whenever there are jets. Without jets the tags are not read: a tagger with no
// jets to tag writes an empty JetTagCollection, keyed on no collection.

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterWildcard.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/transform.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

class JetTagTableProducer final : public edm::global::EDProducer<> {
public:
  explicit JetTagTableProducer(edm::ParameterSet const& cfg)
      : name_(cfg.getParameter<std::string>("name")),
        skipNonExistingSrc_(cfg.getParameter<bool>("skipNonExistingSrc")),
        maxLen_(cfg.existsAs<unsigned int>("maxLen") ? cfg.getParameter<unsigned int>("maxLen")
                                                     : std::numeric_limits<unsigned int>::max()),
        cut_(cfg.getParameter<std::string>("cut")),
        jetsToken_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        tableToken_(produces()) {
    auto const& variables = cfg.getParameter<edm::ParameterSet>("variables");
    for (auto const& name : variables.getParameterNamesForType<edm::ParameterSet>()) {
      auto const& variable = variables.getParameter<edm::ParameterSet>(name);
      columns_.push_back(
          {name,
           variable.getParameter<std::string>("doc"),
           variable.getParameter<int>("precision"),
           edm::vector_transform(variable.getParameter<std::vector<edm::InputTag>>("src"),
                                 [this](edm::InputTag const& tag) { return consumes<reco::JetTagCollection>(tag); })});
    }
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const jets = event.getHandle(jetsToken_);
    // no jets: empty columns (jets->empty() throws for missing jets, unless skipNonExistingSrc)
    const bool noJets = (!jets.isValid() && skipNonExistingSrc_) || jets->empty();
    std::vector<std::size_t> selected;  // the rows of the jet table, as SimpleFlatTableProducer selects them
    for (std::size_t i = 0; !noJets && i < jets->size() && selected.size() < maxLen_; ++i)
      if (cut_((*jets)[i]))
        selected.push_back(i);

    auto table = std::make_unique<nanoaod::FlatTable>(selected.size(), name_, false, true);
    std::vector<float> sum(noJets ? 0 : jets->size()), values(selected.size());
    for (auto const& column : columns_) {
      std::fill(sum.begin(), sum.end(), 0.f);
      for (auto const& token : column.tagTokens) {
        if (noJets)
          break;
        auto const& tags = event.get(token);
        if (tags.empty() || tags.keyProduct().id() != jets.id() || tags.size() != jets->size())
          throw cms::Exception("Configuration")
              << "the " << tags.size() << " jet tags of column " << column.name << " are not those of the "
              << jets->size() << " jets of src " << jets.id();
        for (std::size_t i = 0; i < tags.size(); ++i)
          sum[i] += tags.value(i);  // value(i) is the tag of jet i, without building the (jet, value) pairs
      }
      for (std::size_t row = 0; row < selected.size(); ++row)
        values[row] = sum[selected[row]];
      table->addColumn<float>(column.name, values, column.doc, column.precision);
    }
    event.put(tableToken_, std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("the jets, the collection the jet tags are keyed on");
    desc.add<std::string>("cut", "")->setComment(
        "selection of the jets, on reco::Jet methods: must be that of the jet table");
    desc.addOptional<unsigned int>("maxLen")->setComment(
        "maximum number of selected jets: must be that of the jet table");
    desc.add<std::string>("name")->setComment("name of the jet table, extended with one column per variable");
    desc.add<bool>("skipNonExistingSrc", false)->setComment("an empty table for missing jets instead of failing");
    edm::ParameterSetDescription variable;
    variable.add<std::vector<edm::InputTag>>("src")->setComment("JetTagCollections whose values are summed");
    variable.add<std::string>("doc", "")->setComment("documentation of the column");
    variable.add<int>("precision", -1)->setComment("mantissa bits kept (zero or negative: all)");
    edm::ParameterSetDescription variables;
    variables.addNode(
        edm::ParameterWildcard<edm::ParameterSetDescription>("*", edm::RequireZeroOrMore, true, variable));
    desc.add<edm::ParameterSetDescription>("variables", variables)->setComment("one column per parameter set");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  struct Column {
    std::string name, doc;
    int precision;
    std::vector<edm::EDGetTokenT<reco::JetTagCollection>> tagTokens;
  };

  const std::string name_;
  const bool skipNonExistingSrc_;
  const unsigned int maxLen_;
  const StringCutObjectSelector<reco::Jet> cut_;
  const edm::EDGetTokenT<edm::View<reco::Jet>> jetsToken_;
  const edm::EDPutTokenT<nanoaod::FlatTable> tableToken_;
  std::vector<Column> columns_;
};

DEFINE_FWK_MODULE(JetTagTableProducer);
