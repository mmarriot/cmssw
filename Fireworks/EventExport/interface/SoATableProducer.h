#ifndef Fireworks_EventExport_interface_SoATableProducer_h
#define Fireworks_EventExport_interface_SoATableProducer_h

// FlatTable producer for any PortableHostCollection. It walks the layout's ConstDescriptor (leaf
// layouts depth-first, columns in declaration order), so a new or changed SoA needs only a typedef
// and the column names in the configuration:
//
//   blocks = cms.VPSet(cms.PSet(columns = cms.vstring("name[:doc]", ...), table = cms.string("")), ...)
//
// one PSet per leaf layout. An empty column name skips that column. Leaves go to the main table
// unless they name their own; a leaf whose size differs (e.g. the offsets and content of an
// AssociationMap) must name its own table. If the configuration does not match the layout, the
// exception lists the layout's leaves and column types.
//
// Column conversions: arithmetic types as is (int8 -> int16); enums -> underlying value; types with
// rawId() (DetId) -> uint32; other trivially copyable 1/2/4/8-byte types -> raw bits; Eigen columns
// -> name_0 .. name_{R*C-1} in column-major order. Scalars (SOA_SCALAR) are not supported.

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/TypeDemangler.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace evd::soa {

  template <typename T>
  std::string typeName() {
    return edm::typeDemangle(typeid(T).name());
  }

  // Calls f(leafDescriptor) for every leaf layout, depth-first.
  template <typename D, typename F>
  void forEachLeaf(D const& desc, F&& f) {
    if constexpr (requires { D::blocksNumber; })
      std::apply([&](auto const&... block) { (forEachLeaf(block, f), ...); }, desc.buff);
    else
      f(desc);
  }

  // Calls f.template operator()<LeafDescriptor>() for every leaf layout type, depth-first.
  template <typename D, typename F>
  void forEachLeafType(F&& f) {
    if constexpr (requires { D::blocksNumber; }) {
      [&]<std::size_t... I>(std::index_sequence<I...>) {
        (forEachLeafType<std::tuple_element_t<I, decltype(D::buff)>>(f), ...);
      }(std::make_index_sequence<std::tuple_size_v<decltype(D::buff)>>{});
    } else {
      f.template operator()<D>();
    }
  }

  // Element type of column I of a leaf descriptor (the Eigen type for Eigen columns).
  template <typename D, std::size_t I>
  using ColumnValueType = typename std::tuple_element_t<I, decltype(D::parameterTypes)>::ValueType;

  template <std::size_t N>
  using UnsignedOfSize = std::conditional_t<
      N == 1,
      uint8_t,
      std::conditional_t<N == 2, uint16_t, std::conditional_t<N == 4, uint32_t, uint64_t>>>;

  template <typename T>
  void addColumn(nanoaod::FlatTable& table, std::string const& name, std::span<const T> values, std::string const& doc);

  template <typename S, typename T, typename Convert>
  void addConverted(nanoaod::FlatTable& table,
                    std::string const& name,
                    std::span<const T> values,
                    std::string const& doc,
                    Convert convert) {
    std::vector<S> converted(values.size());
    std::transform(values.begin(), values.end(), converted.begin(), convert);
    addColumn<S>(table, name, std::span<const S>(converted), doc);
  }

  template <typename T>
  void addColumn(nanoaod::FlatTable& table, std::string const& name, std::span<const T> values, std::string const& doc) {
    using U = std::remove_cv_t<T>;
    if constexpr (std::is_same_v<U, bool> || std::is_same_v<U, float> || std::is_same_v<U, double> ||
                  std::is_same_v<U, uint8_t> || std::is_same_v<U, int16_t> || std::is_same_v<U, uint16_t> ||
                  std::is_same_v<U, int32_t> || std::is_same_v<U, uint32_t> || std::is_same_v<U, int64_t> ||
                  std::is_same_v<U, uint64_t>)
      table.addColumn<U>(name, values, doc);
    else if constexpr (std::is_integral_v<U>)  // int8_t, char, long long, ...
      table.addColumn<std::conditional_t<std::is_signed_v<U>, std::conditional_t<sizeof(U) < 2, int16_t, int64_t>, uint64_t>>(
          name, values, doc);
    else if constexpr (std::is_enum_v<U>)
      addConverted<std::underlying_type_t<U>>(
          table, name, values, doc, [](U x) { return static_cast<std::underlying_type_t<U>>(x); });
    else if constexpr (requires(U const& x) {
                         { x.rawId() } -> std::convertible_to<uint32_t>;
                       })
      addConverted<uint32_t>(table, name, values, doc, [](U const& x) { return x.rawId(); });
    else if constexpr (std::is_trivially_copyable_v<U> && std::has_unique_object_representations_v<U> &&
                       (sizeof(U) == 1 || sizeof(U) == 2 || sizeof(U) == 4 || sizeof(U) == 8))
      addConverted<UnsignedOfSize<sizeof(U)>>(
          table, name, values, doc + " (raw bits)", [](U const& x) { return std::bit_cast<UnsignedOfSize<sizeof(U)>>(x); });
    else
      throw cms::Exception("Configuration") << "column '" << name << "' has unsupported type " << typeName<U>()
                                            << "; give it an empty name to skip it";
  }

}  // namespace evd::soa

template <typename Collection>
class SoATableProducer : public edm::global::EDProducer<> {
  using Descriptor = typename Collection::ConstDescriptor;

  struct Column {
    std::string name, doc;
  };
  struct Leaf {
    std::size_t table;  // index into tableNames_
    std::vector<Column> columns;
  };

public:
  explicit SoATableProducer(edm::ParameterSet const& cfg)
      : doc_(cfg.getParameter<std::string>("doc")), srcToken_(consumes(cfg.getParameter<edm::InputTag>("src"))) {
    tableNames_.push_back(cfg.getParameter<std::string>("name"));
    for (auto const& block : cfg.getParameter<std::vector<edm::ParameterSet>>("blocks")) {
      Leaf leaf{0, {}};
      if (auto const table = block.getParameter<std::string>("table"); !table.empty()) {
        auto it = std::find(tableNames_.begin(), tableNames_.end(), table);
        leaf.table = it - tableNames_.begin();
        if (it == tableNames_.end())
          tableNames_.push_back(table);
      }
      for (auto const& entry : block.getParameter<std::vector<std::string>>("columns")) {
        auto const colon = entry.find(':');
        leaf.columns.push_back(
            {entry.substr(0, colon), colon == std::string::npos ? std::string() : entry.substr(colon + 1)});
      }
      leaves_.push_back(std::move(leaf));
    }
    checkConfiguration();
    for (std::size_t t = 0; t < tableNames_.size(); ++t)
      putTokens_.push_back(produces<nanoaod::FlatTable>(t == 0 ? std::string() : tableNames_[t]));
  }

  void produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const override {
    auto const view = event.get(srcToken_).const_view();
    std::vector<std::size_t> sizes;
    if constexpr (requires { Descriptor::blocksNumber; }) {
      for (auto s : view.metadata().size())  // block sizes, flattened depth-first like the descriptor's leaves
        sizes.push_back(s);
    } else {
      sizes.push_back(view.metadata().size());
    }

    std::vector<std::unique_ptr<nanoaod::FlatTable>> tables(tableNames_.size());
    std::size_t i = 0;
    evd::soa::forEachLeaf(Descriptor(view), [&](auto const& leaf) {
      std::size_t const rows = leafRows(leaf, sizes[i]);
      auto& table = tables[leaves_[i].table];
      if (!table) {
        table = std::make_unique<nanoaod::FlatTable>(rows, tableNames_[leaves_[i].table], false);
        table->setDoc(doc_);
      } else if (table->size() != rows) {
        throw cms::Exception("Configuration") << "leaf " << i << " has " << rows << " rows, but table '"
                                              << table->name() << "' has " << table->size()
                                              << "; give the leaf its own table";
      }
      fillLeaf(leaf, leaves_[i], rows, *table);
      ++i;
    });
    for (std::size_t t = 0; t < tables.size(); ++t)
      event.put(putTokens_[t], tables[t] ? std::move(tables[t]) : std::make_unique<nanoaod::FlatTable>(0, tableNames_[t], false));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription block;
    block.add<std::vector<std::string>>("columns")->setComment("one 'name[:doc]' per column; empty name skips it");
    block.add<std::string>("table", "")->setComment("output table; empty means the main table");
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src")->setComment("the PortableHostCollection");
    desc.add<std::string>("name")->setComment("name of the main output table");
    desc.add<std::string>("doc", "")->setComment("table documentation");
    desc.addVPSet("blocks", block)->setComment("one PSet per leaf layout, depth-first");
    descriptions.addDefault(desc);
  }

private:
  // Rows of a leaf: the length of its first plain column, which is the leaf's true allocation (e.g. the
  // n+1 offsets of an AssociationMap, whose block size is reported as n). Leaves with only Eigen columns
  // (padded spans) fall back to the layout's block size.
  template <typename D>
  static std::size_t leafRows(D const& leaf, std::size_t blockSize) {
    std::size_t rows = blockSize;
    bool found = false;
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      ((D::columnTypes[I] == cms::soa::SoAColumnType::column && !found
            ? (rows = std::get<I>(leaf.buff).size(), found = true)
            : false),
       ...);
    }(std::make_index_sequence<D::num_cols>{});
    return rows;
  }

  template <typename D>
  static void fillLeaf(D const& leaf, Leaf const& cfg, std::size_t n, nanoaod::FlatTable& table) {
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      (fillColumn<D, I>(leaf, cfg.columns[I], n, table), ...);
    }(std::make_index_sequence<D::num_cols>{});
  }

  template <typename D, std::size_t I>
  static void fillColumn(D const& leaf, Column const& col, std::size_t n, nanoaod::FlatTable& table) {
    if (col.name.empty())
      return;
    auto const& span = std::get<I>(leaf.buff);
    constexpr auto kind = D::columnTypes[I];
    if constexpr (kind == cms::soa::SoAColumnType::column) {
      evd::soa::addColumn(table, col.name, span.first(n), col.doc);
    } else if constexpr (kind == cms::soa::SoAColumnType::eigen) {
      // Each coefficient is stored as its own column, padded to the alignment.
      using E = evd::soa::ColumnValueType<D, I>;
      constexpr std::size_t coefficients = E::RowsAtCompileTime * E::ColsAtCompileTime;
      std::size_t const stride = span.size() / coefficients;
      for (std::size_t k = 0; k < coefficients; ++k)
        evd::soa::addColumn(table, col.name + "_" + std::to_string(k), span.subspan(k * stride, n), col.doc);
    } else {
      throw cms::Exception("Configuration")
          << "scalar column '" << col.name << "' is not supported; give it an empty name to skip it";
    }
  }

  // Compares the configuration with the layout and, on mismatch, throws a description of the layout.
  void checkConfiguration() const {
    std::ostringstream layout;
    std::size_t nLeaves = 0;
    bool ok = true;
    evd::soa::forEachLeafType<Descriptor>([&]<typename D>() {
      layout << "  leaf " << nLeaves << ":";
      [&]<std::size_t... I>(std::index_sequence<I...>) {
        ((layout << ' '
                 << (D::columnTypes[I] == cms::soa::SoAColumnType::scalar   ? "scalar"
                     : D::columnTypes[I] == cms::soa::SoAColumnType::eigen ? "eigen"
                                                                            : "column")
                 << '<' << evd::soa::typeName<evd::soa::ColumnValueType<D, I>>() << '>'),
         ...);
      }(std::make_index_sequence<D::num_cols>{});
      layout << '\n';
      ok &= nLeaves < leaves_.size() && leaves_[nLeaves].columns.size() == D::num_cols;
      ++nLeaves;
    });
    if (!ok || nLeaves != leaves_.size())
      throw cms::Exception("Configuration")
          << "the 'blocks' configuration of table '" << tableNames_[0] << "' does not match the layout of "
          << evd::soa::typeName<typename Collection::Layout>() << " (" << leaves_.size() << " blocks configured, "
          << nLeaves << " in the layout). Expected one PSet per leaf, with these columns:\n"
          << layout.str();
  }

  std::string const doc_;
  edm::EDGetTokenT<Collection> const srcToken_;
  std::vector<std::string> tableNames_;
  std::vector<Leaf> leaves_;
  std::vector<edm::EDPutTokenT<nanoaod::FlatTable>> putTokens_;
};

#endif
