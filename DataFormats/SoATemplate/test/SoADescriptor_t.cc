#include <cstdlib>
#include <memory>
#include <string_view>
#include <tuple>

#include <Eigen/Core>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

// Test the metadata exposed by the (const) descriptors: column names, block names and number of elements

// The element methods sit between the columns, to check that they do not shift the column names
GENERATE_SOA_LAYOUT(FlatLayout,
                    SOA_COLUMN(int, energy),
                    SOA_CONST_ELEMENT_METHODS(SOA_HOST_DEVICE int twice() const { return 2 * energy(); }),
                    SOA_EIGEN_COLUMN(Eigen::Vector3d, position),
                    SOA_SCALAR(int, id))

// A layout with only an Eigen column, that allocates one more element than requested, like the offsets of
// an AssociationMap: its number of elements cannot be recovered from the (padded) span of its column
GENERATE_SOA_LAYOUT(EigenOnlyLayout, SOA_EIGEN_COLUMN(Eigen::Vector3f, axis))

template <CMS_SOA_BYTE_SIZE_TYPE ALIGNMENT = cms::soa::CacheLineSize::defaultSize,
          bool ALIGNMENT_ENFORCEMENT = cms::soa::AlignmentEnforcement::relaxed>
class ExtendedLayout final : public EigenOnlyLayout<ALIGNMENT, ALIGNMENT_ENFORCEMENT> {
  using Parent = EigenOnlyLayout<ALIGNMENT, ALIGNMENT_ENFORCEMENT>;

public:
  ExtendedLayout() : Parent() {}
  ExtendedLayout(std::byte* mem, cms::soa::size_type elements) : Parent(mem, elements + 1) {}

  static constexpr auto computeDataSize(std::size_t size) { return Parent::computeDataSize(size + 1); }
};

GENERATE_SOA_BLOCKS(InnerBlocks, SOA_BLOCK(flat, FlatLayout), SOA_BLOCK(extended, ExtendedLayout))

GENERATE_SOA_BLOCKS(OuterBlocks, SOA_BLOCK(inner, InnerBlocks), SOA_BLOCK(other, FlatLayout))

using namespace std::literals;

TEST_CASE("SoADescriptor column names") {
  using Flat = FlatLayout<>;
  STATIC_REQUIRE(Flat::ConstDescriptor::num_cols == 3);
  STATIC_REQUIRE(Flat::ConstDescriptor::columnNames == std::array{"energy"sv, "position"sv, "id"sv});
  STATIC_REQUIRE(Flat::Descriptor::columnNames == Flat::ConstDescriptor::columnNames);
}

TEST_CASE("SoADescriptor block names") {
  using Outer = OuterBlocks<>;
  // the layout counts the leaf blocks, the descriptor only its direct blocks
  STATIC_REQUIRE(Outer::blocksNumber == 3);
  STATIC_REQUIRE(Outer::ConstDescriptor::blocksNumber == 2);
  STATIC_REQUIRE(Outer::ConstDescriptor::blockNames == std::array{"inner"sv, "other"sv});
  STATIC_REQUIRE(Outer::Descriptor::blockNames == Outer::ConstDescriptor::blockNames);
  STATIC_REQUIRE(InnerBlocks<>::ConstDescriptor::blockNames == std::array{"flat"sv, "extended"sv});
}

TEST_CASE("SoADescriptor elements") {
  using Outer = OuterBlocks<>;
  const std::array<Outer::size_type, 3> sizes = {10, 20, 30};

  std::unique_ptr<std::byte, decltype(std::free)*> buffer{
      reinterpret_cast<std::byte*>(aligned_alloc(Outer::alignment, Outer::computeDataSize(sizes))), std::free};
  Outer layout{buffer.get(), sizes};
  Outer::View view{layout};
  Outer::ConstView constView{layout};

  SECTION("Descriptor") {
    Outer::Descriptor desc{view};
    auto const& flat = std::get<0>(std::get<0>(desc.buff).buff);
    auto const& extended = std::get<1>(std::get<0>(desc.buff).buff);
    auto const& other = std::get<1>(desc.buff);

    REQUIRE(flat.elements == 10);
    REQUIRE(extended.elements == 21);
    REQUIRE(other.elements == 30);
  }

  SECTION("ConstDescriptor") {
    Outer::ConstDescriptor desc{constView};
    auto const& flat = std::get<0>(std::get<0>(desc.buff).buff);
    auto const& extended = std::get<1>(std::get<0>(desc.buff).buff);
    auto const& other = std::get<1>(desc.buff);

    // the parent reports the requested size, the leaf its own number of elements
    REQUIRE(constView.metadata().size()[1] == 20);
    REQUIRE(extended.elements == 21);
    REQUIRE(flat.elements == 10);
    REQUIRE(other.elements == 30);

    // the span of a plain column has one entry per element, the span of an Eigen column includes the padding
    REQUIRE(std::get<0>(flat.buff).size() == 10);
    REQUIRE(std::get<0>(extended.buff).size() != 21 * 3);
  }
}
