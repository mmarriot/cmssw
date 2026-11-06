#ifndef DataFormats_Portable_interface_SoAConstMultiView_h
#define DataFormats_Portable_interface_SoAConstMultiView_h

#include <array>
#include <cassert>
#include <concepts>
#include <functional>
#include <ranges>
#include <span>

#include "DataFormats/SoATemplate/interface/SoACommon.h"

/**
 * @brief Aggregates multiple ConstViews into a single combined view.
 *
 * `SoAConstMultiView` stores multiple ConstViews within an `std::array`, 
 * accompanied by an offset array to enable access via a global index. 
 * An `SoAConstMultiView` can be passed directly by value to kernels without 
 * requiring device memory copies. As the SoAConstMultiView accepts only ConstViews, 
 * it supports read-only access to the underlying data.
 *
 * Since the underlying SoA memories are not contiguous, cacheline inefficiencies 
 * may arise. Therefore, `SoAConstMultiView` is best suited for use with few ConstViews, 
 * when the underlying buffer is large enough to amortize the overhead of 
 * non-contiguous access patterns when iterating over all elements.
 *
 */

template <typename Collections, typename Getter, typename ConstView>
concept IsValidAViewGetter = requires(Getter getter, std::ranges::range_reference_t<Collections> collection) {
  { std::invoke(getter, collection) } -> std::convertible_to<ConstView>;
};

template <typename ConstView, cms::soa::size_type MaxSize>
class SoAConstMultiView {
public:
  using ConstElement = typename ConstView::const_element;
  using size_type = cms::soa::size_type;
  constexpr static cms::soa::RangeChecking::Mode rangeChecking = ConstView::rangeChecking;

  SoAConstMultiView() = default;

  template <std::ranges::sized_range Collections>
  explicit SoAConstMultiView(Collections const& collections) {
    for (const auto& collection : collections) {
      ConstView view = collection.const_view();
      addView(view, view.metadata().size());
    }
  }

  template <std::ranges::sized_range Collections>
  explicit SoAConstMultiView(Collections const& collections, std::span<const size_type> sizes) {
    assert(std::ranges::size(collections) == sizes.size() &&
           "The number of sizes must match the number of collections");
    for (const auto& collection : collections) {
      ConstView view = collection.const_view();
      assert(n_ < static_cast<size_type>(sizes.size()) && "More collections provided than sizes");
      addView(view, sizes[n_]);
    }
    assert(n_ == static_cast<size_type>(sizes.size()));
  }

  template <std::ranges::sized_range Collections, typename Getter>
    requires IsValidAViewGetter<Collections, Getter, ConstView>
  explicit SoAConstMultiView(Collections const& collections, Getter getter) {
    for (const auto& collection : collections) {
      ConstView view = std::invoke(getter, collection);
      addView(view, view.metadata().size());
    }
  }

  template <std::ranges::sized_range Collections, typename Getter>
    requires IsValidAViewGetter<Collections, Getter, ConstView>
  explicit SoAConstMultiView(Collections const& collections, Getter getter, std::span<const size_type> sizes) {
    assert(std::ranges::size(collections) == sizes.size() &&
           "The number of sizes must match the number of collections");
    for (const auto& collection : collections) {
      assert(n_ < static_cast<size_type>(sizes.size()) && "More collections provided than sizes");
      addView(std::invoke(getter, collection), sizes[n_]);
    }
    assert(n_ == static_cast<size_type>(sizes.size()));
  }

  void addView(ConstView view, const size_type size) {
    assert(n_ < MaxSize && "Exceeded maximum number of views");
    assert(size <= static_cast<size_type>(view.metadata().size()) && "Provided size exceeds elements in the view");
    views_[n_] = view;
    // teh unused views end where the last one does
    const size_type end = (n_ == 0 ? 0 : offsets_[n_ - 1]) + size;
    for (size_type i = n_; i < MaxSize; ++i) {
      offsets_[i] = end;
    }
    n_++;
  }

  [[nodiscard]] SOA_HOST_DEVICE SOA_INLINE ConstElement
  operator[](cms::soa::detail::IndexWithSourceLocation<rangeChecking> globalIndex) const {
    if constexpr (rangeChecking != cms::soa::RangeChecking::disabled) {
      if (globalIndex.value_ >= size() or globalIndex.value_ < 0) {
        SOA_THROW_OUT_OF_RANGE("Index surpasses size total size of SoAConstMultiView", globalIndex, size())
      }
    }

    return find<0>(globalIndex.value_, 0);
  }

  [[nodiscard]] SOA_HOST_DEVICE SOA_INLINE ConstView
  view(cms::soa::detail::IndexWithSourceLocation<rangeChecking> i) const {
    if constexpr (rangeChecking != cms::soa::RangeChecking::disabled) {
      if (i.value_ >= n_) {
        SOA_THROW_OUT_OF_RANGE("Out of range index in SoAConstMultiView::view()", i, n_)
      }
    }
    return views_[i.value_];
  }

  [[nodiscard]] SOA_HOST_DEVICE SOA_INLINE size_type size() const { return offsets_[MaxSize - 1]; }
  [[nodiscard]] SOA_HOST_DEVICE SOA_INLINE size_type numViews() const { return n_; }

private:
  // Compile-time indices: a runtime index into views_ makes every thread copy the kernel arguments to local memory
  template <size_type K>
  SOA_HOST_DEVICE SOA_INLINE ConstElement find(const size_type globalIndex, const size_type first) const {
    if constexpr (K + 1 == MaxSize) {
      return views_[K][globalIndex - first];
    } else {
      if (globalIndex < offsets_[K]) {
        return views_[K][globalIndex - first];
      }
      return find<K + 1>(globalIndex, offsets_[K]);
    }
  }

  std::array<ConstView, MaxSize> views_{};
  std::array<size_type, MaxSize> offsets_{};

  size_type n_{0};
};

#endif  // DataFormats_Portable_interface_SoAConstMultiView_h
