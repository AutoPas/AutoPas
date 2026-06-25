/**
 * @file SoASortedView.h
 *
 * @date 25.06.2026
 * @author Henry Meyran
 */

#pragma once

#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "SoA.h"
#include "SoAView.h"

namespace autopas {

/**
 * Precomputed index bounds for iterating a pre-sorted SoA pair. Produced by CellFunctor and
 * consumed by SoAFunctorPairSorted overrides.
 */
struct SoASortingData {
  size_t start_i;
  std::vector<size_t> maxIndex;
  std::vector<size_t> minIndex;
};

/**
 * A sorted view on a SoA buffer. Particles are projected onto a normalized direction vector, sorted
 * ascending by projection, and packed into a contiguous internal SoA, similar to SortedCellView
 * for AoS.
 *
 * Only the attributes declared by Functor_T::getNeededAttr() are packed from source.
 * Attributes declared by Functor_T::getComputedAttr() are zero-initialized so the functor can
 * accumulate into them. Call scatterBack() after the functor to add the accumulated values back to
 * the source via +=.
 *
 * @note Modifications to the source SoA invalidate the view.
 * @tparam Particle_T Particle type providing AttributeNames and SoAArraysType.
 * @tparam Functor_T Functor type providing getNeededAttr() and getComputedAttr().
 */
template <class Particle_T, class Functor_T>
class SoASortedView {
  using SoAArraysType = typename Particle_T::SoAArraysType;

 public:
  /**
   * Projects particles onto sortingDirection, sorts by projection, and packs needed attributes into
   * an internal SoA. Computed attributes are zero-initialized.
   * @param source View of the source SoA buffer to sort.
   * @param sortingDirection Normalized direction vector along which to sort.
   */
  SoASortedView(SoAView<SoAArraysType> source, const std::array<double, 3> &sortingDirection) : _source(source) {
    const size_t n = source.size();
    if (n == 0) return;

    // Step 1: project each particle onto sortingDirection and sort indices ascending.
    projIdx.resize(n);
    const auto *xPtr = source.template begin<Particle_T::AttributeNames::posX>();
    const auto *yPtr = source.template begin<Particle_T::AttributeNames::posY>();
    const auto *zPtr = source.template begin<Particle_T::AttributeNames::posZ>();
    for (size_t i = 0; i < n; ++i) {
      projIdx[i] = {xPtr[i] * sortingDirection[0] + yPtr[i] * sortingDirection[1] + zPtr[i] * sortingDirection[2], i};
    }
    std::sort(projIdx.begin(), projIdx.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    // Step 2: resize all arrays (value-initializes to 0, covering the computed attributes).
    _sortedSoa.resizeArrays(n);

    // Step 3: pack getNeededAttr() attributes from source in sorted order.
    packNeededImpl(std::make_index_sequence<Functor_T::getNeededAttr().size()>{}, n);

    // Step 4: zero getComputedAttr() attributes so the functor accumulates from 0.
    zeroComputedImpl(std::make_index_sequence<Functor_T::getComputedAttr().size()>{}, n);
  }

  /**
   * Returns a SoAView into the sorted, contiguous internal buffer.
   * @return SoAView of the sorted SoA.
   */
  SoAView<SoAArraysType> getView() { return _sortedSoa.constructView(); }

  /**
   * Scatters computed attribute contributions accumulated in the sorted SoA back to the source
   * via +=. Must be called after the functor to commit writes.
   */
  void scatterBack() {
    scatterBackImpl(std::make_index_sequence<Functor_T::getComputedAttr().size()>{}, projIdx.size());
  }

  /**
   * Sorted (projection value, original index) pairs, ascending by projection.
   */
  std::vector<std::pair<double, size_t>> projIdx;

 private:
  template <size_t AttrIdx>
  void packAttr(size_t n) {
    constexpr size_t attr = static_cast<size_t>(Functor_T::getNeededAttr()[AttrIdx]);
    auto *dst = _sortedSoa.template begin<attr>();
    const auto *src = _source.template begin<attr>();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = src[projIdx[i].second];
    }
  }

  template <size_t... AttrIdxs>
  void packNeededImpl(std::index_sequence<AttrIdxs...>, size_t n) {
    (packAttr<AttrIdxs>(n), ...);
  }

  template <size_t AttrIdx>
  void zeroAttr(size_t n) {
    constexpr size_t attr = static_cast<size_t>(Functor_T::getComputedAttr()[AttrIdx]);
    auto *ptr = _sortedSoa.template begin<attr>();
    std::fill(ptr, ptr + n, 0);
  }

  template <size_t... AttrIdxs>
  void zeroComputedImpl(std::index_sequence<AttrIdxs...>, size_t n) {
    (zeroAttr<AttrIdxs>(n), ...);
  }

  template <size_t AttrIdx>
  void scatterAttr(size_t n) {
    constexpr size_t attr = static_cast<size_t>(Functor_T::getComputedAttr()[AttrIdx]);
    auto *orig = _source.template begin<attr>();
    const auto *sorted = _sortedSoa.template begin<attr>();
    for (size_t i = 0; i < n; ++i) {
      orig[projIdx[i].second] += sorted[i];
    }
  }

  template <size_t... AttrIdxs>
  void scatterBackImpl(std::index_sequence<AttrIdxs...>, size_t n) {
    (scatterAttr<AttrIdxs>(n), ...);
  }

  SoAView<SoAArraysType> _source;
  SoA<SoAArraysType> _sortedSoa;
};

}  // namespace autopas
