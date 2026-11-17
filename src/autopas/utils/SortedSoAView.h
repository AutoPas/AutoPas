/**
 * @file SortedSoAView.h
 *
 * @date 25.06.2026
 * @author hmeyran
 */

#pragma once

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * A sorted view on a SoA buffer. Particles are projected onto a normalized direction vector, sorted
 * ascending by projection, and packed into a contiguous internal SoA, similar to SortedCellView
 * for AoS.
 *
 * Only the attributes declared by Functor_T::getNeededAttr(std::false_type) are packed from source.
 * Attributes declared by Functor_T::getComputedAttr() are zero-initialized so the functor can
 * accumulate into them. Call scatterBack() after the functor to add the accumulated values back to
 * the source via +=.
 *
 * @note Modifications to the source SoA invalidate the view.
 * @tparam Particle_T Particle type providing AttributeNames and SoAArraysType.
 * @tparam Functor_T Functor type providing getNeededAttr(std::false_type) and getComputedAttr().
 */
template <class Particle_T, class Functor_T>
class SortedSoAView {
  using SoAArraysType = Particle_T::SoAArraysType;

 public:
  /**
   * Projects particles onto sortingDirection, sorts by projection, and packs needed attributes into
   * an internal SoA. Computed attributes are zero-initialized.
   * @param source View of the source SoA buffer to sort.
   * @param sortingDirection Normalized direction vector along which to sort.
   * @param cachedSoa Reference to a persistent SoA buffer to avoid allocations.
   * @param cachedProjIdx Reference to a persistent projection vector to avoid allocations.
   */
  SortedSoAView(SoAView<SoAArraysType> source, const std::array<double, 3> &sortingDirection,
                SoA<SoAArraysType> &cachedSoa, std::vector<std::pair<double, size_t>> &cachedProjIdx)
      : _source(source), _sortedSoa(cachedSoa), _projIdx(cachedProjIdx) {
    const size_t n = source.size();
    if (n == 0) {
      resizeArraysImpl(std::make_index_sequence<Functor_T::getNeededAttr(std::false_type()).size()>{},
                       std::make_index_sequence<Functor_T::getComputedAttr().size()>{}, 0);
      return;
    }

    // Step 1: project each particle onto sortingDirection and sort indices ascending.
    _projIdx.resize(n);
    const auto *xPtr = source.template begin<Particle_T::AttributeNames::posX>();
    const auto *yPtr = source.template begin<Particle_T::AttributeNames::posY>();
    const auto *zPtr = source.template begin<Particle_T::AttributeNames::posZ>();
    for (size_t i = 0; i < n; ++i) {
      _projIdx[i] = {xPtr[i] * sortingDirection[0] + yPtr[i] * sortingDirection[1] + zPtr[i] * sortingDirection[2], i};
    }
    std::sort(_projIdx.begin(), _projIdx.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    // Step 2: resize arrays for needed and computed attributes (allocates extra columns and zeroes computed
    // attributes).
    resizeArraysImpl(std::make_index_sequence<Functor_T::getNeededAttr(std::false_type()).size()>{},
                     std::make_index_sequence<Functor_T::getComputedAttr().size()>{}, n);

    // Step 3: pack getNeededAttr(std::false_type) attributes from source in sorted order.
    // std::false_type as this variant of getNeededAttr includes only the read only attributes, preventing race
    // conditions.
    packNeededImpl(std::make_index_sequence<Functor_T::getNeededAttr(std::false_type()).size()>{}, n);

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
    if (_projIdx.empty()) {
      return;
    }
    scatterBackImpl(std::make_index_sequence<Functor_T::getComputedAttr().size()>{}, _projIdx.size());
  }

 private:
  /**
   * Copies the attribute at compile-time index AttrIdx from _source into _sortedSoa in sorted order.
   * @tparam attribute Index into Functor_T::getNeededAttr(std::false_type).
   * @param n Number of particles.
   */
  template <size_t attribute>
  void packAttr(size_t n) {
    // std::false_type as this variant of getNeededAttr includes only the read only attributes, preventing race
    // conditions.
    constexpr auto attr = static_cast<size_t>(Functor_T::getNeededAttr(std::false_type())[attribute]);
    auto *dst = _sortedSoa.template begin<attr>();
    const auto *src = _source.template begin<attr>();
    for (size_t i = 0; i < n; ++i) {
      dst[i] = src[_projIdx[i].second];
    }
  }

  /**
   * Expands the index sequence and calls packAttr for each needed attribute.
   * @tparam attributes Indices into Functor_T::getNeededAttr(std::false_type).
   * @param n Number of particles.
   */
  template <size_t... attributes>
  void packNeededImpl(std::index_sequence<attributes...>, size_t n) {
    (packAttr<attributes>(n), ...);
  }

  /**
   * Zeroes the computed attribute at compile-time index AttrIdx in _sortedSoa.
   * @tparam attribute Index into Functor_T::getComputedAttr().
   * @param n Number of particles.
   */
  template <size_t attribute>
  void zeroAttr(size_t n) {
    constexpr auto attr = static_cast<size_t>(Functor_T::getComputedAttr()[attribute]);
    auto *ptr = _sortedSoa.template begin<attr>();
    std::fill(ptr, ptr + n, 0);
  }

  /**
   * Expands the index sequence and calls zeroAttr for each computed attribute.
   * @tparam attributes Indices into Functor_T::getComputedAttr().
   * @param n Number of particles.
   */
  template <size_t... attributes>
  void zeroComputedImpl(std::index_sequence<attributes...>, size_t n) {
    (zeroAttr<attributes>(n), ...);
  }

  /**
   * Resizes _sortedSoa for all needed and computed attributes of Functor_T.
   * @tparam neededAttributes Indices into Functor_T::getNeededAttr(std::false_type()).
   * @tparam computedAttributes Indices into Functor_T::getComputedAttr().
   * @param n Number of particles.
   */
  template <size_t... neededAttributes, size_t... computedAttributes>
  void resizeArraysImpl(std::index_sequence<neededAttributes...>, std::index_sequence<computedAttributes...>,
                        size_t n) {
    _sortedSoa.template resizeArraysForAttributes<
        static_cast<size_t>(Functor_T::getNeededAttr(std::false_type())[neededAttributes])...,
        static_cast<size_t>(Functor_T::getComputedAttr()[computedAttributes])...>(n);
  }

  /**
   * Adds the computed attribute at compile-time index AttrIdx from _sortedSoa back to _source
   * using the original-index mapping in projIdx.
   * @tparam attribute Index into Functor_T::getComputedAttr().
   * @param n Number of particles.
   */
  template <size_t attribute>
  void scatterAttr(size_t n) {
    constexpr auto attr = static_cast<size_t>(Functor_T::getComputedAttr()[attribute]);
    auto *orig = _source.template begin<attr>();
    const auto *sorted = _sortedSoa.template begin<attr>();
    for (size_t i = 0; i < n; ++i) {
      orig[_projIdx[i].second] += sorted[i];
    }
  }

  /**
   * Expands the index sequence and calls scatterAttr for each computed attribute.
   * @tparam attributes Indices into Functor_T::getComputedAttr().
   * @param n Number of particles.
   */
  template <size_t... attributes>
  void scatterBackImpl(std::index_sequence<attributes...>, size_t n) {
    (scatterAttr<attributes>(n), ...);
  }

  SoAView<SoAArraysType> _source;
  SoA<SoAArraysType> &_sortedSoa;
  /**
   * Sorted (projection value, original index) pairs, ascending by projection.
   */
  std::vector<std::pair<double, size_t>> &_projIdx;
};

}  // namespace autopas
