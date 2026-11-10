/**
 * @file SoAView.h
 * @authors humig
 * @date 07.07.2019
 */

#pragma once

#include <algorithm>
#include <utility>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <class SoAArraysType>
class SoA;

/**
 * View on a fixed part of a SoA between a start index and an end index.
 *
 * It is the user`s responsibility to ensure that modifications to the underlying SoA don't let the start index get
 * bigger than the highest index in the SoA, or the end index higher than the number of entries in the SoA.
 *
 * @tparam SoAArraysType The SoAArrayType of the SoA.
 */
template <class SoAArraysType>
class SoAView {
 public:
  /**
   * Default constructor of SoAView to allow storing it in containers.
   */
  SoAView() : _soa(nullptr), _startIndex(0), _endIndex(0) {}

  /**
   * Constructs a view on \p soa that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `soa->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param soa The SoA to view.
   * @param startIndex The index of the first entry to view of the SoA.
   * @param endIndex The index of the entry after the last entry to view of the SoA.
   */
  SoAView(SoA<SoAArraysType> *soa, size_t startIndex, size_t endIndex)
      : _soa(soa), _startIndex(startIndex), _endIndex(endIndex) {
    if (not(soa->size() >= endIndex and endIndex >= startIndex)) /* @todo C++20 [[unlikely]] */ {
      utils::ExceptionHandler::exception("SoAView: Trying to view particles outside of the SoA.");
    }
  }

  /**
   * Constructs a SoAView on the whole content of \p soa.
   * @param soa The SoA to view.
   */
  explicit SoAView(SoA<SoAArraysType> *soa) : _soa(soa), _startIndex(0), _endIndex(soa->size()) {}

  /**
   * Implicit constructor that converts a SoA to SoAView.
   * @param soa The SoA to view.
   */
  SoAView(SoA<SoAArraysType> &soa) : _soa(&soa), _startIndex(0), _endIndex(soa.size()) {}

  /**
   * Returns a pointer to the given attribute vector.
   * If the underlying SoA column has no elements for this view range, returns nullptr.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector within this view, or nullptr.
   */
  template <size_t attribute>
  [[nodiscard]] auto begin() {
    using ReturnType = decltype(_soa->template begin<attribute>());
    if (_soa == nullptr or _soa->template size<attribute>() <= _startIndex) {
      return static_cast<ReturnType>(nullptr);
    }
    return _soa->template begin<attribute>() + _startIndex;
  }

  /**
   * Returns a pointer to the given attribute vector const.
   * If the underlying SoA column has no elements for this view range, returns nullptr.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector const within this view, or nullptr.
   */
  template <size_t attribute>
  [[nodiscard]] auto begin() const {
    using ReturnType = decltype(std::as_const(*_soa).template begin<attribute>());
    if (_soa == nullptr or _soa->template size<attribute>() <= _startIndex) {
      return static_cast<ReturnType>(nullptr);
    }
    return _soa->template begin<attribute>() + _startIndex;
  }

  /**
   * Returns the number of valid elements for the given attribute in this view.
   * If the underlying SoA column has fewer elements than the view range (e.g. because
   * it was not loaded/resized for this functor), this returns the actual available count.
   * @tparam attribute ID of the desired attribute.
   * @return Available number of elements for this attribute.
   */
  template <size_t attribute>
  [[nodiscard]] size_t size() const {
    if (_soa == nullptr) {
      return 0;
    }
    const size_t underlyingSize = _soa->template size<attribute>();
    if (underlyingSize <= _startIndex) {
      return 0;
    }
    return std::min(underlyingSize - _startIndex, size());
  }

  /**
   * Returns the number of particles in the view.
   *
   * @return Number of particles.
   */
  [[nodiscard]] size_t size() const { return _endIndex - _startIndex; }

  /**
   * Returns the number of extra columns in the underlying SoA.
   * @return Number of extra columns.
   */
  [[nodiscard]] size_t numExtraColumns() const { return _soa ? _soa->numExtraColumns() : 0; }

  /**
   * Returns the available size of the given extra column in this view.
   * @param col Index of the extra column (0-based).
   * @return Number of available elements for this extra column in the view.
   */
  [[nodiscard]] size_t extraColumnSize(size_t col) const {
    if (_soa == nullptr) {
      return 0;
    }
    const size_t total = _soa->extraColumnSize(col);
    if (total <= _startIndex) {
      return 0;
    }
    return std::min(total - _startIndex, size());
  }

  /**
   * Returns a const pointer to the beginning of the given extra column in this view.
   * @param col Index of the extra column (0-based).
   * @return Pointer to data within the view, or nullptr if unavailable.
   */
  [[nodiscard]] const double *extraColumnBegin(size_t col) const {
    if (_soa == nullptr or _soa->extraColumnSize(col) <= _startIndex) {
      return nullptr;
    }
    return _soa->extraColumnBegin(col) + _startIndex;
  }

  /**
   * Returns a pointer to the beginning of the given extra column in this view.
   * @param col Index of the extra column (0-based).
   * @return Pointer to data within the view, or nullptr if unavailable.
   */
  [[nodiscard]] double *extraColumnBegin(size_t col) {
    if (_soa == nullptr or _soa->extraColumnSize(col) <= _startIndex) {
      return nullptr;
    }
    return _soa->extraColumnBegin(col) + _startIndex;
  }

 private:
  /**
   * The underlying SoA.
   */
  SoA<SoAArraysType> *_soa;

  /**
   * The start index of the view in the SoA. (Inclusive)
   */
  size_t _startIndex;

  /**
   * The end index of the view in the SoA. (Exclusive)
   */
  size_t _endIndex;
};

}  // namespace autopas
