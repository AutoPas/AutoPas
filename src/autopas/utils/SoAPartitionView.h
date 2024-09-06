/**
 * @file SoAPartitionView.h
 * @authors humig
 * @date 07.07.2019
 */

#pragma once

namespace autopas {

template <class SoAArraysType>
class SoAPartition;

/**
 * View on a fixed part of a SoAPartition between a start index and an end index.
 *
 * It is the user`s responsibility to ensure that modifications to the underlying SoAPartition don't let the start index get
 * bigger than the highest index in the SoAPartition, or the end index higher than the number of entries in the SoAPartition.
 *
 * @tparam SoAArraysType The SoAArrayType of the SoAPartition.
 */
template <class SoAArraysType>
class SoAPartitionView {
 public:
  /**
   * Default constructor of SoAPartitionView to allow storing it in containers.
   */
  SoAPartitionView() : _soa(nullptr), _startIndex(0), _endIndex(0) {}

  /**
   * Constructs a view on \p soa that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `soa->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param soa The SoAPartition to view.
   * @param startIndex The index of the first entry to view of the SoAPartition.
   * @param endIndex The index of the entry after the last entry to view of the SoAPartition.
   */
  SoAPartitionView(SoAPartition<SoAArraysType> *soa, size_t startIndex, size_t endIndex)
      : _soa(soa), _startIndex(startIndex), _endIndex(endIndex) {
    if (not(soa->size() >= endIndex and endIndex >= startIndex)) /* @todo C++20 [[unlikely]] */ {
      utils::ExceptionHandler::exception("SoAPartitionView: Trying to view particles outside of the SoAPartition.");
    }
  }

  /**
   * Constructs a SoAPartitionView on the whole content of \p soa.
   * @param soa The SoAPartition to view.
   */
  explicit SoAPartitionView(SoAPartition<SoAArraysType> *soa) : _soa(soa), _startIndex(0), _endIndex(soa->size()) {}

  /**
   * Implicit constructor that converts a SoAPartition to SoAPartitionView.
   * @param soa The SoAPartition to view.
   */
  SoAPartitionView(SoAPartition<SoAArraysType> &soa) : _soa(&soa), _startIndex(0), _endIndex(soa.size()) {}

  /**
   * Returns a pointer to the given attribute vector.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector
   */
  template <size_t attribute>
  [[nodiscard]] auto begin() {
    return _soa->template begin<attribute>() + _startIndex;
  }

  /**
   * Returns a pointer to the given attribute vector const.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector const
   */
  template <size_t attribute>
  [[nodiscard]] auto begin() const {
    return _soa->template begin<attribute>() + _startIndex;
  }

  /**
   * Returns the number of particles in the view.
   *
   * @return Number of particles.
   */
  [[nodiscard]] size_t size() const { return _endIndex - _startIndex; }

 private:
  /**
   * The underlying SoAPartition.
   */
  SoAPartition<SoAArraysType> *_soa;

  /**
   * The start index of the view in the SoAPartition. (Inclusive)
   */
  size_t _startIndex;

  /**
   * The end index of the view in the SoAPartition. (Exclusive)
   */
  size_t _endIndex;
};

}  // namespace autopas
