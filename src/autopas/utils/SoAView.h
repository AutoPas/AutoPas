/**
 * @file SoAView.h
 * @authors S. J. Newcome
 * @date 06/09/2024
 */

#pragma once

#include <tuple>
#include <vector>
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <class SoAType>
class SoA;

/**
 * View on a fixed part of a SoA between a start index and an end index.
 *
 * It is the user`s responsibility to ensure that modifications to the underlying SoA don't result in the start index
 * becoming higher than the highest index in the SoA, or the end index higher than the length of the SoA.
 *
 * @tparam SoAType The SoAType of the SoA
 */
template <class SoAType>
class SoAView {
 public:
  /**
   * Default constructor of SoAView to allow storing it in containers.
   */
  SoAView() : _soa(nullptr), _startIndex(0), _endIndex(0) {}

  /**
   * Constructs a view on \p SoA that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `SoA->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param soa The SoA to view.
   * @param startIndex The index of the first entry of the SoAView.
   * @param endIndex The index of the entry after the last entry of the SoAView.
   */
  SoAView(SoAType *soa, std::size_t startIndex, std::size_t endIndex)
      : _soa(soa), _startIndex(startIndex), _endIndex(endIndex) {
    if (not(soa->size() >= endIndex and endIndex >= startIndex)) /* @todo C++20 [[unlikely]] */ {
      autopas::utils::ExceptionHandler::exception("SoAView: Trying to view particles outside of the SoA.");
    }
  }

  /**
   * Constructs a SoAView on the whole content of \p SoA.
   * @param soa The SoA to view.
   */
  explicit SoAView(SoAType *soa) : _soa(soa), _startIndex(0), _endIndex(soa->size()) {}

  /**
   * Implicit constructor that converts a SoA to SoAView.
   * @param soa The SoA to view.
   */
  SoAView(SoAType &soa) : _soa(&soa), _startIndex(0), _endIndex(soa.size()) {}

  /**
   * Returns a pointer to the given attribute vector in the main SoA partition.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector const
   */
  template <size_t attribute>
  auto begin() {
    return _soa->template begin<attribute>();
  }

  /**
   * Returns a pointer to the given attribute vector in the additional SoA partitions.
   * @tparam additionalPartitionType index corresponding to the desired additional SoA partition's type.
   * @tparam attribute desired attribute index
   * @param depth depth of desired SoA partition
   * @return Pointer to the beginning of the attribute vector const
   */
  template <size_t additionalPartitionType, size_t attribute>
  [[nodiscard]] auto begin(size_t depth) {
    return _soa->template begin<additionalPartitionType, attribute>(depth);
  }

  /**
   * Returns the number of particles in the view.
   *
   * @return Number of particles.
   */
  [[nodiscard]] size_t size() const { return _endIndex - _startIndex; }

 private:
  /**
   * The underlying SoA.
   */
  SoAType *_soa;

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
