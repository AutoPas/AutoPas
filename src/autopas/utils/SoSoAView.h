/**
 * @file SoSoAView.h
 * @authors S. J. Newcome
 * @date 06/09/2024
 */

#pragma once

#include <tuple>
#include <vector>
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

template <typename... SoAtypes>
class SoSoA;

/**
 * View on a fixed part of a SoSoA between a start index and an end index.
 *
 * It is the user`s responsibility to ensure that modifications to the underlying SoSoA don't result in the start index
 * becoming higher than the highest index in the SoSoA, or the end index higher than the number of entries in the SoSoA.
 *
 * @tparam SoSoAType The SoSoAType of the SoSoA
 */
template <class SoSoAType>
class SoSoAView {
 public:
  /**
   * Default constructor of SoSoAView to allow storing it in containers.
   */
  SoSoAView() : _sosoa(nullptr), _startIndex(0), _endIndex(0) {}

  /**
   * Constructs a view on \p sosoa that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `sosoa->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param sosoa The SoSoA to view.
   * @param startIndex The index of the first entry of the SoSoAView.
   * @param endIndex The index of the entry after the last entry of the SoSoAView.
   */
  SoSoAView(SoSoAType *sosoa, std::size_t startIndex, std::size_t endIndex)
      : _sosoa(sosoa), _startIndex(startIndex), _endIndex(endIndex) {
    if (not(sosoa->size() >= endIndex and endIndex >= startIndex)) /* @todo C++20 [[unlikely]] */ {
      autopas::utils::ExceptionHandler::exception("SoSoAView: Trying to view particles outside of the SoSoA.");
    }
  }

  /**
   * Constructs a SoSoAView on the whole content of \p sosoa.
   * @param sosoa The SoSoA to view.
   */
  explicit SoSoAView(SoSoAType *sosoa) : _sosoa(sosoa), _startIndex(0), _endIndex(sosoa->size()) {}

  /**
   * Implicit constructor that converts a SoSoA to SoSoAView.
   * @param sosoa The SoSoA to view.
   */
  SoSoAView(SoSoAType &sosoa) : _sosoa(&sosoa), _startIndex(0), _endIndex(sosoa.size()) {}

  /**
   * Returns a pointer to the given attribute vector const.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector const
   */
  template <size_t soaTypeIndex, size_t attribute>
  [[nodiscard]] auto begin(size_t depth) {
    return _sosoa->template begin<soaTypeIndex, attribute>(depth);
  }

  /**
   * Returns the number of particles in the view.
   *
   * @return Number of particles.
   */
  [[nodiscard]] size_t size() const { return _endIndex - _startIndex; }

 private:
  /**
   * The underlying SoSoA.
   */
  SoSoAType *_sosoa;

  /**
   * The start index of the view in the SoSoA. (Inclusive)
   */
  size_t _startIndex;

  /**
   * The end index of the view in the SoSoA. (Exclusive)
   */
  size_t _endIndex;
};

}  // namespace autopas
