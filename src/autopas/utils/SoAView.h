/**
 * @file SoAView.h
 * @authors humig
 * @date 07.07.2019
 */

#pragma once

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
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `soa->getNumberOfParticles()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param soa The SoA to view.
   * @param startIndex The index of the first entry to view of the SoA.
   * @param endIndex The index of the entry after the last entry to view of the SoA.
   */
  SoAView(SoA<SoAArraysType> *soa, size_t startIndex, size_t endIndex)
      : _soa(soa), _startIndex(startIndex), _endIndex(endIndex) {
    if (not(soa->getNumberOfParticles() >= endIndex and endIndex >= startIndex)) /* @todo C++20 [[unlikely]] */ {
      utils::ExceptionHandler::exception("SoAView: Trying to view particles outside of the SoA.");
    }
  }

  /**
   * Constructs a SoAView on the whole content of \p soa.
   * @param soa The SoA to view.
   */
  explicit SoAView(SoA<SoAArraysType> *soa) : _soa(soa), _startIndex(0), _endIndex(soa->getNumberOfParticles()) {}

  /**
   * Implicit constructor that converts a SoA to SoAView.
   * @param soa The SoA to view.
   */
  SoAView(SoA<SoAArraysType> &soa) : _soa(&soa), _startIndex(0), _endIndex(soa.getNumberOfParticles()) {}

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
  [[nodiscard]] size_t getNumberOfParticles() const { return _endIndex - _startIndex; }

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
