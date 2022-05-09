/**
 * @file SoA.h
 * @authors tchipevn, seckler
 * @date 18.01.2018
 */

#pragma once

#include <algorithm>
#include <initializer_list>
#include <map>
#include <tuple>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoAStorage.h"
#include "autopas/utils/SoAType.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * Structur of the array class.
 * @tparam SoAArraysType The SoAArrayType to be used for storage.
 */
template <class SoAArraysType>
class SoA {
 public:
  /**
   * Default constructor.
   */
  SoA() = default;

  /**
   * Copy constructor.
   * @param soa SoA to copy.
   */
  SoA(const SoA &soa) = default;

  /**
   * Resizes all Vectors to the given length.
   * @param length new length.
   */
  void resizeArrays(size_t length) {
    soaStorage.apply([=](auto &list) { list.resize(length); });
  }

  /**
   * Pushes a given value to the desired attribute array.
   * @tparam attribute Index of array to push to.
   * @param value Value to push.
   */
  template <std::size_t attribute>
  void push(const double value) {
    soaStorage.template get<attribute>().push_back(value);
  }

  /**
   * Writes / updates the value of an attribute for a specific particle.
   * @tparam attribute Attribute to update.
   * @tparam ValueType type of the attribute
   * @param particleId Particle to update.
   * @param value New value.
   */
  template <int attribute, class ValueType>
  void write(size_t particleId, const ValueType &value) {
    soaStorage.template get<attribute>().at(particleId) = value;
  }

  /**
   * Appends the other SoA buffer to this.
   * @param other Other buffer.
   */
  void append(const SoA<SoAArraysType> &other) {
    if (other.getNumberOfParticles() > 0) {
      append_impl(other.soaStorage, std::make_index_sequence<std::tuple_size<SoAArraysType>::value>{});
    }
  }

  /**
   * Appends the other SoA buffer to this.
   * @param other Other buffer.
   */
  void append(const SoAView<SoAArraysType> &other) {
    if (other.getNumberOfParticles() > 0) {
      append_impl(other, std::make_index_sequence<std::tuple_size<SoAArraysType>::value>{});
    }
  }

  /**
   * Appends the specified attributes from the other SoA buffer to this.
   * @tparam attributes Attributes to append.
   * @param other Other buffer.
   */
  template <int... attributes>
  void append(const SoA<SoAArraysType> &other) {
    if (other.getNumberOfParticles() > 0) {
      const auto newSize = getNumberOfParticles() + other.getNumberOfParticles();
      append_impl(other.soaStorage, std::index_sequence<attributes...>{});
      resizeArrays(newSize);
    }
  }

  /**
   * Writes or updates values of attributes for a specific particle.
   * @tparam attributes Array of attributes to update.
   * @tparam ValueArrayType type of the array
   * @param particleId Particle to update.
   * @param values New value.
   */
  template <int... attributes, class ValueArrayType>
  void writeMultiple(size_t particleId, const ValueArrayType &values) {
    write_impl<attributes...>(particleId, values);
  }

  /**
   * Specialized version to pass arrays without specifying it directly.
   * @tparam attributes
   * @tparam N
   * @param particleId
   * @param values
   */
  template <int... attributes, size_t N = sizeof...(attributes)>
  inline void writeMultiple(size_t particleId, std::array<double, N> values) {
    write_impl<attributes...>(particleId, values);
  }

  /**
   * Reads from all given attribute arrays at position `particleId`.
   * @tparam ArrayLength length of the returned array. Should be equal
   * attributes.size().
   * @tparam attributes Attributes to read from.
   * @param particleId Position to read from.
   * @return Array of attributes ordered by given attribute order.
   */
  template <int... attributes>
  std::array<double, sizeof...(attributes)> readMultiple(size_t particleId) const {
    std::array<double, sizeof...(attributes)> retArray;
    if (particleId >= getNumberOfParticles()) {
      autopas::utils::ExceptionHandler::exception(
          "SoA::read: requested particle id ({}) is bigger than number of particles ({})", particleId,
          getNumberOfParticles());
      return retArray;
    }
    read_impl<attributes...>(particleId, retArray);
    return retArray;
  }

  /**
   * Reads the value of a given attribute of a given particle.
   * @tparam attribute Attribute to read from.
   * @param particleId Position to read from.
   * @return Attribute value.
   */
  template <std::size_t attribute>
  auto read(size_t particleId) const {
    return soaStorage.template get<attribute>().at(particleId);
  }

  /**
   * Returns a pointer to the given attribute vector.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector
   */
  template <std::size_t attribute>
  auto begin() {
    return soaStorage.template get<attribute>().data();
  }

  /**
   * Returns the number of particles.
   *
   * This function only checks the size of the first array since it is assumed
   * that the user manages the arrays responsibly.
   *
   * @return Number of particles.
   */
  inline size_t getNumberOfParticles() const {
    size_t maxLength = 0;
    utils::TupleUtils::for_each(soaStorage.getTuple(), [&](auto &v) { maxLength = std::max(maxLength, v.size()); });
    return maxLength;
  }

  /**
   * delete all particles in the soa
   */
  void clear() {
    soaStorage.apply([](auto &list) { list.clear(); });
  }

  /**
   * swap the position of two particles in the soa
   * @param a position of the first particle
   * @param b position of the second particle
   */
  void swap(std::size_t a, std::size_t b) {
    soaStorage.apply([=](auto &list) { std::swap(list[a], list[b]); });
  }

  /**
   * Delete the last particle in the SoA.
   */
  void pop_back() {
    soaStorage.apply([](auto &list) { list.pop_back(); });
  }

  /**
   * Constructs a SoAView for the whole SoA and returns it.
   * @return the constructed SoAView on the whole SoA.
   */
  SoAView<SoAArraysType> constructView() { return {this, 0, getNumberOfParticles()}; }

  /**
   * Constructs a view that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `this->getNumberOfParticles()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param startIndex The index of the first entry to view.
   * @param endIndex The index of the entry after the last entry to view.
   * @return the constructed SoAView from \p startIndex (inclusive) to \p endIndex (exclusive).
   */
  SoAView<SoAArraysType> constructView(size_t startIndex, size_t endIndex) { return {this, startIndex, endIndex}; }

 private:
  // actual implementation of read
  template <int attribute, int... attributes, class ValueArrayType>
  void read_impl(size_t particleId, ValueArrayType &values, int _current = 0) const {
    values[_current] = soaStorage.template get<attribute>().at(particleId);
    read_impl<attributes...>(particleId, values, _current + 1);
  }

  // stop of recursive read call
  template <class ValueArrayType>
  void read_impl(size_t particleId, ValueArrayType &values, int _current = 0) const {}

  // actual implementation of the write function.
  // uses a recursive call.
  template <int attribute, int... attributes, class ValueArrayType>
  void write_impl(size_t particleId, const ValueArrayType &values, int _current = 0) {
    soaStorage.template get<attribute>().at(particleId) = values[_current];
    write_impl<attributes...>(particleId, values, _current + 1);
  }

  // Stop of the recursive write_impl call
  template <class ValueArrayType>
  void write_impl(size_t particleId, const ValueArrayType &values, int _current = 0) {}

  // helper function to append a single array
  template <std::size_t attribute>
  void appendSingleArray(const utils::SoAStorage<SoAArraysType> &valArrays) {
    auto &currentVector = soaStorage.template get<attribute>();
    auto &otherVector = valArrays.template get<attribute>();
    currentVector.insert(currentVector.end(), otherVector.begin(), otherVector.end());
  }

  // helper function to append a single array
  template <std::size_t attribute>
  void appendSingleArray(const SoAView<SoAArraysType> &valArrays) {
    auto &currentVector = soaStorage.template get<attribute>();
    auto otherVectorIterator = valArrays.template begin<attribute>();
    currentVector.insert(currentVector.end(), otherVectorIterator,
                         otherVectorIterator + valArrays.getNumberOfParticles());
  }

  // actual implementation of append
  template <std::size_t... Is>
  void append_impl(const utils::SoAStorage<SoAArraysType> &valArrays, std::index_sequence<Is...>) {
    // fold expression
    (appendSingleArray<Is>(valArrays), ...);
  }

  // actual implementation of append
  template <std::size_t... Is>
  void append_impl(const SoAView<SoAArraysType> &valArrays, std::index_sequence<Is...>) {
    // fold expression
    (appendSingleArray<Is>(valArrays), ...);
  }

  // ------------- members ---------------

  // storage container for the SoA's
  utils::SoAStorage<SoAArraysType> soaStorage;
};
}  // namespace autopas
