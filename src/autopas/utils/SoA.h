/**
 * @file SoA.h
 * @authors tchipevn, seckler
 * @date 18.01.2018
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/AutoPasMacros.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoAStorage.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * Structure of the array class.
 * @tparam SoAArraysType The SoAArrayType to be used for storage.
 */
template <class SoAArraysType>
class SoA {
 public:
  /**
   * Number of base attributes defined in SoAArraysType.
   */
  static constexpr size_t NumBaseAttrs = std::tuple_size_v<SoAArraysType>;

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
   * Copy assignment operator.
   * @param soa SoA to copy.
   * @return *this
   */
  SoA &operator=(const SoA &soa) = default;

  /**
   * Move constructor.
   * @param soa SoA to move from.
   */
  SoA(SoA &&soa) noexcept = default;

  /**
   * Move assignment operator.
   * @param soa SoA to move from.
   * @return *this
   */
  SoA &operator=(SoA &&soa) noexcept = default;

  /**
   * Resizes all base vectors and all existing extra columns to the given length.
   * @param length new length.
   */
  void resizeArrays(size_t length) {
    _size = length;
    soaStorage.apply([=](auto &list) { list.resize(length); });
    for (auto &col : _extraColumns) {
      col.resize(length);
    }
  }

  /**
   * Resizes only the arrays for the specified attributes, allocating extra columns as needed.
   * @tparam attributes Attributes to resize.
   * @param length new length.
   */
  template <size_t... attributes>
  void resizeArraysForAttributes(const size_t length) {
    _size = length;
    if constexpr (sizeof...(attributes) > 0) {
      constexpr size_t maxAttr = std::max({size_t{0}, attributes...});
      if constexpr (maxAttr >= NumBaseAttrs) {
        ensureExtraColumnCount(maxAttr - NumBaseAttrs + 1);
      }
      (resizeSingleAttributeUnchecked<attributes>(length), ...);
    }
  }

  /**
   * Resizes a single attribute array to the given length.
   * @tparam attribute Attribute to resize.
   * @param length new length.
   */
  template <size_t attribute>
  void resizeSingleAttribute(size_t length) {
    getOrCreateColumnVector<attribute>().resize(length);
  }

  /**
   * Pushes a given value to the desired attribute array.
   * @tparam attribute Index of array to push to.
   * @param value Value to push.
   */
  template <size_t attribute>
  void push(const double value) {
    auto &column = getOrCreateColumnVector<attribute>();
    column.push_back(value);
    _size = std::max(_size, column.size());
  }

  /**
   * Writes / updates the value of an attribute for a specific particle.
   * @tparam attribute Attribute to update.
   * @tparam ValueType type of the attribute
   * @param particleId Particle to update.
   * @param value New value.
   */
  template <size_t attribute, class ValueType>
  void write(size_t particleId, const ValueType &value) {
    if constexpr (attribute < NumBaseAttrs) {
      soaStorage.template get<attribute>().at(particleId) = value;
    } else {
      _extraColumns.at(attribute - NumBaseAttrs).at(particleId) = value;
    }
  }

  /**
   * Appends the other SoA buffer to this.
   * @param other Other buffer.
   */
  void append(const SoA<SoAArraysType> &other) {
    if (other.size() > 0) {
      _size += other.size();
      append_impl(other.soaStorage, std::make_index_sequence<NumBaseAttrs>{});
      ensureExtraColumnCount(other._extraColumns.size());
      for (size_t col = 0; col < other._extraColumns.size(); ++col) {
        _extraColumns[col].insert(_extraColumns[col].end(), other._extraColumns[col].begin(),
                                  other._extraColumns[col].end());
      }
    }
  }

  /**
   * Appends the other SoA buffer to this.
   * Handles partially sized columns safely by only appending columns that contain data in the view.
   * @param other Other buffer.
   */
  void append(const SoAView<SoAArraysType> &other) {
    if (other.size() > 0) {
      _size += other.size();
      append_impl(other, std::make_index_sequence<NumBaseAttrs>{});
      const size_t numExtra = other.numExtraColumns();
      ensureExtraColumnCount(numExtra);
      for (size_t col = 0; col < numExtra; ++col) {
        if (const auto otherCol = other.extraColumn(col); not otherCol.empty()) {
          _extraColumns[col].insert(_extraColumns[col].end(), otherCol.begin(), otherCol.end());
        }
      }
    }
  }

  /**
   * Appends the specified attributes from the other SoA buffer to this.
   * @tparam attributes Attributes to append.
   * @param other Other buffer.
   */
  template <size_t... attributes>
  void append(const SoA<SoAArraysType> &other) {
    if (other.size() > 0) {
      const auto newSize = size() + other.size();
      (appendSingleArray<attributes>(other), ...);
      _size = newSize;
    }
  }

  /**
   * Appends the specified attributes from the other SoAView to this.
   * @tparam attributes Attributes to append.
   * @param other Other buffer.
   */
  template <size_t... attributes>
  void append(const SoAView<SoAArraysType> &other) {
    if (other.size() > 0) {
      const auto newSize = size() + other.size();
      (appendSingleArray<attributes>(other), ...);
      _size = newSize;
    }
  }

  /**
   * Writes or updates values of attributes for a specific particle.
   * @tparam attributes Array of attributes to update.
   * @tparam ValueArrayType type of the array
   * @param particleId Particle to update.
   * @param values New value.
   */
  template <size_t... attributes, class ValueArrayType>
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
  template <size_t... attributes, size_t N = sizeof...(attributes)>
  void writeMultiple(size_t particleId, const std::array<double, N> &values) {
    write_impl<attributes...>(particleId, values);
  }

  /**
   * Reads from all given attribute arrays at position `particleId`.
   * @tparam attributes Attributes to read from.
   * @param particleId Position to read from.
   * @return Array of attributes ordered by given attribute order.
   */
  template <size_t... attributes>
  std::array<double, sizeof...(attributes)> readMultiple(size_t particleId) const {
    std::array<double, sizeof...(attributes)> retArray;
    if (particleId >= size()) {
      utils::ExceptionHandler::exception(
          "SoA::read: requested particle id ({}) is bigger than number of particles ({})", particleId, size());
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
  template <size_t attribute>
  auto read(size_t particleId) const {
    if constexpr (attribute < NumBaseAttrs) {
      return soaStorage.template get<attribute>().at(particleId);
    } else {
      return _extraColumns.at(attribute - NumBaseAttrs).at(particleId);
    }
  }

  /**
   * Returns a pointer to the given attribute vector.
   * Caller must ensure any extra-column attribute has been allocated.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector.
   */
  template <size_t attribute>
  [[nodiscard]] AUTOPAS_ALWAYS_INLINE auto begin() {
    if constexpr (attribute < NumBaseAttrs) {
      return soaStorage.template get<attribute>().data();
    } else {
      assert(_extraColumns.size() > attribute - NumBaseAttrs);
      return _extraColumns[attribute - NumBaseAttrs].data();
    }
  }

  /**
   * Returns a const pointer to the given attribute vector.
   * Caller must ensure any extra-column attribute has been allocated.
   * @tparam attribute ID of the desired attribute.
   * @return Pointer to the beginning of the attribute vector.
   */
  template <size_t attribute>
  [[nodiscard]] AUTOPAS_ALWAYS_INLINE auto begin() const {
    if constexpr (attribute < NumBaseAttrs) {
      return soaStorage.template get<attribute>().data();
    } else {
      assert(_extraColumns.size() > attribute - NumBaseAttrs);
      return _extraColumns[attribute - NumBaseAttrs].data();
    }
  }

  /**
   * Returns the actual allocated size of the array for the given attribute.
   * @tparam attribute ID of the desired attribute.
   * @return Number of elements in the attribute's array.
   */
  template <size_t attribute>
  [[nodiscard]] size_t size() const {
    if constexpr (attribute < NumBaseAttrs) {
      return soaStorage.template get<attribute>().size();
    } else {
      if (_extraColumns.size() <= attribute - NumBaseAttrs) {
        return 0;
      }
      return _extraColumns[attribute - NumBaseAttrs].size();
    }
  }

  /**
   * Returns the number of particles.
   * @return Number of particles.
   */
  [[nodiscard]] AUTOPAS_ALWAYS_INLINE size_t size() const { return _size; }

  /**
   * Returns the number of dynamically allocated extra columns.
   * @return Number of extra columns.
   */
  [[nodiscard]] size_t numExtraColumns() const { return _extraColumns.size(); }

  /**
   * Returns a const span over the elements of the given extra column.
   * @param col Index of the extra column (0-based).
   * @return Span of elements, or an empty span if col is out of range.
   */
  [[nodiscard]] std::span<const double> extraColumn(size_t col) const {
    if (col < _extraColumns.size()) {
      return _extraColumns[col];
    }
    return {};
  }

  /**
   * Returns the allocated size of the given extra column.
   * @param col Index of the extra column (0-based).
   * @return Number of elements in the extra column, or 0 if col is out of range.
   */
  [[nodiscard]] size_t extraColumnSize(size_t col) const { return extraColumn(col).size(); }

  /**
   * Returns a const pointer to the beginning of the given extra column.
   * @param col Index of the extra column (0-based).
   * @return Pointer to data, or nullptr if empty or out of range.
   */
  [[nodiscard]] const double *extraColumnBegin(size_t col) const {
    const auto colSpan = extraColumn(col);
    return colSpan.empty() ? nullptr : colSpan.data();
  }

  /**
   * Returns a pointer to the beginning of the given extra column.
   * @param col Index of the extra column (0-based).
   * @return Pointer to data, or nullptr if empty or out of range.
   */
  [[nodiscard]] double *extraColumnBegin(size_t col) {
    return const_cast<double *>(std::as_const(*this).extraColumnBegin(col));
  }

  /**
   * delete all particles in the soa
   */
  void clear() {
    _size = 0;
    soaStorage.apply([](auto &list) { list.clear(); });
    for (auto &col : _extraColumns) {
      col.clear();
    }
  }

  /**
   * swap the position of two particles in the soa
   * @param a position of the first particle
   * @param b position of the second particle
   */
  void swap(size_t a, size_t b) {
    soaStorage.apply([=](auto &list) {
      if (not list.empty()) {
        std::swap(list[a], list[b]);
      }
    });
    for (auto &col : _extraColumns) {
      if (not col.empty()) {
        std::swap(col[a], col[b]);
      }
    }
  }

  /**
   * Delete the last particle in the SoA.
   */
  void pop_back() {
    if (_size > 0) {
      --_size;
    }
    soaStorage.apply([](auto &list) {
      if (not list.empty()) {
        list.pop_back();
      }
    });
    for (auto &col : _extraColumns) {
      if (not col.empty()) {
        col.pop_back();
      }
    }
  }

  /**
   * Constructs a SoAView for the whole SoA and returns it.
   * @return the constructed SoAView on the whole SoA.
   */
  AUTOPAS_ALWAYS_INLINE SoAView<SoAArraysType> constructView() { return SoAView<SoAArraysType>{this}; }

  /**
   * Constructs a view that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `this->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param startIndex The index of the first entry to view.
   * @param endIndex The index of the entry after the last entry to view.
   * @return the constructed SoAView from \p startIndex (inclusive) to \p endIndex (exclusive).
   */
  AUTOPAS_ALWAYS_INLINE SoAView<SoAArraysType> constructView(size_t startIndex, size_t endIndex) {
    return {this, startIndex, endIndex};
  }

 private:
  // actual implementation of read
  template <size_t attribute, size_t... attributes, class ValueArrayType>
  void read_impl(size_t particleId, ValueArrayType &values, int _current = 0) const {
    if constexpr (attribute < NumBaseAttrs) {
      values[_current] = soaStorage.template get<attribute>().at(particleId);
    } else {
      values[_current] = _extraColumns.at(attribute - NumBaseAttrs).at(particleId);
    }
    read_impl<attributes...>(particleId, values, _current + 1);
  }

  // stop of recursive read call
  template <class ValueArrayType>
  void read_impl(size_t particleId, ValueArrayType &values, int _current = 0) const {}

  // actual implementation of the write function.
  // uses a recursive call.
  template <size_t attribute, size_t... attributes, class ValueArrayType>
  void write_impl(size_t particleId, const ValueArrayType &values, int _current = 0) {
    if constexpr (attribute < NumBaseAttrs) {
      soaStorage.template get<attribute>().at(particleId) = values[_current];
    } else {
      _extraColumns.at(attribute - NumBaseAttrs).at(particleId) = values[_current];
    }
    write_impl<attributes...>(particleId, values, _current + 1);
  }

  // Stop of the recursive write_impl call
  template <class ValueArrayType>
  void write_impl(size_t particleId, const ValueArrayType &values, int _current = 0) {}

  // ensures _extraColumns holds at least minColumns vectors
  void ensureExtraColumnCount(const size_t minColumns) {
    if (_extraColumns.size() < minColumns) {
      _extraColumns.resize(minColumns);
    }
  }

  // returns reference to the vector backing attribute, allocating extra columns if needed
  template <size_t attribute>
  auto &getOrCreateColumnVector() {
    if constexpr (attribute < NumBaseAttrs) {
      return soaStorage.template get<attribute>();
    } else {
      ensureExtraColumnCount(attribute - NumBaseAttrs + 1);
      return _extraColumns[attribute - NumBaseAttrs];
    }
  }

  // resizes a single attribute array, assuming any required extra column already exists
  template <size_t attribute>
  void resizeSingleAttributeUnchecked(size_t length) {
    if constexpr (attribute < NumBaseAttrs) {
      soaStorage.template get<attribute>().resize(length);
    } else {
      assert(_extraColumns.size() > attribute - NumBaseAttrs);
      _extraColumns[attribute - NumBaseAttrs].resize(length);
    }
  }

  // helper function to append a single base array from SoAStorage
  template <size_t attribute>
  void appendSingleBaseArray(const utils::SoAStorage<SoAArraysType> &valArrays) {
    auto &currentVector = soaStorage.template get<attribute>();
    auto &otherVector = valArrays.template get<attribute>();
    currentVector.insert(currentVector.end(), otherVector.begin(), otherVector.end());
  }

  // helper function to append a single array from a SoA or SoAView
  template <size_t attribute, class SoASource_T>
  void appendSingleArray(const SoASource_T &other) {
    const size_t count = other.template size<attribute>();
    if (count == 0) {
      return;
    }
    auto &currentVector = getOrCreateColumnVector<attribute>();
    const auto otherBegin = other.template begin<attribute>();
    currentVector.insert(currentVector.end(), otherBegin, otherBegin + count);
  }

  // actual implementation of append
  template <size_t... Is>
  void append_impl(const utils::SoAStorage<SoAArraysType> &valArrays, std::index_sequence<Is...>) {
    (appendSingleBaseArray<Is>(valArrays), ...);
  }

  // actual implementation of append
  template <size_t... Is>
  void append_impl(const SoAView<SoAArraysType> &valArrays, std::index_sequence<Is...>) {
    (appendSingleArray<Is>(valArrays), ...);
  }

  // ------------- members ---------------

  // storage container for base SoA's
  utils::SoAStorage<SoAArraysType> soaStorage;

  // storage container for dynamically allocated extra double columns
  std::vector<std::vector<double, AlignedAllocator<double>>> _extraColumns;

  // tracked size of the SoA
  size_t _size{0};
};
}  // namespace autopas
