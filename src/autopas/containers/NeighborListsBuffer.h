/**
 * @file NeighborListsBuffer.h
 * @author F. Gratl
 * @date 24.08.23
 */

#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {
/**
 *
 * @tparam Key Preferably something that is cheap to destroy.
 * @tparam Value Has to be default constructable.
 * @tparam Hash
 * @tparam KeyEqual
 */
template <class Key, class Value, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key>>
class NeighborListsBuffer {
 public:
  template <bool throwIfIndexOutOfBounds = false>
  std::vector<Value> &getNeighborListRef(size_t index) {
    if constexpr (throwIfIndexOutOfBounds) {
      if (index > _lastValidListIndex) {
        autopas::utils::ExceptionHandler::exception(
            "NeighborListsBuffer::getNeighborListRef() index out of bounds: {} > {}", index, _lastValidListIndex);
      }
    }
    return _neighborLists[index];
  }

  template <bool throwIfKeyIsUnknown = false>
  std::vector<Value> &getNeighborListRef(const Key &key) {
    if constexpr (throwIfKeyIsUnknown) {
      if (_keyMap.find(key) == _keyMap.end()) {
        autopas::utils::ExceptionHandler::exception("NeighborListsBuffer::getNeighborListRef() unknown key: {}", key);
      }
    }
    return _neighborLists[_keyMap[key]];
  }

  size_t addNeighborList() {
    // TODO: find a better solution than critical. Maybe id pool per thread?
//#pragma omp critical
//    {
      // if the buffer is saturated...
      if (_lastValidListIndex == _neighborLists.size()) {
        // grow it and initialize all new lists
        _neighborLists.resize(_neighborLists.size() * _growthFactor, std::vector<Value>(_defaultListLength));
      }
      return ++_lastValidListIndex;
//    }
  }

  //  size_t addNeighborList(const std::conditional<std::is_pointer_v<Key>, Key, Key &> key) {
  size_t addNeighborList(const Key &key) {
    const auto newIndex = addNeighborList();
    _keyMap.emplace(key, newIndex);
    return newIndex;
  }

  void clear() {
    _keyMap.clear();
    _lastValidListIndex = 0;
  }

  void reserveNeighborLists(size_t n) { _neighborLists.reserve(n); }

  void setDefaultListLength(size_t defaultListLength) { NeighborListsBuffer::_defaultListLength = defaultListLength; }

  void setGrowthFactor(double growthFactor) { NeighborListsBuffer::_growthFactor = growthFactor; }

 private:
  /**
   * Map to associate keys with indexes of _neighborLists.
   */
  std::unordered_map<Key, size_t, Hash, KeyEqual> _keyMap{};
  /**
   * Continuous but unordered collection of neighbor lists that is never supposed to be deallocated during the
   * lifetime of this object.
   */
  std::vector<std::vector<Value>> _neighborLists{};
  /**
   * Anything beyond this index in _neighborLists is considered deleted.
   */
  size_t _lastValidListIndex{0ul};
  /**
   * Number of elements allocated by default in every new neighbor list.
   */
  size_t _defaultListLength{10};
  /**
   * Growth factor for the _neighborLists collection if more memory is necessary.
   */
  double _growthFactor{2};
};
}  // namespace autopas