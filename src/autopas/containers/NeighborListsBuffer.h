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
 * Class for manual memory management of neighbor lists.
 *
 * The key, value system behaves similar to a std::unordered_map, so see the
 * [official doc](https://en.cppreference.com/w/cpp/container/unordered_map) regarding these details.
 *
 * A internal buffer of lists is kept in memory, never deleted, and only its content cleared and reassigned.
 * Get access to new lists via the index returned from getNewNeighborList() and work on them via getNeighborListRef()
 * .
 * @tparam Key Preferably something that is cheap to destroy.
 * @tparam Value Has to be default constructable.
 * @tparam Hash
 * @tparam KeyEqual
 */
template <class Key, class Value, class Hash = std::hash<Key>, class KeyEqual = std::equal_to<Key>>
class NeighborListsBuffer {
 public:
  /**
   * Getter for a reference to a neighbor list by index.
   * @tparam throwIfIndexOutOfBounds
   * @param index
   * @return
   */
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

  /**
   * Getter for a reference to a neighbor list by key.
   * @tparam throwIfKeyIsUnknown
   * @param key
   * @return
   */
  template <bool throwIfKeyIsUnknown = false>
  std::vector<Value> &getNeighborListRef(const Key &key) {
    if constexpr (throwIfKeyIsUnknown) {
      if (_keyMap.find(key) == _keyMap.end()) {
        autopas::utils::ExceptionHandler::exception("NeighborListsBuffer::getNeighborListRef() unknown key: {}", key);
      }
    }
    return _neighborLists[_keyMap[key]];
  }

  /**
   * Reserves a neighbor list for use.
   * Grows the internal buffer if there are no spare lists left.
   *
   * @note This function is NOT thread safe.
   *
   * @return Index of the newly reserved list.
   */
  size_t getNewNeighborList() {
    // if the buffer is saturated...
    if (_lastValidListIndex >= _neighborLists.size()) {
      // ...grow it and initialize all new lists. Make sure to grow to at least 10 lists.
      _neighborLists.resize(std::max(10ul, static_cast<size_t>(_neighborLists.size() * _growthFactor)),
                            std::vector<Value>(_defaultListLength));
    }
    return ++_lastValidListIndex;
  }

  /**
   * Assigns a neighbor list to the given key.
   *
   * @note This function is NOT thread safe.
   *
   * @param key
   * @return Index of the newly reserved list.
   */
  size_t getNewNeighborList(const Key &key) {
    const auto newIndex = getNewNeighborList();
    _keyMap.emplace(key, newIndex);
    return newIndex;
  }

  /**
   * Clears the internal key map and moves _lastValidListIndex to indicate an empty buffer.
   * This function does not touch the buffer itself.
   */
  void clear() {
    _keyMap.clear();
    _lastValidListIndex = std::numeric_limits<size_t>::max();
  }

  /**
   * Resize the internal buffer so that there are new spare lists.
   * If n is not larger than the current capacity nothing happens.
   *
   * @param n Number of lists to allocate space for.
   */
  void reserveNeighborLists(size_t n) { _neighborLists.resize(n, std::vector<Value>(_defaultListLength)); }

  /**
   * Set the initial length of new neighbor lists.
   * @param defaultListLength
   */
  void setDefaultListLength(size_t defaultListLength) { NeighborListsBuffer::_defaultListLength = defaultListLength; }

  /**
   * Set the growth factor for the internal buffer.
   * @param growthFactor
   */
  void setGrowthFactor(double growthFactor) { NeighborListsBuffer::_growthFactor = growthFactor; }

 private:
  /**
   * Map to associate keys with indexes of _neighborLists.
   */
  std::unordered_map<Key, size_t, Hash, KeyEqual> _keyMap{};
  /**
   * Continuous but unordered collection buffer of neighbor lists that is never supposed to be deallocated during the
   * lifetime of this object.
   * The idea is to always have capacity == size to have spare vectors to assign. The actual size is controlled via
   * _lastValidListIndex.
   */
  std::vector<std::vector<Value>> _neighborLists{};
  /**
   * Anything beyond this index in _neighborLists is considered deleted or unassigned.
   */
  size_t _lastValidListIndex{std::numeric_limits<size_t>::max()};
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