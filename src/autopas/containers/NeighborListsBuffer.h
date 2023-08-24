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
  template <bool throwOnError>
  std::vector<Value> &getNeighborListRef(size_t index) {
    if constexpr (throwOnError) {
      if (index > lastValidListIndex) {
        autopas::utils::ExceptionHandler::exception(
            "NeighborListsBuffer::getNeighborListRef() index out of bounds: {} > {}", index, lastValidListIndex);
      }
    }
    return neighborLists[index];
  }

  template <bool throwOnError>
  //  std::vector<Value> &getNeighborListRef(const std::conditional<std::is_pointer_v<Key>, Key, Key &> key) {
  std::vector<Value> &getNeighborListRef(const Key &key) {
    if constexpr (throwOnError) {
      if (keyMap.find(key) == keyMap.end()) {
        autopas::utils::ExceptionHandler::exception("NeighborListsBuffer::getNeighborListRef() unknown key: {}", key);
      }
    }
    return neighborLists[keyMap[key]];
  }

  size_t addNeighborList() {
    // if the buffer is saturated...
    if (lastValidListIndex == neighborLists.size()) {
      // grow it and initialize all new lists
      neighborLists.resize(neighborLists.size() * growthFactor, std::vector<Value>(defaultListLength));
    }
    return ++lastValidListIndex;
  }

  //  size_t addNeighborList(const std::conditional<std::is_pointer_v<Key>, Key, Key &> key) {
  size_t addNeighborList(const Key &key) {
    const auto newIndex = addNeighborList();
    keyMap.emplace(key, newIndex);
    return newIndex;
  }

  void clear() {
    keyMap.clear();
    lastValidListIndex = 0;
  }

 private:
  /**
   * Map to associate keys with indexes of neighborLists
   */
  std::unordered_map<Key, size_t, Hash, KeyEqual> keyMap{};
  /**
   * Continuous but unordered collection of neighbor lists that is never supposed to be deallocated during the
   * lifetime of this object.
   */
  std::vector<std::vector<Value>> neighborLists{};
  /**
   * Anything beyond this index in neighborLists is considered deleted.
   */
  size_t lastValidListIndex{0ul};
  /**
   * Number of elements allocated by default in every new neighbor list.
   */
  size_t defaultListLength;
  /**
   * Growth factor for the neighborLists collection if more memory is necessary.
   */
  int growthFactor;
};
}  // namespace autopas