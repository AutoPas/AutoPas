/**
 * @file SoAStorage.h
 * @author seckler
 * @date 11.06.18
 */

#pragma once

#include <tuple>

#include "autopas/utils/TupleUtils.h"

namespace autopas::utils {

/**
 * SoAStorage is a helper to access the stored SoA's.
 * @tparam SoAArraysType the type of the storage arrays. should be a tuple of aligned vectors
 */
template <class SoAArraysType>
class SoAStorage {
 public:
  /**
   * Apply the specific function to all vectors.
   * This can e.g. be resize operations, ...
   * The functor func should not require input parameters. The returns of the functor are ignored.
   * @tparam FunctorT the type of the functor
   * @param func a functor, that should be applied on all vectors (e.g. lambda functions, should take `auto& list` as an
   * argument)
   * @todo c++20: replace with expansion statement: `for... (auto& elem : tup) {}`
   */
  template <typename FunctorT>
  void apply(FunctorT func) {
    utils::TupleUtils::for_each(soaStorageTuple, func);
  }

  /**
   * Get the vector at the specific entry of the storage
   * @tparam soaAttribute the attribute for which the vector should be returned
   * @return a reference to the vector for the specific attribute
   */
  template <size_t soaAttribute>
  inline constexpr auto &get() {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  /**
   * @copydoc get()
   * @note const variant
   */
  template <size_t soaAttribute>
  inline constexpr const auto &get() const {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  /**
   * Getter for the internal tuple.
   * @return
   */
  inline constexpr const auto &getTuple() const { return soaStorageTuple; }

 private:
  SoAArraysType soaStorageTuple;
};

}  // namespace autopas::utils