/**
 * @file SoAStorage.h
 * @author seckler
 * @date 11.06.18
 */

#pragma once

#include <tuple>

namespace autopas {
namespace utils {

/**
 * SoAStorage is a helper to access the stored SoA's.
 * @tparam SoAArraysType the type of the storage arrays. should be a tuple of aligned vectors
 */
template <class SoAArraysType>
class SoAStorage {
 private:
  // End of iteration.
  template <std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if<I == sizeof...(Tp), void>::type for_each(std::tuple<Tp...>&, FuncT) {}

  template <std::size_t I = 0, typename FuncT, typename... Tp>
      inline typename std::enable_if < I<sizeof...(Tp), void>::type for_each(std::tuple<Tp...>& t, FuncT f) {
    f(std::get<I>(t));
    for_each<I + 1, FuncT, Tp...>(t, f);
  }

 public:
  /**
   * Apply the specific function to all vectors.
   * This can e.g. be resize operations, ...
   * The functor func should not require input parameters. The returns of the functor are ignored.
   * @tparam TFunctor the type of the functor
   * @param func a functor, that should be applied on all vectors (e.g. lambda functions)
   */
  template <typename TFunctor>
  void apply(TFunctor func) {
    for_each(soaStorageTuple, func);
  }

  /**
   * Get the vector at the specific entry of the storage
   * @tparam soaAttribute the attribute for which the vector should be returned
   * @return a reference to the vector for the specific attribute
   */
  template <size_t soaAttribute>
  auto& get() {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  /**
   * @copydoc auto& get()
   * @note const variant
   */
  template <size_t soaAttribute>
  const auto& get() const {
    return std::get<soaAttribute>(soaStorageTuple);
  }

 private:
  SoAArraysType soaStorageTuple;
};

}  // namespace utils
}  // namespace autopas