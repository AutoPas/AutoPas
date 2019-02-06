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
  // End of iteration/recursion.
  template <std::size_t I, std::size_t END, typename FunctorT, typename... Tp>
  inline typename std::enable_if<I == END, void>::type for_each(std::tuple<Tp...>&, FunctorT) {}

  template <std::size_t I, std::size_t END, typename FunctorT, typename... Tp>
  inline typename std::enable_if<(I < END), void>::type for_each(std::tuple<Tp...>& t, FunctorT f) {
    f(std::get<I>(t));
    for_each<I + 1, END, FunctorT, Tp...>(t, f);
  }

  template <typename FunctorT, typename... Tp>
  inline void for_each_wrap(std::tuple<Tp...>& t, FunctorT f) {
    for_each<0, sizeof...(Tp), FunctorT, Tp...>(t, f);
  }

  template <std::size_t END, typename FunctorT, typename... Tp>
  inline void for_each_wrap(std::tuple<Tp...>& t, FunctorT f) {
    for_each<0, END, FunctorT, Tp...>(t, f);
  }

 public:
  /**
   * Apply the specific function to all vectors.
   * This can e.g. be resize operations, ...
   * The functor func should not require input parameters. The returns of the functor are ignored.
   * @tparam FunctorT The type of the functor
   * @param func A functor, that should be applied on all vectors (e.g. lambda functions, should take `auto& list` as an
   * argument).
   */
  template <typename FunctorT>
  void apply(FunctorT func) {
    for_each_wrap(soaStorageTuple, func);
  }

  /**
   * Apply the specific function to all vectors.
   * This can e.g. be resize operations, ...
   * The functor func should not require input parameters. The returns of the functor are ignored.
   * @tparam FunctorT The type of the functor.
   * @tparam numElements Number of elements to apply func to.
   * @param func A functor, that should be applied on all vectors (e.g. lambda functions, should take `auto& list` as an
   * argument).
   */
  template <std::size_t numElements, typename FunctorT>
  void apply(FunctorT func) {
    for_each_wrap<numElements>(soaStorageTuple, func);
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
   * @copydoc get()
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