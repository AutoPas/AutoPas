/**
 * @file SoAHelpers.h
 * @author seckler
 * @date 11.06.18
 */

#pragma once

#include <tuple>

namespace autopas {
namespace utils {

template <class SoAArraysType>
class SoAStorage {
 private:
  // Unused arguments are given no names.
  template <std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if<I == sizeof...(Tp), void>::type for_each(std::tuple<Tp...>&, FuncT) {}

  template <std::size_t I = 0, typename FuncT, typename... Tp>
      inline typename std::enable_if < I<sizeof...(Tp), void>::type for_each(std::tuple<Tp...>& t, FuncT f) {
    f(std::get<I>(t));
    for_each<I + 1, FuncT, Tp...>(t, f);
  }

 public:
  template <typename TFunctor>
  void apply(TFunctor func) {
    for_each(soaStorageTuple, func);
  }

  template <size_t soaAttribute>
  auto& get() {
    return std::get<soaAttribute>(soaStorageTuple);
  }

  template <size_t soaAttribute>
  const auto& get() const {
    return std::get<soaAttribute>(soaStorageTuple);
  }

 private:
  SoAArraysType soaStorageTuple;
};

}  // namespace utils
}  // namespace autopas