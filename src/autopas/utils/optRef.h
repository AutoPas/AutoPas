/**
 * @file optRef.h
 * @date 11/02/2026
 * @author S. J. Newcome
 */
#pragma once

#include <functional>
#include <optional>

namespace autopas::utils {
/**
 * Short alias for std::optional<std::reference_wrapper<T>>
 */
template <typename T>
using optRef = std::optional<std::reference_wrapper<T>>;

template <typename T>
auto begin(optRef<T> &opt) {
  using Iter = decltype(std::declval<T &>().begin());
  return opt ? opt->get().begin() : Iter{};
}

template <typename T>
auto begin(const optRef<T> &opt) {
  using Iter = decltype(std::declval<T &>().begin());
  return opt ? opt->get().begin() : Iter{};
}

template <typename T>
auto end(optRef<T> &opt) {
  using Iter = decltype(std::declval<T &>().end());
  return opt ? opt->get().end() : Iter{};
}

template <typename T>
auto end(const optRef<T> &opt) {
  using Iter = decltype(std::declval<T &>().end());
  return opt ? opt->get().end() : Iter{};
}

template <typename T>
auto size(optRef<T> &opt) {
  return opt ? opt->get().size() : typename T::size_type{0};
}

template <typename T>
auto size(const optRef<T> &opt) {
  return opt ? opt->get().size() : typename T::size_type{0};
}

}  // namespace autopas::utils
