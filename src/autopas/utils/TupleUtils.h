/**
 * @file TupleUtils.h
 * @author F. Gratl
 * @date 26/05/2020
 */

#pragma once

#include <tuple>

namespace autopas::utils::TupleUtils {

/**
 * Applies a function f on every element of the tuple. Elements are processed in the order they are declared.
 * @note Needs not be std::tuple but anything that satisfies std::apply.
 * @tparam T Type of the tuple
 * @tparam F Type of the function
 * @param tuple The tuple to iterate
 * @param f The function to apply. Typically [&](auto &elem) {...}
 */
template <class T, class F>
void for_each(T &&tuple, F &&f) {
  return std::apply(
      // function that takes all elements of the tuple as variadic argument
      [&](auto &... t) {
        // unpack the tuple using a fold (... op pack) with op as the comma operator to apply a lambda on every element
        // This will cause problems if anyone overloads the comma operator for elements of tuple!
        (..., f(t));
      },
      tuple);
}

}  // namespace autopas::utils::TupleUtils