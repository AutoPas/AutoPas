/**
 * @file StaticBoolSelector.h
 * @author seckler
 * @date 02.05.19
 */

#pragma once

namespace autopas::utils {
/**
 * Function to execute the code passed in the lambda with a static bool.
 * The static value of the boolean will be passed using the first argument of function.
 * @param theBool The bool that should be used statically.
 * @param func Function to be called, should accept a static bool type, e.g., [&](auto theBool){};
 * @return Returns whatever func returns.
 */
template <typename F>
decltype(auto) withStaticBool(bool theBool, F &&func) {
  if (theBool) {
    std::true_type t;
    return func(t);
  } else {
    std::false_type f;
    return func(f);
  }
}
}  // namespace autopas::utils