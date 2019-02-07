/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 07.02.19
 */

#pragma once

/**
 * Will execute the passed body (=variadic argument) with the static value of the bool.
 * @param theBool The bool to be used.
 * @note The second Argument is variadic such that commas pose no problem.
 */
#define AUTOPAS_WITH_STATIC_BOOL(theBool, ...) \
  {                                            \
    if (theBool) {                             \
      constexpr bool c_##theBool = true;       \
      __VA_ARGS__                              \
    } else {                                   \
      constexpr bool c_##theBool = false;      \
      __VA_ARGS__                              \
    }                                          \
  }
