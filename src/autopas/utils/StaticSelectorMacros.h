/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 02.05.19
 */

#pragma once

/**
 * Macro to execute the code defined in ... with a static bool.
 * The bool will be renamed to be called c_BOOLNAME,
 * where BOOLNAME is the name of the bool passed in this macro.
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
