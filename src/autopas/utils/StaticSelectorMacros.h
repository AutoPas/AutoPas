/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 02.05.19
 */

#pragma once

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
