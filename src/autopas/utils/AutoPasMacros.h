/**
 * @file AutoPasMacros.h
 * @author seckler
 * @date 14.05.19
 */

#pragma once

// We try not to use precompiler macros so if you feel like adding something here
// think if it really can not be done in a different way!

/**
 * @def AUTOPAS_ALWAYS_INLINE
 * Hints/forces the compiler to inline performance-critical functions in optimized builds.
 */
#ifndef AUTOPAS_ALWAYS_INLINE
#if defined(_MSC_VER)
#define AUTOPAS_ALWAYS_INLINE __forceinline
#elif defined(__OPTIMIZE__)
#if defined(__has_attribute)
#if __has_attribute(always_inline)
#define AUTOPAS_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define AUTOPAS_ALWAYS_INLINE inline
#endif
#elif defined(__GNUC__) || defined(__clang__) || defined(__FUJITSU) || defined(__CLANG_FUJITSU)
#define AUTOPAS_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define AUTOPAS_ALWAYS_INLINE inline
#endif
#else
#define AUTOPAS_ALWAYS_INLINE inline
#endif
#endif