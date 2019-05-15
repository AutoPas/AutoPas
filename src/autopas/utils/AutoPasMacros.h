/**
 * @file AutoPasMacros.h
 * @author seckler
 * @date 14.05.19
 */

#pragma once

#if defined(__GNUC__) or defined(__clang__)
/**
 * For gcc + clang: warn if return value is not used.
 * @todo c++17: __attribute__((warn_unused_result)) should be replaced with [[nodiscard]]
 */
#define AUTOPAS_WARN_UNUSED_RESULT __attribute__((warn_unused_result))
#else
/**
 * empty if not gnu
 */
#define AUTOPAS_WARN_UNUSED_RESULT
#endif