/**
 * @file AutoPasMacros.h
 * @author seckler
 * @date 14.05.19
 */

#pragma once

#ifdef __GNUC__
/**
 * @todo c++17: __attribute__((warn_unused_result)) should be replaced with [[nodiscard]]
 */
#define AUTOPAS_WARN_UNUSED_RESULT __attribute__((warn_unused_result))
#endif