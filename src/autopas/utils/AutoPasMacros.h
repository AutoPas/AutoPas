/**
 * @file AutoPasMacros.h
 * @author seckler
 * @date 14.05.19
 */

#pragma once

/**
 * Highly encourages the compiler to produce a warning if the return value is ignored.
 * @note: this is c++17
 */
#define AUTOPAS_WARN_UNUSED_RESULT [[nodiscard]]
