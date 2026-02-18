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
}  // namespace autopas::utils
