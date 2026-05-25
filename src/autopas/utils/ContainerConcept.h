/**
 * @file ContainerConcept.h
 * @author J. Schuhmacher
 * @date 11.02.26
 */

#pragma once

#include <array>
#include <concepts>
#include <set>
#include <type_traits>
#include <vector>

namespace autopas::utils {

/**
 * Collection of structs that define what we consider a container.
 */
namespace is_container_impl {
/**
 * Default case: T is not a container.
 * @tparam T
 */
template <typename T>
struct is_container : std::false_type {};
/**
 * Specialization to allow std::array.
 * @tparam T
 * @tparam N
 */
template <typename T, std::size_t N>
struct is_container<std::array<T, N>> : std::true_type {};

/**
 * Specialization to allow std::vector.
 * @tparam Args
 */
template <typename... Args>
struct is_container<std::vector<Args...>> : std::true_type {};

/**
 * Specialization to allow std::set.
 * @tparam Args
 */
template <typename... Args>
struct is_container<std::set<Args...>> : std::true_type {};

}  // namespace is_container_impl

/**
 * Concept to check if something is a container (vector, array or set).
 * @tparam T Type to check.
 */
template <typename T>
concept ContainerType = is_container_impl::is_container<std::decay_t<T>>::value;

}  // namespace autopas::utils