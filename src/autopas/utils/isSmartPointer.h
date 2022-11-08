/**
 * @file isSmartPointer.h
 * @author F. Gratl
 * @date 08.11.22
 */

#pragma once

#include "memory"

namespace autopas::utils {

/**
 * Type trait to determine if a type is a shared pointer at compile time.
 * @tparam T
 */
template <typename T>
struct is_shared_ptr : std::false_type {};
/**
 * True specialization for shared_ptr
 * @tparam T
 */
template <typename T>
struct is_shared_ptr<std::shared_ptr<T>> : std::true_type {};

/**
 * Type trait to determine if a type is a unique pointer at compile time.
 * @tparam T
 */
template <typename T>
struct is_unique_ptr : std::false_type {};
/**
 * True specialization for unique_ptr
 * @tparam T
 */
template <typename T>
struct is_unique_ptr<std::unique_ptr<T>> : std::true_type {};

/**
 * Type trait to determine if a type is a smart pointer at compile time.
 * @tparam T
 */
template <typename T>
struct is_smart_ptr : std::false_type {};
/**
 * True specialization for unique_ptr
 * @tparam T
 */
template <typename T>
struct is_smart_ptr<std::unique_ptr<T>> : std::true_type {};
/**
 * True specialization for shared_ptr
 * @tparam T
 */
template <typename T>
struct is_smart_ptr<std::shared_ptr<T>> : std::true_type {};
/**
 * True specialization for weak_ptr
 * @tparam T
 */
template <typename T>
struct is_smart_ptr<std::weak_ptr<T>> : std::true_type {};
}  // namespace autopas::utils