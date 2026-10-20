/**
 * @file NumericConcept.h
 * @author Samuel J. Newcome
 * @date 30.04.2026
 */

#pragma once

#include <concepts>

/**
 * Concept to check if any type is arithmetic.
 */
template <typename T>
concept Numeric = std::is_arithmetic_v<T>;
