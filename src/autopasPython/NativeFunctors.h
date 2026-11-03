/**
 * @file NativeFunctors.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/pybind11.h>

namespace autopasPython {

/**
 * Binds the functor base class and the functors implemented in C++.
 * @param module
 */
void bindNativeFunctors(pybind11::module_ &module);

}  // namespace autopasPython
