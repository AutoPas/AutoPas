/**
 * @file PythonFunctors.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/pybind11.h>

namespace autopasPython {

/**
 * Binds the base classes for functors implemented in Python, as well as the Particle they operate on.
 * @param module
 */
void bindPythonFunctors(pybind11::module_ &module);

}  // namespace autopasPython
