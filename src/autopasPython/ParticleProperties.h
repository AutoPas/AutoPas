/**
 * @file ParticleProperties.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/pybind11.h>

namespace autopasPython {

/**
 * Binds the ParticlePropertiesLibrary.
 * @param module
 */
void bindParticleProperties(pybind11::module_ &module);

}  // namespace autopasPython
