/**
 * @file Generators.h
 * @date 03.11.2026
 */
#pragma once

#include <pybind11/pybind11.h>

#include <vector>

#include "src/TypeDefinitions.h"

namespace autopasPython {

/**
 * Binds the submodule with the particle generators.
 * @param module
 */
void bindGenerators(pybind11::module_ &module);

/**
 * Converts particles to a dictionary of NumPy arrays, sorted by the ID of the particles.
 * @param particles
 * @return dictionary with the keys positions, velocities, forces, ids, and type_ids.
 */
pybind11::dict particlesToArrays(std::vector<ParticleType> particles);

}  // namespace autopasPython
