/**
 * @file SortingThreshold.h
 * @author Henry Meyran
 * @date 07.07.2026
 */

#pragma once

#include <cstddef>

namespace autopas {
/**
 * Default number of particles in two cells from which point on AoS sorting should be performed by CellFunctor(3B).
 * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
 */
constexpr size_t defaultAoSSortingThreshold = 8;

/**
 * Default sum of the SoA buffer sizes of two cells from which point on SoA sorting should be performed.
 */
constexpr size_t defaultSoASortingThreshold = 25;
}  // namespace autopas
