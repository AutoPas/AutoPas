/**
 * @file TwoCellsInteractionRatioGenerator.h
 * @author H. Meyran
 */

#pragma once

#include <array>
#include <cmath>
#include <random>

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopasTools::generators::TwoCellsInteractionRatioGenerator {

/**
 * Fills two cells adjacent in x with n particles each such that exactly
 * round(n * ratio) particles per cell have a guaranteed in-range partner
 * in the other cell (3D distance < cutoff), and the remaining particles
 * are placed far enough that no cross-cell interaction is possible.
 *
 * Preconditions (checked via ExceptionHandler):
 *   - boxMax1[0] == boxMin2[0]                      — shared boundary
 *   - (boxMax1[0] - boxMin1[0]) > 1.5 * cutoff      — cell1 x-extent sufficient for far zone
 *   - (boxMax2[0] - boxMin2[0]) > 1.5 * cutoff      — cell2 x-extent sufficient for far zone
 *   - ratio in [0.0, 1.0]
 *
 * Particle IDs are assigned sequentially starting from defaultParticle.getID():
 *   cell1 particles: [id0,     id0 + n)
 *   cell2 particles: [id0 + n, id0 + 2n)
 *
 * Interaction guarantees:
 *   - The first k = round(n * ratio) particles of cell1 are each within cutoff of
 *     the corresponding particle (same index) in cell2 (3D distance = 2·dx < cutoff).
 *   - The remaining n-k far particles have an x-distance >= cutoff to every particle
 *     in the other cell, so they produce no cross-cell interactions.
 *
 * @tparam Container Any container supporting addParticle().
 * @tparam Particle Particle type.
 * @param cell1 First cell (x < boundary).
 * @param cell2 Second cell (x >= boundary).
 * @param boxMin1 Minimum coordinates of cell1.
 * @param boxMax1 Maximum coordinates of cell1. boxMax1[0] must equal boxMin2[0].
 * @param boxMin2 Minimum coordinates of cell2. boxMin2[0] must equal boxMax1[0].
 * @param boxMax2 Maximum coordinates of cell2.
 * @param n Number of particles per cell.
 * @param ratio Fraction of particles per cell with a guaranteed in-range cross-cell partner. Must be in [0, 1].
 * @param cutoff Interaction cutoff radius.
 * @param defaultParticle Template particle; IDs are assigned sequentially starting from defaultParticle.getID().
 * @param seed Random seed for reproducibility.
 */
template <class Container, class Particle>
void fillWithParticles(Container &cell1, Container &cell2, const std::array<double, 3> &boxMin1,
                       const std::array<double, 3> &boxMax1, const std::array<double, 3> &boxMin2,
                       const std::array<double, 3> &boxMax2, std::size_t n, double ratio, double cutoff,
                       const Particle &defaultParticle = Particle{}, unsigned int seed = 42) {
  if (std::abs(boxMax1[0] - boxMin2[0]) > 1e-10) {
    autopas::utils::ExceptionHandler::exception(
        "TwoCellsInteractionRatioGenerator: cells must share an x-boundary "
        "(boxMax1[0]={} != boxMin2[0]={})",
        boxMax1[0], boxMin2[0]);
  }
  if ((boxMax1[0] - boxMin1[0]) <= 1.5 * cutoff) {
    autopas::utils::ExceptionHandler::exception(
        "TwoCellsInteractionRatioGenerator: cell1 x-extent ({}) must be > 1.5 * cutoff ({})",
        boxMax1[0] - boxMin1[0], 1.5 * cutoff);
  }
  if ((boxMax2[0] - boxMin2[0]) <= 1.5 * cutoff) {
    autopas::utils::ExceptionHandler::exception(
        "TwoCellsInteractionRatioGenerator: cell2 x-extent ({}) must be > 1.5 * cutoff ({})",
        boxMax2[0] - boxMin2[0], 1.5 * cutoff);
  }
  if (ratio < 0.0 || ratio > 1.0) {
    autopas::utils::ExceptionHandler::exception(
        "TwoCellsInteractionRatioGenerator: ratio ({}) must be in [0, 1]", ratio);
  }

  const double boundary = boxMax1[0];
  const auto k = static_cast<std::size_t>(std::round(static_cast<double>(n) * ratio));
  const std::size_t farCount = n - k;
  const auto idBase = static_cast<std::size_t>(defaultParticle.getID());

  std::mt19937 rng(seed);

  // Small offset so paired particles are strictly inside the cutoff sphere and their cell bounds.
  const double eps = cutoff * 1e-6;
  std::uniform_real_distribution<double> dxDist(eps, cutoff / 2.0 - eps);

  auto uniform = [&](double lo, double hi) { return std::uniform_real_distribution<double>(lo, hi)(rng); };

  // k paired (in-range) particles: 3D distance between partners = 2·dx < cutoff
  for (std::size_t i = 0; i < k; ++i) {
    const double dx = dxDist(rng);
    const double y = uniform(boxMin1[1], boxMax1[1]);
    const double z = uniform(boxMin1[2], boxMax1[2]);

    Particle p1(defaultParticle);
    p1.setR({boundary - dx, y, z});
    p1.setID(idBase + i);
    p1.setOwnershipState(autopas::OwnershipState::owned);
    cell1.addParticle(p1);

    Particle p2(defaultParticle);
    p2.setR({boundary + dx, y, z});
    p2.setID(idBase + n + i);
    p2.setOwnershipState(autopas::OwnershipState::owned);
    cell2.addParticle(p2);
  }

  // n-k far particles in cell1: x-distance to any cell2 particle >= cutoff
  for (std::size_t i = 0; i < farCount; ++i) {
    Particle p(defaultParticle);
    p.setR({uniform(boxMin1[0], boundary - cutoff), uniform(boxMin1[1], boxMax1[1]),
            uniform(boxMin1[2], boxMax1[2])});
    p.setID(idBase + k + i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    cell1.addParticle(p);
  }

  // n-k far particles in cell2: x-distance to any cell1 particle >= cutoff
  for (std::size_t i = 0; i < farCount; ++i) {
    Particle p(defaultParticle);
    p.setR({uniform(boundary + cutoff, boxMax2[0]), uniform(boxMin2[1], boxMax2[1]),
            uniform(boxMin2[2], boxMax2[2])});
    p.setID(idBase + n + k + i);
    p.setOwnershipState(autopas::OwnershipState::owned);
    cell2.addParticle(p);
  }
}

}  // namespace autopasTools::generators::TwoCellsInteractionRatioGenerator
