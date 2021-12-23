/**
 * @file BoundaryConditions.h
 * @author F. Gratl
 * @date 2/8/19
 */

#pragma once
#include <array>
#include <vector>

#include "autopas/AutoPasDecl.h"

/**
 * This namespace implements functions for periodic boundaries using AutoPas.
 */
namespace BoundaryConditions {
/**
 * Anonymous namespace for internal helpers.
 */
namespace {
/**
 * Convert the leaving particle to entering particles.
 * This function updates the positions of the given particles to periodic images inside the domain.
 * @tparam Particle
 * @param autoPas
 * @param particles Particles that are outside the domain box.
 */
template <class Particle>
void wrapPositionsAroundBoundaries(autopas::AutoPas<Particle> &autoPas, std::vector<Particle> &particles) {
  const auto &boxMin = autoPas.getBoxMin();
  const auto &boxMax = autoPas.getBoxMax();
  const auto boxLength = std::array<double, 3>{
      boxMax[0] - boxMin[0],
      boxMax[1] - boxMin[1],
      boxMax[2] - boxMin[2],
  };

  for (auto &p : particles) {
    auto pos = p.getR();
    for (auto dim = 0ul; dim < pos.size(); dim++) {
      if (pos[dim] < boxMin[dim]) {
        // has to be smaller than boxMax
        pos[dim] = std::min(std::nextafter(boxMax[dim], boxMin[dim]), pos[dim] + boxLength[dim]);
      } else if (pos[dim] >= boxMax[dim]) {
        // should at least be boxMin
        pos[dim] = std::max(boxMin[dim], pos[dim] - boxLength[dim]);
      }
    }
    p.setR(pos);
  }
}

/**
 * Identifies particles that are near the edge of the domain and returns properly shifted periodic copies of them as
 * halo particles.
 *
 * @tparam Particle
 * @param autoPas
 * @return vector of particles that are already shifted for the next process.
 */
template <class Particle>
std::vector<Particle> identifyNewHaloParticles(autopas::AutoPas<Particle> &autoPas) {
  std::vector<Particle> haloParticles;

  const auto &boxMin = autoPas.getBoxMin();
  const auto &boxMax = autoPas.getBoxMax();
  const auto boxLength = std::array<double, 3>{
      boxMax[0] - boxMin[0],
      boxMax[1] - boxMin[1],
      boxMax[2] - boxMin[2],
  };

  // look in every direction
  for (auto x : {-1, 0, 1}) {
    for (auto y : {-1, 0, 1}) {
      for (auto z : {-1, 0, 1}) {
        if (x == 0 and y == 0 and z == 0) continue;
        std::array<int, 3> direction{x, y, z};
        // Define a region where to search for particles that need to be copied to the halo.
        std::array<double, 3> min{}, max{}, shiftVec{};
        for (auto dim = 0ul; dim < direction.size(); ++dim) {
          // Flag for the periodic wrap. Will occur in at least one direction.
          bool needsShift = false;
          // The search domain has to be enlarged by the skin as the position of the particles is not certain.
          if (direction[dim] == -1) {
            // just inside of boxMin
            min[dim] = boxMin[dim] - autoPas.getVerletSkin();
            max[dim] = boxMin[dim] + autoPas.getCutoff() + autoPas.getVerletSkin();
            needsShift = true;
          } else if (direction[dim] == 1) {
            // just inside of boxMax
            min[dim] = boxMax[dim] - autoPas.getCutoff() - autoPas.getVerletSkin();
            max[dim] = boxMax[dim] + autoPas.getVerletSkin();
            needsShift = true;
          } else {  // direction[dim] == 0
            min[dim] = boxMin[dim] - autoPas.getVerletSkin();
            max[dim] = boxMax[dim] + autoPas.getVerletSkin();
          }
          if (needsShift) {
            // particles left of min (direction == -1) need a shift by +boxLength. For right of max the opposite.
            shiftVec[dim] = -boxLength[dim] * direction[dim];
          } else {
            shiftVec[dim] = 0;
          }
        }
        // Find non-halo particles in the designated region, create periodic copies of them and insert those as halos.
        autoPas.forEachInRegion([&] (Particle &p) {
                    auto particleCopy = p;
                    particleCopy.addR(shiftVec);
                    haloParticles.push_back(particleCopy);
        }, autopas::IteratorBehavior::owned);
      }
    }
  }
  return haloParticles;
}

/**
 * Insertes the given entering particles into the AutoPas object.
 * @tparam Particle
 * @param autoPas
 * @param enteringParticles
 */
template <class Particle>
void addEnteringParticles(autopas::AutoPas<Particle> &autoPas, std::vector<Particle> &enteringParticles) {
  for (auto &p : enteringParticles) {
    autoPas.addParticle(p);
  }
}

/**
 * Inserts the given halo particles into the AutoPas object.
 * @tparam Particle
 * @param autoPas
 * @param haloParticles
 */
template <class Particle>
void addHaloParticles(autopas::AutoPas<Particle> &autoPas, std::vector<Particle> &haloParticles) {
  for (auto &p : haloParticles) {
    autoPas.addHaloParticle(p);
  }
}

}  // namespace

/**
 * Applies periodic boundary conditions to inner and halo particles.
 *
 * Inner particles that leave the domain are re-inserted on the opposing side.
 * For inner particles that are near the edge of the domain, copies are placed in the halo on the opposing side.
 *
 * @tparam Particle
 * @param autoPas
 * @param forceUpdate If a container update should be forced or left to AutoPas.
 */
template <class Particle>
void applyPeriodic(autopas::AutoPas<Particle> &autoPas, bool forceUpdate) {
  // 1. update Container; return value is a vector of the particles leaving the domain box.

  // auto leavingParticles = autoPas.updateContainer();

  // 2. apply periodic wrap by shifting positions of leaving particles to positions of periodic images.
  // wrapPositionsAroundBoundaries(autoPas, leavingParticles);
  // 2b. re-insert shifted particles
  // addEnteringParticles(autoPas, leavingParticles);

  // 3. identify inner particles for which a periodic copy in the opposing halo region is needed.
  auto haloParticles = identifyNewHaloParticles(autoPas);
  addHaloParticles(autoPas, haloParticles);
}
}  // namespace BoundaryConditions
