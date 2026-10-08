/**
 * @file PeriodicBoundariesHelpers.h
 * @author S. Newcome
 * @date 08.10.2026
 */

#pragma once

#include <array>
#include <vector>

/**
 * Wraps a position that left the box back into it, as if all boundaries were periodic.
 * @param pos
 * @param boxMin
 * @param boxMax
 * @return position inside [boxMin, boxMax).
 */
inline std::array<double, 3> wrapIntoBox(std::array<double, 3> pos, const std::array<double, 3> &boxMin,
                                         const std::array<double, 3> &boxMax) {
  for (size_t d = 0; d < 3; ++d) {
    const double boxLength = boxMax[d] - boxMin[d];
    if (pos[d] < boxMin[d]) {
      pos[d] += boxLength;
    } else if (pos[d] >= boxMax[d]) {
      pos[d] = pos[d] - boxLength;
    }
  }
  return pos;
}

/**
 * Generates periodic images of a particle within the interaction length of the boundaries.
 * Images keep the id of the particle they are copied from.
 * @tparam Particle_T
 * @param particle Particle inside the box.
 * @param boxMin
 * @param boxMax
 * @param interactionLength
 * @return The images, which all lie outside the box.
 */
template <class Particle_T>
std::vector<Particle_T> generatePeriodicImages(const Particle_T &particle, const std::array<double, 3> &boxMin,
                                               const std::array<double, 3> &boxMax, const double interactionLength) {
  std::vector<Particle_T> images;
  const auto &pos = particle.getR();
  // Check every face, edge, and corner of the box, given as a direction from the box center. The particle has an
  // image on the opposite side if it is within the interaction length of all boundaries of this direction.
  for (const int x : {-1, 0, 1}) {
    for (const int y : {-1, 0, 1}) {
      for (const int z : {-1, 0, 1}) {
        if (x == 0 and y == 0 and z == 0) {
          continue;
        }
        const std::array<int, 3> direction{x, y, z};
        bool nearBoundaries = true;
        std::array<double, 3> shift{};
        for (size_t d = 0; d < 3; ++d) {
          const double boxLength = boxMax[d] - boxMin[d];
          if (direction[d] == -1) {
            nearBoundaries = nearBoundaries and pos[d] < boxMin[d] + interactionLength;
            shift[d] = boxLength;
          } else if (direction[d] == 1) {
            nearBoundaries = nearBoundaries and pos[d] >= boxMax[d] - interactionLength;
            shift[d] = -boxLength;
          }
        }
        if (nearBoundaries) {
          auto image = particle;
          image.addR(shift);
          images.push_back(image);
        }
      }
    }
  }
  return images;
}