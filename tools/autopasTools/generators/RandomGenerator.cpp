/**
 * @file RandomGenerator.cpp
 * @author seckler
 * @date 22.05.18
 */

#include "RandomGenerator.h"

double autopasTools::generators::RandomGenerator::fRand(double fMin, double fMax) {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

std::array<double, 3> autopasTools::generators::RandomGenerator::randomPosition(const std::array<double, 3> &boxMin,
                                                                                const std::array<double, 3> &boxMax) {
  std::array<double, 3> r{0, 0, 0};
  for (int d = 0; d < 3; ++d) {
    r[d] = fRand(boxMin[d], boxMax[d]);
  }
  return r;
}
