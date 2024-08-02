/**
 * @file UniformGenerator.cpp
 * @author seckler
 * @date 22.05.18
 */

#include <random>

#include "UniformGenerator.h"

namespace autopasTools::generators {

std::array<double, 3> UniformGenerator::randomPosition(std::mt19937 &generator, const std::array<double, 3> &boxMin,
                                                      const std::array<double, 3> &boxMax) {
  std::array<std::uniform_real_distribution<double>, 3> distributions = {
      std::uniform_real_distribution<double>{boxMin[0], boxMax[0]},
      std::uniform_real_distribution<double>{boxMin[1], boxMax[1]},
      std::uniform_real_distribution<double>{boxMin[2], boxMax[2]},
  };
  std::array<double, 3> r{};
  for (int d = 0; d < 3; ++d) {
    r[d] = distributions[d](generator);
  }
  return r;
}
}  // namespace autopasTools::generators
