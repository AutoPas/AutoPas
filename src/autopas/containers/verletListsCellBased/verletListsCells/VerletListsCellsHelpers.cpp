/**
 * @file VerletListsCellsHelpers.cpp
 * @author F.Gratl
 * @date 29.05.24
 */

#include "VerletListsCellsHelpers.h"

#include "cmath"

size_t autopas::VerletListsCellsHelpers::estimateListLength(size_t numParticles, const std::array<double, 3> &boxSize,
                                                            double interactionLength, double correctionFactor) {
  const auto boxVolume = boxSize[0] * boxSize[1] * boxSize[2];
  constexpr double sphereConstant = 4. / 3. * M_PI;
  const auto listVolume = sphereConstant * interactionLength * interactionLength * interactionLength;
  const auto volumeFraction = listVolume / boxVolume;
  // ceil because if we need space for e.g. 3.5 particles better reserve for 4.
  return static_cast<size_t>(std::ceil(numParticles * volumeFraction * correctionFactor));
}
