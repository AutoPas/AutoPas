/**
 * @file SPHParticle.cpp
 * @author seckler
 * @date 19.01.18
 */

#include "autopas/sph/SPHParticle.h"

#include <cmath>

#include "autopas/utils/ArrayMath.h"

using namespace autopas::sph;

void SPHParticle::addAcceleration(const std::array<double, 3> &acc) {
  using namespace autopas::utils::ArrayMath::literals;
  _acc = _acc + acc;
}

void SPHParticle::subAcceleration(const std::array<double, 3> &acc) {
  using namespace autopas::utils::ArrayMath::literals;
  _acc = _acc - acc;
}

void SPHParticle::calcPressure() {
  const double hcr = 1.4;
  _pressure = (hcr - 1.0) * _density * _energy;
  _snds = sqrt(hcr * _pressure / _density);
}
