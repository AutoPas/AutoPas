/**
 * @file SPHParticle.cpp
 * @author seckler
 * @date 19.01.18
 */

#include "autopas/sph/SPHParticle.h"

#include <cmath>

#include "autopas/utils/ArrayMath.h"

using namespace autopas::sph;

KOKKOS_FUNCTION
void SPHParticle::addAcceleration(const std::array<double, 3> &acc) { _acc = utils::ArrayMath::add(_acc, acc); }

KOKKOS_FUNCTION
void SPHParticle::subAcceleration(const std::array<double, 3> &acc) { _acc = utils::ArrayMath::sub(_acc, acc); }

KOKKOS_FUNCTION
void SPHParticle::calcPressure() {
  const double hcr = 1.4;
  _pressure = (hcr - 1.0) * _density * _energy;
  _snds = sqrt(hcr * _pressure / _density);
}
