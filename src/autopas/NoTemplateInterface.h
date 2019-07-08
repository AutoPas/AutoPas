/**
 * @file NoTemplateInterface.h
 * An interface to avoid templates for TuningStrategies
 * @author Candas
 * @date 07.07.2019
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include "autopas/utils/ArrayMath.h"

namespace autopas {
class NoTemplateInterface {
 public:
  ~NoTemplateInterface() = default;

  std::array<double, 3> getBoxMin() { return {0, 0, 0}; }
  std::array<double, 3> getBoxMax() { return {0, 0, 0}; }

  unsigned long getNumberOfParticles() { return 0; }

  double getCutoff() { return 0; }

  double getVerletSkin() { return 0; }

  unsigned int getVerletRebuildFrequency() { return 0; }
};
}