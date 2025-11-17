/**
 *@file KokkosDataLayoutConverter.h
 *@author Luis Gall
 *@date 13.11.2025
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"

namespace autopas::utils {

class KokkosDataLayoutConverter {
 public:

  explicit KokkosDataLayoutConverter(DataLayoutOption dataLayout)
      : _dataLayout(dataLayout) {}

  template <class Input, class Output, typename Particle>
  void loadDataLayout(Input &srcParticles, Output &dstParticles, size_t numParticles) {

    // AoS to SoA
    for (size_t i = 0; i < numParticles; ++i) {
      (dstParticles. template set<0>(srcParticles(i).template get<Particle::AttributeNames::id>(), i));
    }
  }

  template <class Input, class Output, typename Particle>
  void storeDataLayout(Input &srcParticles, Output &dstParticles, size_t numParticles) {

    // SoA to AoS
    for (size_t i = 0; i < numParticles; ++i) {
      (dstParticles(i). template set<Particle::AttributeNames::id>(srcParticles. template get<0>(i)));
    }
  }

 private:

  DataLayoutOption _dataLayout;
};

}  // namespace autopas::utils
