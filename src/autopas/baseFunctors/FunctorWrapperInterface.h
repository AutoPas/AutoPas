//
// Created by markus on 9/25/25.
//

#pragma once
#include "autopas/utils/SoA.h"

template <typename Particle_T>
class FunctorWrapperInterface {
public:
  virtual ~FunctorWrapperInterface() = default;

  virtual void initTraversal() = 0;

  virtual void endTraversal(bool newton3) = 0;

  virtual void SoALoader(void *cell, void *soa, size_t offset, bool skipSoAResize) = 0;
  virtual void SoAExtractor(void *cell, void *soa, size_t offset) = 0;

  virtual bool allowsNewton3() = 0;
  virtual bool allowsNonNewton3() = 0;

  virtual bool isRelevantForTuning() = 0;

  virtual std::string getName() = 0;

  [[nodiscard]] virtual double getCutoff() const = 0;

  [[nodiscard]] virtual size_t getNumFlops() const = 0;

  [[nodiscard]] virtual double getHitRate() const = 0;
};