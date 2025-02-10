/**
 * @file HGC01Traversal.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC01Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGC01Traversal : public HGTraversalBase<ParticleCell>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGC01Traversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell>(dataLayout, useNewton3), _functor(functor){};

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    //_functor->setCutoff(this->_cutoffs[level]);
    return std::make_unique<LCC01Traversal<ParticleCell, Functor>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  };

  void traverseParticles() override {
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Currently only AoS and No Newton3 is supported on hgrid_c01 traversal.");
    }
    // computeInteractions for each level independently first
    for (size_t level = 0; level < this->_numLevels; level++) {
      auto traversal = generateNewTraversal(level);
      this->_levels->at(level)->computeInteractions(traversal.get());
    }
    // computeInteractions across different levels
    for (size_t level = 0; level < this->_numLevels; level++) {
      for (size_t innerLevel = 0; innerLevel < this->_numLevels; innerLevel++) {
        if (level == innerLevel) {
          continue;
        }
        // calculate cutoff for level pair
        const double cutoff = (this->_cutoffs[level] + this->_cutoffs[innerLevel]) / 2;
        //_functor->setCutoff(cutoff);
        const double interactionLength = cutoff + this->_skin;
        AUTOPAS_OPENMP(parallel) {
          using utils::ArrayMath::operator-;
          using utils::ArrayMath::operator+;
          auto it = this->_levels->at(level)->begin(IteratorBehavior::owned);
          for (; it.isValid(); ++it) {
            const std::array<double, 3> &pos = it->getR();
            const std::array<double, 3> dir = {cutoff, cutoff, cutoff};
            const auto posMin = pos - dir, posMax = pos + dir;
            this->_levels->at(innerLevel)
                ->forEachInRegion([&it, &functor = this->_functor, &useNewton3 = this->_useNewton3](
                                      Particle &j) { functor->AoSFunctor(*it, j, useNewton3); },
                                  posMin, posMax, IteratorBehavior::ownedOrHalo);
          }
        }
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_c01_iterator; };

  [[nodiscard]] bool isApplicable() const override {
    return this->_useNewton3 == false && this->_dataLayout == DataLayoutOption::aos;
  };

  void initTraversal() override { _storedCutoff = _functor->getCutoff(); };

  void endTraversal() override { _functor->setCutoff(_storedCutoff); };

 protected:
  Functor *_functor;
  double _storedCutoff = 0;
};
}  // namespace autopas
