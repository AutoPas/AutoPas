/**
 * @file HGColorTraversal.h
 * @author atacann
 * @date 02.01.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGColorTraversal : public HGTraversalBase<ParticleCell>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGColorTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell>(dataLayout, useNewton3), _functor(functor){};

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    _functor->setCutoff(this->_cutoffs[level]);
    return std::make_unique<LCC08Traversal<ParticleCell, Functor>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  };

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Currently only AoS and newton3 is supported on hgrid_color traversal.");
    }
    // computeInteractions for each level independently first
    // TODO: For SoA, do not call endTraversal() here (do not load SoA to Particle), do them all together after cross level interactions
    // TODO: If cross level SoA is implemented only, otherwise not needed
    for (size_t level = 0; level < this->_numLevels; level++) {
      auto traversal = generateNewTraversal(level);
      this->_levels->at(level)->computeInteractions(traversal.get());
    }
    // computeInteractions across different levels
    for (size_t upperLevel = 1; upperLevel < this->_numLevels; upperLevel++) {
      // only look top-down if newton3 is enabled, both ways otherwise
      size_t levelLimit = this->_useNewton3 ? upperLevel : this->_numLevels;
      for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
        if (lowerLevel == upperLevel) {
          continue;
        }
        // calculate cutoff for level pair
        const double cutoff = (this->_cutoffs[upperLevel] + this->_cutoffs[lowerLevel]) / 2;
        const std::array<double, 3> upperLength = this->_levels->at(upperLevel)->getTraversalSelectorInfo().cellLength;
        //const std::array<double, 3> lowerLength = this->_levels->at(lowerLevel)->getTraversalSelectorInfo().cellLength;
        _functor->setCutoff(cutoff);
        const double interactionLength = cutoff + this->_skin;

        // find out the stride so that cells we check on lowerLevel do not intersect
        std::array<unsigned long, 3> stride{};
        for (size_t i = 0; i < 3; i++) {
          //stride[i] = 1 + std::ceil(std::ceil(interactionLength / lowerLength[i]) * 2 * lowerLength[i] / upperLength[i]);
          if (this->_useNewton3) {
            stride[i] = 1 + static_cast<unsigned long>(std::ceil(interactionLength  * 2 / upperLength[i]));
          }
          else {
            // do c01 traversal if newton3 is disabled
            stride[i] = 1;
          }
        }
        using utils::ArrayUtils::operator<<;
        // std::ostringstream text;
        // text << "Stride: " << stride << " lowerLength: " << lowerLength << " upperLength: " << upperLength <<
        //   " \nlowerLevel " << lowerLevel << " upperLevel " << upperLevel << " interactionLength " << interactionLength;
        // AutoPasLog(INFO, text.str());

        const std::array<double, 3> dir = {interactionLength, interactionLength, interactionLength};

        // TODO: SoAView?
        const auto traverse = [&lowerLevel = *(this->_levels->at(lowerLevel)), &dir, &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock(),
          &functor = this->_functor, &useNewton3 = this->_useNewton3](unsigned long x, unsigned long y, unsigned long z) {
          auto &cell = upperLevelCB.getCell({x, y, z});
          for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
            const std::array<double, 3> &pos = p1Ptr->getR();
            const auto posMin = pos - dir, posMax = pos + dir;
            lowerLevel.forEachInRegion([&p1Ptr, &functor, &useNewton3](Particle &j) { functor->AoSFunctor(*p1Ptr, j, useNewton3); },
                                  posMin, posMax, IteratorBehavior::ownedOrHalo);
          }
        };

        const auto end = this->getTraversalSelectorInfo(upperLevel).cellsPerDim;
        const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

        // do the colored traversal
        const unsigned long numColors = stride[0] * stride[1] * stride[2];
        for (unsigned long col = 0; col < numColors; ++col) {

          const std::array<unsigned long, 3> start(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

          const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
          const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];

          AUTOPAS_OPENMP(parallel for schedule(dynamic, 1) collapse(3))
          for (unsigned long x = start_x; x < end_x; x += stride_x) {
            for (unsigned long y = start_y; y < end_y; y += stride_y) {
              for (unsigned long z = start_z; z < end_z; z += stride_z) {
                traverse(x, y, z);
              }
            }
          }
        }
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_color; };

  [[nodiscard]] bool isApplicable() const override {
    return true;
  };

  void initTraversal() override {
    // store the initial cutoff of functor so that it does not change during traversal (creates issues with some tests
    // otherwise)
    _storedCutoff = _functor->getCutoff();
  };

  void endTraversal() override {
    // restore the initial cutoff of functor
    _functor->setCutoff(_storedCutoff);
  };

 protected:
  Functor *_functor;
  double _storedCutoff = 0;
};
}  // namespace autopas
