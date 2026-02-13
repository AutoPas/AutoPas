/**
 * @file PsVLC08Traversal.h
 * @date 20.12.2025
 * @author Lars Doll
 */

#pragma once

#include "PsVLTraversalInterface.h"
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler.h"

namespace autopas {

/**
 * This class provides the psvl_c08 traversal.
 *
 * The traversal uses the c08 base step performed on every single cell.
 * \image html C08.png "C08 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eight colors is applied.
 * \image html C08_domain.png "C08 domain coloring in 2D. 4 colors are required."
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor_T>
class PsVLC08Traversal : public C08BasedTraversal<ParticleCell_T, PairwiseFunctor_T>,
                         public PsVLTraversalInterface<ParticleCell_T> {
 public:
  /**
   * Constructor of the psvl_c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit PsVLC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor_T *pairwiseFunctor,
                            double interactionLength, const std::array<double, 3> &cellLength,
                            DataLayoutOption dataLayout, bool useNewton3)
      : C08BasedTraversal<ParticleCell_T, PairwiseFunctor_T>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                             dataLayout, useNewton3),

        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap,
                     dataLayout, useNewton3) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::psvl_c08; }

  /**
   * C08 traversals are always usable. Supports only AoS.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
    if (this->_dataLayout == DataLayoutOption::aos && this->_overlap[0] == 1 && this->_overlap[1] == 1 &&
        this->_overlap[2] == 1) {
      return true;
    }
    return false;
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

  /**
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) override;

 private:
  LCC08CellHandler<ParticleCell_T, PairwiseFunctor_T> _cellHandler;
};

template <class ParticleCell_T, class PairwiseFunctor_T>
void PsVLC08Traversal<ParticleCell_T, PairwiseFunctor_T>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) {
  _cellHandler.setOrientationList(list);
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void PsVLC08Traversal<ParticleCell_T, PairwiseFunctor_T>::traverseParticles() {
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(baseIndex);
  });
}

}  // namespace autopas
