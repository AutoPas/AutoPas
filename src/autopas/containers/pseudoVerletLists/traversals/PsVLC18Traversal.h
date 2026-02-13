/**
 * @file PsVLC18Traversal.h
 * @date 06.12.2025
 * @author Lars Doll
 */

#pragma once

#include "PsVLTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/cellTraversals/C18BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the psVl_c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell.
 * \image html C18.png "C18 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eighteen colors is applied.
 * \image html C18_domain.png "C18 domain coloring in 2D. 6 colors are required."
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor_T>
class PsVLC18Traversal : public C18BasedTraversal<ParticleCell_T, PairwiseFunctor_T>,
                         public PsVLTraversalInterface<ParticleCell_T> {
 public:
  /**
   * Constructor of the psvl_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit PsVLC18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor_T *pairwiseFunctor,
                            const double interactionLength, const std::array<double, 3> &cellLength,
                            DataLayoutOption dataLayout, bool useNewton3)
      : C18BasedTraversal<ParticleCell_T, PairwiseFunctor_T>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                             dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength, dataLayout, useNewton3) {
    this->computeOffsets();
  }

  void traverseParticles() override;

  /**
   * Computes all interactions between the base cell and adjacent cells with greater ID.
   * @param x X-index of the base cell.
   * @param y Y-index of the base cell.
   * @param z Z-index of the base cell.
   */
  void processBaseCell(unsigned long x, unsigned long y, unsigned long z);

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::psvl_c18; }

  /**
   * C18 traversal is always usable. Supports only AoS.
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
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) override;

  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell_T, PairwiseFunctor_T,
                        /*bidirectional*/ true>
      _cellFunctor;
};

template <class ParticleCell_T, class PairwiseFunctor_T>
void PsVLC18Traversal<ParticleCell_T, PairwiseFunctor_T>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) {
  _cellFunctor.setOrientationList(list);
}

template <class ParticleCell_T, class PairwiseFunctor_T>
void PsVLC18Traversal<ParticleCell_T, PairwiseFunctor_T>::processBaseCell(unsigned long x, unsigned long y,
                                                                          unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  const unsigned long xArray = this->getIndex(x, 0);
  const unsigned long yArray = this->getIndex(y, 1);

  auto &offsets = this->_cellOffsets[yArray][xArray];
  for (auto const &[offset, r] : offsets) {
    unsigned long otherIndex = baseIndex + offset;

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseIndex);
    } else {
      this->_cellFunctor.processCellPair(baseIndex, otherIndex, r);
    }
  }
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void PsVLC18Traversal<ParticleCell_T, PairwiseFunctor_T>::traverseParticles() {
  this->template c18Traversal</*allCells*/ false>(
      [&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(x, y, z); });
}

}  // namespace autopas
