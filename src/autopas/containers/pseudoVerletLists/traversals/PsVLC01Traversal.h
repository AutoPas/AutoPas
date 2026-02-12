/**
 * @file PsVLC01Traversal.h
 * @author Lars Doll
 * @date 12.01.2026
 */

#pragma once

#include "PsVLTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides the psVl_c01 traversal.
 * @tparam ParticleCell_T the type of cells.
 * @tparam PairwiseFunctor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor_T>
class PsVLC01Traversal : public C01BasedTraversal<ParticleCell_T, PairwiseFunctor_T, 3>,
                         public PsVLTraversalInterface<ParticleCell_T> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param pairwiseFunctor The functor that defines the interaction of particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length in CellBlock3D
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * in that case the interactionLength is needed!
   */
  explicit PsVLC01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor_T *pairwiseFunctor,
                            const double interactionLength, const std::array<double, 3> &cellLength,
                            DataLayoutOption dataLayout, bool useNewton3)
      : C01BasedTraversal<ParticleCell_T, PairwiseFunctor_T, 3>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                            dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength, dataLayout, useNewton3) {
    this->computeOffsets();
  }

  void traverseParticles() override;

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
    if (this->_useNewton3 == false && this->_dataLayout == DataLayoutOption::aos) {
      return true;
    }
    return false;
  }
  /**
   * Getter.
   * @return
   */
  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::psvl_c01; }

  /**
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) override;

  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  using CellOffsetsType = std::vector<std::vector<std::pair<long, std::array<double, 3>>>>;

  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  inline void processBaseCell(unsigned long x, unsigned long y, unsigned long z);

  /**
   * Pairwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellPairwise(unsigned long x, unsigned long y, unsigned long z);

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell_T, PairwiseFunctor_T,
                            /*bidirectional*/ false>
      _cellFunctor;
};

template <class ParticleCell_T, class PairwiseFunctor_T>
void PsVLC01Traversal<ParticleCell_T, PairwiseFunctor_T>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) {
  PsVLTraversalInterface<ParticleCell_T>::setOrientationList(list);
  _cellFunctor.setOrientationList(list);
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void PsVLC01Traversal<ParticleCell_T, PairwiseFunctor_T>::processBaseCell(unsigned long x, unsigned long y,
                                                                             unsigned long z) {
  processBaseCellPairwise(x, y, z);
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void PsVLC01Traversal<ParticleCell_T, PairwiseFunctor_T>::processBaseCellPairwise(unsigned long x, unsigned long y,
                                                                                     unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  for (const auto &slice : this->_cellOffsets) {
    for (auto const &[offset, r] : slice) {
      const unsigned long otherIndex = baseIndex + offset;

      if (baseIndex == otherIndex) {
        this->_cellFunctor.processCell(baseIndex);
      } else {
        this->_cellFunctor.processCellPair(baseIndex, otherIndex, r);
      }
    }
  }
}

template <class ParticleCell_T, class PairwiseFunctor_T>
inline void PsVLC01Traversal<ParticleCell_T, PairwiseFunctor_T>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(x, y, z); });
}
}  // namespace autopas
