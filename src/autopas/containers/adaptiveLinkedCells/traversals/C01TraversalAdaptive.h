/**
 * @file C01TraversalAdaptive.h
 * @author C.Menges
 * @date 25.06.2018
 */

#pragma once

#include "autopas/containers/adaptiveLinkedCells/Octree.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/selectors/TraversalSelectorInfoAdaptive.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the adaptive c01 traversal.
 *
 * The traversal uses the c01 base step performed on every single cell inside an octree structure.
 * Due to geometrical reasons, c01 cannot apply newton3 optimization!
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 * @tparam combineSoA
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class C01TraversalAdaptive : public CellPairTraversal<ParticleCell, dataLayout, useNewton3>,
                             public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff radius
   * @param cellLength min. cell length
   */
  explicit C01TraversalAdaptive(PairwiseFunctor *pairwiseFunctor,
                                std::unique_ptr<TraversalSelectorInfo<ParticleCell>> info)
      : CellPairTraversal<ParticleCell, dataLayout, useNewton3>(info->dims),
        _cellFunctor(pairwiseFunctor, info->cutoff),
        _pairwiseFunctor(pairwiseFunctor),
        _dataLayoutConverter(pairwiseFunctor) /*,
         _octree(reinterpret_cast<TraversalSelectorInfoAdaptive<ParticleCell>*>(info.get())->octree)*/
  {}

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * load Data Layouts required for this Traversal.
   * @param cells where the data should be loaded
   */
  void initTraversal(std::vector<ParticleCell> &cells) override {
#ifdef AUTOPAS_OPENMP
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < cells.size(); ++i) {
      _dataLayoutConverter.loadDataLayout(cells[i]);
    }
  }

  /**
   * write Data to AoS.
   * @param cells for which the data should be written back
   */
  void endTraversal(std::vector<ParticleCell> &cells) override {
#ifdef AUTOPAS_OPENMP
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < cells.size(); ++i) {
      _dataLayoutConverter.storeDataLayout(cells[i]);
    }
  }

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   * @return
   */
  bool isApplicable() const override { return dataLayout != DataLayoutOption::cuda and not useNewton3; }

  TraversalOption getTraversalType() const override { return TraversalOption::c01Adaptive; }

 private:
  /**
   * Computes all interactions between the base
   * cell, which is represented by an octree node, and adjacent cells.
   * @param cells vector of all cells.
   * @param node The current node in the octree.
   */
  inline void processBaseCell(std::vector<ParticleCell> &cells,
                              internal::OctreeNode<typename ParticleCell::ParticleType, ParticleCell> &node);

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, false, false>
      _cellFunctor;

  PairwiseFunctor *_pairwiseFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
inline void C01TraversalAdaptive<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, internal::OctreeNode<typename ParticleCell::ParticleType, ParticleCell> &node) {
  /// @todo add implementation
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
inline void C01TraversalAdaptive<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception("The adaptive C01 traversal cannot work with enabled newton3");
  }
  /// @todo add implementation

  _cellFunctor.processCell(cells[0]);
}

}  // namespace autopas
