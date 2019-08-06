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
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class C01TraversalAdaptive : public CellPairTraversal<ParticleCell>, public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param info
   */
  explicit C01TraversalAdaptive(PairwiseFunctor *pairwiseFunctor, std::unique_ptr<TraversalSelectorInfo> info)
      : CellPairTraversal<ParticleCell>(info->dims),
        _cellFunctor(pairwiseFunctor, info->interactionLength),
        _dataLayoutConverter(pairwiseFunctor),
        _octree(*static_cast<Octree<typename ParticleCell::ParticleType, ParticleCell> *>(
            reinterpret_cast<TraversalSelectorInfoAdaptive *>(info.get())->octree)) {
    std::cout << static_cast<std::string>(_octree) << std::endl;
  }

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseParticlePairs() override;

  /**
   * load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    if (this->_cells) {
      auto &cells = *(this->_cells);
#ifdef AUTOPAS_OPENMP
      // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   * @return
   */
  bool isApplicable() const override { return dataLayout != DataLayoutOption::cuda and not useNewton3; }

  TraversalOption getTraversalType() const override { return TraversalOption::c01Adaptive; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

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

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;

  Octree<typename ParticleCell::ParticleType, ParticleCell> &_octree;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
inline void C01TraversalAdaptive<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, internal::OctreeNode<typename ParticleCell::ParticleType, ParticleCell> &node) {
  auto &leaf = dynamic_cast<internal::OctreeExternalNode<typename ParticleCell::ParticleType, ParticleCell> &>(node);
  auto &baseCell{leaf.getCell()};
  for (auto const &[neighborIndex, r] : leaf.getNeighbors()) {
    _cellFunctor.processCellPair(baseCell, cells[neighborIndex], r);
  }
  _cellFunctor.processCell(baseCell);
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
inline void C01TraversalAdaptive<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception("The adaptive C01 traversal cannot work with enabled newton3");
  }
  _octree.apply([&](auto &node) { this->processBaseCell(cells, node); }, internal::ExecutionPolicy::par);
}

}  // namespace autopas
