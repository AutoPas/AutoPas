/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include <autopas/utils/WrapOpenMP.h>
#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The traversal is defined in the function c01Traversal and uses 1 color. Interactions between two cells are allowed
 * only if particles of the first cell are modified. This means that newton3 optimizations are NOT allowed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C01BasedTraversal : public CBasedTraversal<ParticleCell> {
  using ParticleFloatType = typename ParticleCell::ParticleType::ParticleFloatingPointType;

 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3>& dims, PairwiseFunctor* pairwiseFunctor,
                             ParticleFloatType cutoff = 1.0,
                             const std::array<ParticleFloatType, 3>& cellLength = {1.0, 1.0, 1.0})
      : CBasedTraversal<ParticleCell>(dims, cutoff, cellLength), _dataLayoutConverter(pairwiseFunctor) {}

  void initTraversal(std::vector<ParticleCell>& cells) override {
#ifdef AUTOPAS_OPENMP
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < cells.size(); ++i) {
      _dataLayoutConverter.loadDataLayout(cells[i]);
    }
  }

  void endTraversal(std::vector<ParticleCell>& cells) override {
#ifdef AUTOPAS_OPENMP
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < cells.size(); ++i) {
      _dataLayoutConverter.storeDataLayout(cells[i]);
    }
  }

 protected:
  /**
   * The main traversal of the C01Traversal.
   * This provides the structure of the loops and its parallelization.
   * @tparam LoopBody
   * @param loopBody The body of the loop as a function. Normally a lambda function, that takes as as parameters
   * (x,y,z). If you need additional input from outside, please use captures (by reference).
   */
  template <typename LoopBody>
  inline void c01Traversal(LoopBody&& loopBody);

 private:
  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, DataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
template <typename LoopBody>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::c01Traversal(
    LoopBody&& loopBody) {
  const auto offset = this->_overlap;
  const auto end = ArrayMath::sub(this->_cellsPerDimension, this->_overlap);
  this->cTraversal(std::forward<LoopBody>(loopBody), end, {1ul, 1ul, 1ul}, offset);
}
}  // namespace autopas
