/**
 * @file OTNaiveTraversal.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class OTNaiveTraversal : public CellPairTraversal<ParticleCell> {
 public:
  // TODO(johannes): The TraversalSelector passes the interactionLength as the cutoff value: Keep in mind when implementing...
  /**
   * Constructor for the Octree traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the octree traversal, please don't use the interaction length here.)
   */
  explicit OTNaiveTraversal(PairwiseFunctor *pairwiseFunctor, double cutoff)
    : CellPairTraversal<ParticleCell>({2, 1, 1}),
      _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
      _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_naive; }

  [[nodiscard]] bool isApplicable() const override {
    int nDevices = 0;
    #if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
    #endif
    if (dataLayout == DataLayoutOption::cuda)
      return nDevices > 0;
    else
      return true;
  }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    printf("Johannes' OTNaiveTraversal::initTraversal\n");
  }

  void endTraversal() override {
    printf("Johannes' OTNaiveTraversal::endTraversal\n");
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override {
    printf("Johannes' OTNaiveTraversal::traverseParticlePairs\n");
  }

private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
      true>
      _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

}