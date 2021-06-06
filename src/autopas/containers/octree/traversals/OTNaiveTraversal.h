/**
 * @file OTNaiveTraversal.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/traversals/OTTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class OTNaiveTraversal : public CellPairTraversal<OctreeLeafNode<Particle>>,
                         public OTTraversalInterface<OctreeNodeWrapper<Particle>> {
 public:
  //using ParticleCell = OctreeNodeWrapper<Particle>;
  using ParticleCell = OctreeLeafNode<Particle>;

  // TODO(johannes): The TraversalSelector passes the interactionLength as the cutoff value: Keep in mind when
  // implementing...
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

  void initTraversal() override { printf("Johannes' OTNaiveTraversal::initTraversal\n"); }

  void endTraversal() override { printf("Johannes' OTNaiveTraversal::endTraversal\n"); }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override {
    // Gather all leaves
    std::vector<OctreeLeafNode<Particle> *> leaves;
    auto *wrapper = dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[0]);
    wrapper->appendAllLeaves(leaves);

    // Get neighboring cells for each leaf
    for (OctreeLeafNode<Particle> *leaf : leaves) {
      // Get all face neighbors
      for (Face *face = getFaces(); *face != O; ++face) {
        OctreeNodeInterface<Particle> *neighbor = leaf->GTEQ_FACE_NEIGHBOR(*face);
        if (neighbor) {
          auto neighborLeaf = dynamic_cast<OctreeLeafNode<Particle> *>(neighbor);
          OctreeLeafNode<Particle> &leafRef = *leaf;
          OctreeLeafNode<Particle> &neighborLeafRef = *neighborLeaf;
          _cellFunctor.processCellPair(leafRef, neighborLeafRef);
        }
      }

      // Get all edge neighbors
      for (Edge *edge = getEdges(); *edge != OO; ++edge) {
        OctreeNodeInterface<Particle> *neighbor = leaf->GTEQ_EDGE_NEIGHBOR(*edge);
        if (neighbor) {
          auto neighborLeaf = dynamic_cast<OctreeLeafNode<Particle> *>(neighbor);
          OctreeLeafNode<Particle> &leafRef = *leaf;
          OctreeLeafNode<Particle> &neighborLeafRef = *neighborLeaf;
          _cellFunctor.processCellPair(leafRef, neighborLeafRef);
        }
      }

      // Get all face neighbors
      for (Vertex *vertex = VERTICES(); *vertex != OOO; ++vertex) {
        OctreeNodeInterface<Particle> *neighbor = leaf->GTEQ_VERTEX_NEIGHBOR(*vertex);
        if (neighbor) {
          auto neighborLeaf = dynamic_cast<OctreeLeafNode<Particle> *>(neighbor);
          OctreeLeafNode<Particle> &leafRef = *leaf;
          OctreeLeafNode<Particle> &neighborLeafRef = *neighborLeaf;
          _cellFunctor.processCellPair(leafRef, neighborLeafRef);
        }
      }
    }
  }

  void setCells(std::vector<OctreeNodeWrapper<Particle>> *cells) override { _cells = cells; }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<Particle, ParticleCell, PairwiseFunctor, dataLayout, useNewton3, true> _cellFunctor;
  /*internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                        true>
      _cellFunctor;*/

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;

  std::vector<OctreeNodeWrapper<Particle>> *_cells;
};

}  // namespace autopas