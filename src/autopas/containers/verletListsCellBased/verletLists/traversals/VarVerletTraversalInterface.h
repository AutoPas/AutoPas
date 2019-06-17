/**
 * @file VarVerletTraversalInterface.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

namespace autopas {

/**
 *
 * @tparam ParticleCell Needed because all traversals have to be cell pair traversals.
 * @tparam NeighborList
 */
template <class ParticleCell, class NeighborList>
class VarVerletTraversalInterface : public CellPairTraversal<ParticleCell> {
 public:
  VarVerletTraversalInterface() : CellPairTraversal<ParticleCell>({0, 0, 0}) {}

  ~VarVerletTraversalInterface() override = default;

  virtual void iterateVerletLists(NeighborList &neighborList) = 0;

  virtual bool usesNewton3() = 0;

  void initTraversal(std::vector<ParticleCell> &cells) override {}

  void endTraversal(std::vector<ParticleCell> &cells) override {}

  virtual void initVerletTraversal(NeighborList &neighborList) = 0;

  virtual void endVerletTraversal(NeighborList &neighborList) = 0;

  virtual DataLayoutOption getDataLayout() = 0;
};

}  // namespace autopas
