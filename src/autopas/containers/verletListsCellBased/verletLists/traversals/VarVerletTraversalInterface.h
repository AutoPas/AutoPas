/**
 * @file VarVerletTraversalInterface.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

namespace autopas {

/**
 * Interface for all traversals for VarVerletLists containers.
 * @tparam ParticleCell Needed because all traversals have to be cell pair traversals.
 * @tparam NeighborList The neighbor list type of the container to traverse.
 */
template <class ParticleCell, class NeighborList>
class VarVerletTraversalInterface : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of VarVerletTraversalInterface.
   */
  VarVerletTraversalInterface() : CellPairTraversal<ParticleCell>({0, 0, 0}) {}

  /**
   * Virtual default destructor.
   */
  ~VarVerletTraversalInterface() override = default;

  /**
   * Iterates over the given neighbor list.
   * @param neighborList The neighbor list to iterate over.
   */
  virtual void iterateVerletLists(NeighborList &neighborList) = 0;

  /**
   * Returns if this traversal uses newton 3.
   * @return True, if this traversal uses newton 3, false otherwise.
   */
  virtual bool usesNewton3() = 0;

  /**
   * Empty body. Just there to fulfill CellPairTraversal interface!
   * @param cells
   */
  void initTraversal(std::vector<ParticleCell> &cells) override {}

  /**
   * Empty body. Just there to fulfill CellPairTraversal interface!
   * @param cells
   */
  void endTraversal(std::vector<ParticleCell> &cells) override {}

  /**
   * Initializes the traversal using the neighbor list to iterate over. Should always be called before
   * iterateVerletLists().
   * @param neighborList The neighbor list to initialize from.
   */
  virtual void initVerletTraversal(NeighborList &neighborList) = 0;

  /**
   * Finalizes the traversal using the neighbor list that was iterated over. Should always be called after
   * iterateVerletLists().
   * @param neighborList The neighbor list to initialize from.
   */
  virtual void endVerletTraversal(NeighborList &neighborList) = 0;

  /**
   * Returns the data layout used by this traversal.
   * @return the data layout used by this traversal.
   */
  virtual DataLayoutOption getDataLayout() = 0;
};

}  // namespace autopas
