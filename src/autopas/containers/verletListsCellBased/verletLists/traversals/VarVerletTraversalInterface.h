/**
 * @file VarVerletTraversalInterface.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

namespace autopas {

/**
 * Interface for all traversals for VarVerletLists containers.
 * @tparam NeighborList The neighbor list type of the container to traverse.
 */
template <class NeighborList>
class VarVerletTraversalInterface {
 public:
  /**
   * Constructor of VarVerletTraversalInterface.
   */
  VarVerletTraversalInterface() = default;

  /**
   * Virtual default destructor.
   */
  virtual ~VarVerletTraversalInterface() = default;

  /**
   * Iterates over the given neighbor list.
   * @param neighborList The neighbor list to iterate over.
   */
  virtual void iterateVerletLists(NeighborList &neighborList) = 0;

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
};

}  // namespace autopas
