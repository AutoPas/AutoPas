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
   * Sets the neighbor list to traverse.
   * @param neighborList the neighbor list to traverse.
   */
  virtual void setNeighborListToTraverse(NeighborList &neighborList) { _neighborList = &neighborList; }

 protected:
  /**
   * The neighbor list to traverse.
   */
  NeighborList *_neighborList;
};

}  // namespace autopas
