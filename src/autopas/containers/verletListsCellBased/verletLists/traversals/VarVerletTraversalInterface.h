/**
 * @file VarVerletTraversalInterface.h
 * @author humig
 * @date 21.05.19
 */

#pragma once

namespace autopas {

template<class NeighborList>
class VarVerletTraversalInterface {
  // TODO: Have SoA as an option
 public:
  virtual ~VarVerletTraversalInterface() = default;

  virtual void iterateVerletLists(NeighborList &neighborList) = 0;

  virtual bool usesNewton3() = 0;
};

} //namespace autopas
