/**
 * @file OctreeStaticNodeSelector.h
 * @author Johannes Spies
 * @date 01.12.2021
 */

#pragma once

#include <memory>

#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/OctreeNodeWrapper.h"

namespace autopas {
template <typename ParticleT>
class OctreeInnerNode;

/**
 * Will execute the passed function on the given root node.
 *
 * @tparam ParticleT
 * @tparam FunctionType
 * @param root
 * @param function
 * @return
 */
template <typename ParticleT, typename FunctionType>
decltype(auto) withStaticNodeType(const std::unique_ptr<OctreeNodeInterface<ParticleT>> &root, FunctionType &&function) {
  OctreeNodeInterface<ParticleT> *nodePtr = root.get();
  // TODO: These should be static casts because we already do the runtime check via hasChildren()
  if (nodePtr->hasChildren()) {
    return function(dynamic_cast<OctreeInnerNode<ParticleT> *>(nodePtr));
  } else {
    return function(dynamic_cast<OctreeLeafNode<ParticleT> *>(nodePtr));
  }
}
}  // namespace autopas
