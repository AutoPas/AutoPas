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
template <typename Particle>
class OctreeInnerNode;

/**
 * Will execute the passed function on the given root node.
 *
 * @tparam Particle
 * @tparam FunctionType
 * @param root
 * @param function
 * @return
 */
template <typename Particle, typename FunctionType>
decltype(auto) withStaticNodeType(const std::unique_ptr<OctreeNodeInterface<Particle>> &root, FunctionType &&function) {
  OctreeNodeInterface<Particle> *nodePtr = root.get();
  if (nodePtr->hasChildren()) {
    return function(dynamic_cast<OctreeInnerNode<Particle> *>(nodePtr));
  } else {
    return function(dynamic_cast<OctreeLeafNode<Particle> *>(nodePtr));
  }
}
}  // namespace autopas
