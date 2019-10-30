/**
 * @file AdaptiveOctree.h
 * @date 15.10.19
 * @author Joachim Marin
 */

#pragma once

#include <memory>
#include <set>
#include <string>
#include <vector>
#include "FmmParticle.h"
#include "Math3D.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "AdaptiveOctreeNode.h"

class AdaptiveOctree {
 public:
  /**
   * Create a new Octree.
   * @param domain AutoPas object used as domain.
   * @param maxParticlesPerLeaf Maximum number of particles a leaf is supposed to contain.
   * @param orderOfExpansion Order of the multipole and local expansions.
   * @param minDepth Minimum depth of the tree.
   * @param maxDepth Maximum depth of the tree. Use -1 to have no maximum depth.
   */
  AdaptiveOctree(AutoPasCont &domain, int maxParticlesPerLeaf, int orderOfExpansion, int minDepth, int maxDepth);

  [[nodiscard]] AutoPasCont *getDomain() const { return domain; }
  [[nodiscard]] int getMaxParticlesPerNode() const { return maxParticlesPerNode; }
  [[nodiscard]] const std::array<double, 3> &getDomainMinCorner() const { return domainMinCorner; }
  [[nodiscard]] const std::array<double, 3> &getDomainMaxCorner() const { return domainMaxCorner; }
  [[nodiscard]] const std::array<double, 3> &getDomainSize() const { return domainSize; }
  [[nodiscard]] int getOrderOfExpansion() const { return orderOfExpansion; }
  [[nodiscard]] AdaptiveOctreeNode *getRoot() const { return &*root; }
  [[nodiscard]] int getMinDepth() const { return minDepth; }
  [[nodiscard]] int getMaxDepth() const { return maxDepth; }

  /**
   * Return the number of particles within a region.
   * @param minCorner
   * @param maxCorner
   * @return
   */
  [[nodiscard]] int getNumberOfParticlesInRegion(std::array<double, 3> minCorner,
                                                 std::array<double, 3> maxCorner) const;

  // Statistical information:
  int currentMaxDepth = 0;
  int numberOfLeaves = 0;
  int numberOfNodes = 0;
  int totalInteractionNodes = 0;
  int totalNearFieldNodes = 0;

 private:
  AutoPasCont *domain;
  int maxParticlesPerNode;
  int orderOfExpansion;
  std::array<double, 3> domainMinCorner;
  std::array<double, 3> domainMaxCorner;
  std::array<double, 3> domainSize;
  std::unique_ptr<AdaptiveOctreeNode> root;
  int minDepth;
  int maxDepth;
};

