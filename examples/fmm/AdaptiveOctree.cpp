/**
 * @file AdaptiveOctree.cpp
 * @date 15.10.19
 * @author Joachim Marin
 */

#include "AdaptiveOctree.h"

AdaptiveOctree::AdaptiveOctree(AutoPasCont &domain, int maxParticlesPerNode, int orderOfExpansion, int minDepth,
                               int maxDepth)
    : domain(&domain),
      maxParticlesPerNode(maxParticlesPerNode),
      orderOfExpansion(orderOfExpansion),
      domainMinCorner(domain.getBoxMin()),
      domainMaxCorner(domain.getBoxMax()),
      domainSize(autopas::ArrayMath::sub(domainMaxCorner, domainMinCorner)),
      minDepth(minDepth),
      maxDepth(maxDepth) {
  root = std::make_unique<AdaptiveOctreeNode>(*this, nullptr, 0, domainMinCorner, domainMaxCorner);
  root->initNeighbourList();
  root->initNearFieldList();
  root->initInteractionList();
}

int AdaptiveOctree::getNumberOfParticlesInRegion(std::array<double, 3> minCorner,
                                                 std::array<double, 3> maxCorner) const {
  int count = 0;
  for (auto iter = this->domain->getRegionIterator(minCorner, maxCorner); iter.isValid(); ++iter) {
    count++;
  }
  return count;
}
