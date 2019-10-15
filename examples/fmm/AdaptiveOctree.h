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

using AutoPasCont = autopas::AutoPas<FmmParticle, autopas::FullParticleCell<FmmParticle>>;
using Complex = std::complex<double>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

class AdaptiveOctreeNode;

class AdaptiveOctree {
 public:
  AdaptiveOctree(AutoPasCont &domain, int maxParticlesPerNode, int orderOfExpansion);
  [[nodiscard]] AutoPasCont *getDomain() const { return domain; }
  [[nodiscard]] int getMaxParticlesPerNode() const { return maxParticlesPerNode; }
  [[nodiscard]] const std::array<double, 3> &getDomainMinCorner() const { return domainMinCorner; }
  [[nodiscard]] const std::array<double, 3> &getDomainMaxCorner() const { return domainMaxCorner; }
  [[nodiscard]] const std::array<double, 3> &getDomainSize() const { return domainSize; }
  [[nodiscard]] int getNumberOfParticlesInRegion(std::array<double, 3> minCorner,
                                                 std::array<double, 3> maxCorner) const;
  [[nodiscard]] int getOrderOfExpansion() const { return orderOfExpansion; }
  [[nodiscard]] const AdaptiveOctreeNode *getRoot() const { return &*root; }

 private:
  AutoPasCont *domain;
  int maxParticlesPerNode;
  int orderOfExpansion;

 private:
  std::array<double, 3> domainMinCorner;
  std::array<double, 3> domainMaxCorner;
  std::array<double, 3> domainSize;
  std::unique_ptr<AdaptiveOctreeNode> root;
};

class AdaptiveOctreeNode {
 public:
  AdaptiveOctreeNode(AdaptiveOctree &tree, AdaptiveOctreeNode *parent, int childIndex, std::array<double, 3> minCorner,
                     std::array<double, 3> maxCorner);
  /**
   * Returns a pointer to the leaf node containing position.
   * @param position
   * @return
   */
  AdaptiveOctreeNode *findLeaf(const std::array<double, 3> &position) const;
  AdaptiveOctreeNode *findNeighbour(int x, int y, int z) const;
  AdaptiveOctreeNode *findNode(const std::array<double, 3> &position, int maxDepth) const;
  [[nodiscard]] AdaptiveOctree *getTree() const { return tree; }
  [[nodiscard]] AdaptiveOctreeNode *getParent() const { return parent; }
  [[nodiscard]] bool isLeaf() const { return _isLeaf; }
  [[nodiscard]] const std::array<double, 3> &getNodeMinCorner() const { return nodeMinCorner; }
  [[nodiscard]] const std::array<double, 3> &getNodeCenter() const { return nodeCenter; }
  [[nodiscard]] const std::array<double, 3> &getNodeMaxCorner() const { return nodeMaxCorner; }
  [[nodiscard]] const std::array<double, 3> &getNodeSize() const { return nodeSize; }
  Complex getM(int m, int n);
  void setM(int m, int n, Complex value);
  Complex getL(int m, int n);
  void setL(int m, int n, Complex value);
  [[nodiscard]] bool isZeroM() const { return _isZeroM; }
  [[nodiscard]] bool isZeroL() const { return _isZeroL; }

  void initNeighbourList();
  void initInteractionList();
  std::set<AdaptiveOctreeNode *> &getNeighbourList() { return neighbourList; }
  std::set<AdaptiveOctreeNode *> &getInteractionList() { return interactionList; }

  [[nodiscard]] AdaptiveOctreeNode &getChild(int i) const { return *child[i]; }
  [[nodiscard]] std::string getName() const {return name;}

 private:
  AdaptiveOctree *tree;
  AdaptiveOctreeNode *parent;
  int depth;
  std::string name;
  std::vector<std::unique_ptr<AdaptiveOctreeNode>> child;
  bool _isLeaf;
  std::set<AdaptiveOctreeNode *> neighbourList;
  std::set<AdaptiveOctreeNode *> interactionList;
  std::array<double, 3> nodeMinCorner;
  std::array<double, 3> nodeCenter;
  std::array<double, 3> nodeMaxCorner;
  std::array<double, 3> nodeSize;

  ComplexMatrix fmmM;
  ComplexMatrix fmmL;
  bool _isZeroM = true;
  bool _isZeroL = true;
};