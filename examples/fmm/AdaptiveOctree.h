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

using AutoPasCont = autopas::AutoPas<FmmParticle, autopas::FullParticleCell<FmmParticle>>;
using Complex = std::complex<double>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

class AdaptiveOctreeNode;

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

class AdaptiveOctreeNode {
 public:
  /**
   * Create a new node.
   * @param tree The tree the node belongs to.
   * @param parent The parent node of this node.
   * @param childIndex Which child of the parent node this node is. The name is 'parent name'+"->"
   * @param minCorner
   * @param maxCorner
   */
  AdaptiveOctreeNode(AdaptiveOctree &tree, AdaptiveOctreeNode *parent, int childIndex, std::array<double, 3> minCorner,
                     std::array<double, 3> maxCorner);

  /**
   * Returns a pointer to a neighbour of the current node. The neighbour's depth is not greater than the node's depth.
   * @param x X-Direction of the neighbour. Must be in {-1,0,1}.
   * @param y Y-Direction of the neighbour. Must be in {-1,0,1}.
   * @param z Z-Direction of the neighbour. Must be in {-1,0,1}.
   * @return A pointer to the neighbour node in the given direction.
   */
  [[nodiscard]] AdaptiveOctreeNode *findNeighbour(int x, int y, int z) const;

  /**
   * Returns a pointer to a descendant of this node containing position. If this node does not contain the position, a
   * pointer to the closest node is returned instead.
   * @param position The position the descendant should contain.
   * @param maxDepth The maximum absolute depth of the descendant.
   * @returnA A pointer to the descendant.
   */
  [[nodiscard]] AdaptiveOctreeNode *findNode(const std::array<double, 3> &position, int maxDepth) const;

  /**
   * Returns a pointer to the tree this node belongs to.
   * @return A pointer to the tree this node belongs to.
   */
  [[nodiscard]] AdaptiveOctree *getTree() const { return tree; }

  /**
   * Returns a pointer to the node's parent node. May be nullptr, if the node is the root of the tree.
   * @return A pointer to the node's parent node.
   */
  [[nodiscard]] AdaptiveOctreeNode *getParent() const { return parent; }

  /**
   * Returns true, if the node has no child nodes. Returns false otherwise.
   * @return
   */
  [[nodiscard]] bool isLeaf() const { return _isLeaf; }

  [[nodiscard]] const std::array<double, 3> &getNodeMinCorner() const { return nodeMinCorner; }
  [[nodiscard]] const std::array<double, 3> &getNodeCenter() const { return nodeCenter; }
  [[nodiscard]] const std::array<double, 3> &getNodeMaxCorner() const { return nodeMaxCorner; }
  [[nodiscard]] const std::array<double, 3> &getNodeSize() const { return nodeSize; }

  /**
   * Get the multipole expansion coefficient.
   * @param m
   * @param n
   * @return
   */
  Complex getM(int m, int n);

  /**
   * Set the multipole expansion coefficient.
   * @param m
   * @param n
   * @param value
   */
  void setM(int m, int n, Complex value);

  /**
   * Get the local expansion coefficient.
   * @param m
   * @param n
   * @return
   */
  Complex getL(int m, int n);

  /**
   * Get the local expansion coefficient.
   * @param m
   * @param n
   * @param value
   */
  void setL(int m, int n, Complex value);
  [[nodiscard]] bool isZeroM() const { return _isZeroM; }
  [[nodiscard]] bool isZeroL() const { return _isZeroL; }
  [[nodiscard]] int getDepth() const { return depth; }

  void initNeighbourList();
  void initNearFieldList();
  void initInteractionList();
  std::set<AdaptiveOctreeNode *> &getNeighbourList() { return neighbourList; }
  std::set<AdaptiveOctreeNode *> &getNearFieldList() { return nearFieldList; }
  std::set<AdaptiveOctreeNode *> &getInteractionList() { return interactionList; }

  // String representations for the lists.
  std::string neighbourListString;
  std::string nearFieldListString;
  std::string interactionListString;

  /**
   * Return the i-th child node.
   * @param i
   * @return
   */
  [[nodiscard]] AdaptiveOctreeNode *getChild(int i) const { return &*child[i]; }

  /**
   * Return the name of the node.
   * @return
   */
  [[nodiscard]] std::string getName() const { return name; }

 private:
  AdaptiveOctree *tree;
  AdaptiveOctreeNode *parent;
  int depth;
  std::string name;
  std::vector<std::unique_ptr<AdaptiveOctreeNode>> child;
  bool _isLeaf;
  std::set<AdaptiveOctreeNode *> neighbourList;
  std::set<AdaptiveOctreeNode *> nearFieldList;
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