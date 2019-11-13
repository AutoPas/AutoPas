/**
 * @file FmmTreeNode.h
 * @author Joachim Marin
 * @date 22.10.2019
 */
#pragma once

#include <array>
#include <complex>
#include <memory>
#include <set>
#include <vector>
#include "autopas/utils/ArrayMath.h"

namespace autopas::fmm {

// forward declaration
template <class ParticleCell>
class FmmTree;

template <class ParticleCell>
class FmmTreeNode {
 public:
  FmmTreeNode(FmmTree<ParticleCell> &tree, FmmTreeNode *parent, std::array<double, 3> boxMin,
              std::array<double, 3> boxMax)
      : tree(&tree), parent(parent), boxMin(boxMin), boxMax(boxMax) {
    if (parent == nullptr) {
      depth = 0;
    } else {
      depth = parent->depth + 1;
    }
    boxCenter = autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::add(boxMin, boxMax), 0.5);
    _isOctreeNode = (depth % 3 == 0);
  }
  void split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin) {
    if (children.empty()) {
      children.emplace_back(*this->tree, this, boxMin, firstBoxMax);
      children.emplace_back(*this->tree, this, boxMin, secondBoxMin);
    } else {
      autopas::utils::ExceptionHandler::exception("trying to split an already split FmmTreeNode");
    }
  }

  FmmTreeNode &getChild(int index) { return children[index]; }
  [[nodiscard]] std::array<double, 3> getBoxMin() const { return boxMin; }
  [[nodiscard]] std::array<double, 3> getBoxMax() const { return boxMax; }
  [[nodiscard]] std::array<double, 3> getBoxCenter() const { return boxCenter; }
  [[nodiscard]] bool isLeaf() const { return _isLeaf; }
  [[nodiscard]] bool isOctreeNode() const { return _isOctreeNode; }
  [[nodiscard]] bool isOctreeLeaf() const { return _isOctreeLeaf; }
  [[nodiscard]] long getDepth() const { return depth; }
  [[nodiscard]] FmmTree<ParticleCell> &getTree() const { return *tree; }

  void makeLeaf() {
    _isLeaf = true;
    makeOctreeLeaf();
  }
  FmmTreeNode &getOctreeChild(unsigned long index) {
    return getChild((index & 0x4u) > 1u).getChild((index & 0x2u) > 1u).getChild((index & 0x1u) > 1u);
  }
  FmmTreeNode &getOctreeParent() { return *parent->parent->parent; }

  void setM(long m, long n, std::complex<double> value) {
    // Todo
  }
  void setL(long m, long n, std::complex<double> value) {
    // Todo
  }
  std::complex<double> getL(long m, long n) {
    // Todo
    return 0.0;
  }
  std::complex<double> getM(long m, long n) {
    // Todo
    return 0.0;
  }
  std::set<FmmTreeNode *> &getInteractionList() {
    // current warning: return sth. proper that is not on the stack.
    // Todo
    static std::set<FmmTreeNode *> set;
    return set;
  }

 private:
  void makeOctreeLeaf() {
    if (_isOctreeNode) {
      _isOctreeLeaf = true;
    } else {
      if (parent != nullptr) {
        parent->makeOctreeLeaf();
      }
    }
  }
  FmmTree<ParticleCell> *tree;
  FmmTreeNode *parent;
  std::array<double, 3> boxMin;
  std::array<double, 3> boxMax;
  std::array<double, 3> boxCenter;
  long depth;
  std::vector<FmmTreeNode> children;
  bool _isLeaf = false;
  bool _isOctreeNode = false;
  bool _isOctreeLeaf = false;
};
}  // namespace autopas::fmm
