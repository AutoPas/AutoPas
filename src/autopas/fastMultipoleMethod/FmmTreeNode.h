/**
 * @file FmmTreeNode.h
 * @author Joachim Marin
 * @date 22.10.2019
 */
#pragma once

#include <array>
#include <complex>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>
#include "autopas/utils/ArrayMath.h"

namespace autopas::fmm {

// forward declaration
class FmmTree;

class FmmTreeNode {
 public:
  std::string name;
  FmmTreeNode(FmmTree &tree, FmmTreeNode *parent, std::array<double, 3> boxMin, std::array<double, 3> boxMax)
      : tree(&tree),
        parent(parent),
        boxMin(boxMin),
        boxMax(boxMax),
        boxCenter(autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::add(boxMin, boxMax), 0.5)),
        sphereRadius(autopas::utils::ArrayMath::L2Norm(autopas::utils::ArrayMath::sub(boxMax, boxCenter))) {
    if (parent == nullptr) {
      depth = 0;
    } else {
      depth = parent->depth + 1;
    }

    std::cout << "sphereRadius = " << sphereRadius << std::endl;
  }
  void split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin) {
    if (children.empty()) {
      children.emplace_back(*this->tree, this, boxMin, firstBoxMax);
      children.emplace_back(*this->tree, this, secondBoxMin, boxMax);
    } else {
      autopas::utils::ExceptionHandler::exception("trying to split an already split FmmTreeNode");
    }
  }

  void initRec(FmmTreeNode &otherNode) {
    double r1 = sphereRadius;
    double r2 = otherNode.sphereRadius;
    double maxRadius = std::max(r1, r2);
    double distance = autopas::utils::ArrayMath::L2Norm(autopas::utils::ArrayMath::sub(boxCenter, otherNode.boxCenter));

    // Well separated, if distance > 3 * maxRadius

    if (distance > 3 * maxRadius) {
      // Check if already well separated for parent
      if (parent != nullptr) {
        double rParent = parent->sphereRadius;
        double maxRadiusParent = std::max(rParent, r2);
        double distanceParent =
            autopas::utils::ArrayMath::L2Norm(autopas::utils::ArrayMath::sub(parent->boxCenter, otherNode.boxCenter));
        if (distanceParent > 3 * maxRadiusParent) {
          // Well separated for parent
        } else {
          // Well separated for this node

          // Check if a descendant is well separated for the parent
          // ...

          interactionList.insert(&otherNode);
        }
      } else {
        // Well separated for this node
        interactionList.insert(&otherNode);
      }
    } else {
      {
        // Node is not a leaf, so try to find more refined nodes to add to the interaction list.
        if (!otherNode.isLeaf()) {
          initRec(otherNode.getChild(0));
          initRec(otherNode.getChild(1));
        } else {
          if (isLeaf()) {
            nearFieldList.insert(&otherNode);
          }
        }
      }
    }
  }

  void init(FmmTreeNode &otherNode) {
    initRec(otherNode);
    // Debug print near field list and interaction list
    std::cout << "Node = " << name << std::endl;
    std::cout << "NearFieldList = " << std::endl;
    for (auto node : nearFieldList) {
      std::cout << node->name << " , ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "InteractionList = " << std::endl;
    for (auto node : interactionList) {
      std::cout << node->name << " , ";
    }
    std::cout << std::endl << std::endl;
    if (!isLeaf()) {
      getChild(0).init(otherNode);
      getChild(1).init(otherNode);
    }
  }

  FmmTreeNode &getChild(int index) { return children[index]; }
  [[nodiscard]] std::array<double, 3> getBoxMin() const { return boxMin; }
  [[nodiscard]] std::array<double, 3> getBoxMax() const { return boxMax; }
  [[nodiscard]] std::array<double, 3> getBoxCenter() const { return boxCenter; }
  [[nodiscard]] bool isLeaf() const { return _isLeaf; }
  [[nodiscard]] long getDepth() const { return depth; }
  [[nodiscard]] FmmTree &getTree() const { return *tree; }
  [[nodiscard]] FmmTreeNode *getParent() const { return parent; }
  [[nodiscard]] double getSphereRadius() const { return sphereRadius; }

  void makeLeaf() { _isLeaf = true; }

  void initCoefficients(long orderOfExpansion) {
    coefficientM = std::vector<std::vector<std::complex<double>>>(
        orderOfExpansion * 2 + 1, std::vector<std::complex<double>>(orderOfExpansion + 1, 0));
    coefficientL = std::vector<std::vector<std::complex<double>>>(
        orderOfExpansion * 2 + 1, std::vector<std::complex<double>>(orderOfExpansion + 1, 0));
  }

  void setM(long m, long n, std::complex<double> value) {
    assert(!__isnan(value.real()) && !__isnan(value.imag()));
    long offset = m <= 0 ? 0 : -1;
    coefficientM.at(2 * std::abs(m) + offset).at(n) = value;
  }
  void setL(long m, long n, std::complex<double> value) {
    assert(!__isnan(value.real()) && !__isnan(value.imag()));
    long offset = m <= 0 ? 0 : -1;
    coefficientL.at(2 * std::abs(m) + offset).at(n) = value;
  }
  std::complex<double> getM(long m, long n) {
    long offset = m <= 0 ? 0 : -1;
    long indexM = 2 * std::abs(m) + offset;
    if (indexM < 0 || indexM >= static_cast<int>(coefficientM.size())) {
      return 0.0;
    }
    if (n < 0 || n >= static_cast<int>(coefficientM[0].size())) {
      return 0.0;
    }
    return coefficientM.at(indexM).at(n);
  }
  std::complex<double> getL(long m, long n) {
    long offset = m <= 0 ? 0 : -1;
    long indexL = 2 * std::abs(m) + offset;
    if (indexL < 0 || indexL >= static_cast<int>(coefficientL.size())) {
      return 0.0;
    }
    if (n < 0 || n >= static_cast<int>(coefficientL[0].size())) {
      return 0.0;
    }
    return coefficientL.at(indexL).at(n);
  }
  std::unordered_set<FmmTreeNode *> &getInteractionList() { return interactionList; }
  std::unordered_set<FmmTreeNode *> &getNearFieldList() { return nearFieldList; }

 private:
  FmmTree *tree;
  FmmTreeNode *parent;
  std::array<double, 3> boxMin;
  std::array<double, 3> boxMax;
  std::array<double, 3> boxCenter;
  long depth;
  std::vector<FmmTreeNode> children;
  bool _isLeaf = false;
  // The radius of the sphere that fully contains the node.
  double sphereRadius;
  std::unordered_set<FmmTreeNode *> interactionList;
  std::unordered_set<FmmTreeNode *> nearFieldList;

  std::vector<std::vector<std::complex<double>>> coefficientM;
  std::vector<std::vector<std::complex<double>>> coefficientL;
};
}  // namespace autopas::fmm
