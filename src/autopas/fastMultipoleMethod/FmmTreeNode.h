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
  FmmTreeNode(FmmTreeNode *parent, std::array<double, 3> boxMin, std::array<double, 3> boxMax)
      : parent(parent),
        boxMin(boxMin),
        boxMax(boxMax),
        boxCenter(autopas::utils::ArrayMath::mulScalar(autopas::utils::ArrayMath::add(boxMin, boxMax), 0.5)),
        sphereRadius(autopas::utils::ArrayMath::L2Norm(autopas::utils::ArrayMath::sub(boxMax, boxCenter))) {
    if (parent == nullptr) {
      depth = 0;
    } else {
      depth = parent->depth + 1;
    }

    // std::cout << "sphereRadius = " << sphereRadius << std::endl;
  }
  FmmTreeNode(const FmmTreeNode &) = delete;
  FmmTreeNode &operator=(const FmmTreeNode &) = delete;
  FmmTreeNode(FmmTreeNode &&) = default;
  FmmTreeNode &operator=(FmmTreeNode &&) = default;

  void split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin) {
    if (children.empty()) {
      children.emplace_back(this, boxMin, firstBoxMax);
      children.emplace_back(this, secondBoxMin, boxMax);
    } else {
      autopas::utils::ExceptionHandler::exception("trying to split an already split FmmTreeNode");
    }
  }

  void insertLeaves(FmmTreeNode &otherNode) {
    if (otherNode.isLeaf()) {
      if (ancestorInteractionList.find(&otherNode) == ancestorInteractionList.end()) {
        interactionList.insert(&otherNode);
      }
    } else {
      insertLeaves(otherNode.getChild(0));
      insertLeaves(otherNode.getChild(1));
    }
  }

  bool checkInteract(FmmTreeNode &otherNode) {
    /*std::cout << "checkInteract:" << std::endl;
    std::cout << name << std::endl;
    std::cout << otherNode.name << std::endl;*/
    double r1 = sphereRadius;
    double r2 = otherNode.sphereRadius;
    double maxRadius = std::max(r1, r2);
    double distance = autopas::utils::ArrayMath::L2Norm(autopas::utils::ArrayMath::sub(boxCenter, otherNode.boxCenter));
    // std::cout << "maxRadius = " << maxRadius << std::endl;
    // std::cout << "distance = " << distance << std::endl;
    return distance > 3 * maxRadius;
  }
  double getLargestSize() {
    auto boxSize = utils::ArrayMath::sub(boxMax, boxMin);
    double largestSize = boxSize[0];
    for (size_t i = 1; i < 3; ++i) {
      if (boxSize[i] > largestSize) {
        largestSize = boxSize[i];
      }
    }
    return largestSize;
  }

  void initLists(FmmTreeNode &otherNode) {
    if (checkInteract(otherNode)) {
      interactionList.insert(&otherNode);
    } else {
      if (getLargestSize() > otherNode.getLargestSize()) {
        if (not isLeaf()) {
          getChild(0).initLists(otherNode);
          getChild(1).initLists(otherNode);
        } else {
          nearFieldList.insert(&otherNode);
        }
      } else {
        if (not otherNode.isLeaf()) {
          initLists(otherNode.getChild(0));
          initLists(otherNode.getChild(1));
        } else {
          nearFieldList.insert(&otherNode);
        }
      }
    }
  }

  void printLists() {
    if (sphereRadius < 2) {
      return;
    }
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
    if (not isLeaf()) {
      getChild(0).printLists();
      getChild(1).printLists();
    }
  }

  FmmTreeNode &getChild(int index) { return children[index]; }
  [[nodiscard]] std::array<double, 3> getBoxMin() const { return boxMin; }
  [[nodiscard]] std::array<double, 3> getBoxMax() const { return boxMax; }
  [[nodiscard]] std::array<double, 3> getBoxCenter() const { return boxCenter; }
  [[nodiscard]] bool isLeaf() const { return _isLeaf; }
  [[nodiscard]] long getDepth() const { return depth; }
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
    assert(not __isnan(value.real()) && not __isnan(value.imag()));
    long offset = m <= 0 ? 0 : -1;
    coefficientM.at(2 * std::abs(m) + offset).at(n) = value;
  }
  void setL(long m, long n, std::complex<double> value) {
    assert(not __isnan(value.real()) && not __isnan(value.imag()));
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
  std::unordered_set<FmmTreeNode *> ancestorInteractionList;

  std::vector<std::vector<std::complex<double>>> coefficientM;
  std::vector<std::vector<std::complex<double>>> coefficientL;
};
}  // namespace autopas::fmm
