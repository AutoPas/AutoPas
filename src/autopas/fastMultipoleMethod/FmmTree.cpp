/**
 * @file FmmTree.cpp
 * @author Joachim Marin
 * @date 22.10.2019
 */

#include "FmmTree.h"

autopas::fmm::FmmTreeNode &autopas::fmm::FmmTree::setRoot(std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
  root = std::make_unique<autopas::fmm::FmmTreeNode>(*this, nullptr, boxMin, boxMax);
  return *root;
}
