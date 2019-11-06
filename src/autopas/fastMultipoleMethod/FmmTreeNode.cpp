/**
 * @file FmmTreeNode.cpp
 * @author Joachim Marin
 * @date 22.10.2019
 */

#include "FmmTreeNode.h"
#include "FmmTree.h"

autopas::fmm::FmmTreeNode::FmmTreeNode(autopas::fmm::FmmTree &tree, autopas::fmm::FmmTreeNode *parent, std::array<double, 3> boxMin, std::array<double, 3> boxMax)
    : tree(&tree), parent(parent), boxMin(boxMin), boxMax(boxMax) {}

void autopas::fmm::FmmTreeNode::split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin) {
  child = std::vector<std::unique_ptr<autopas::fmm::FmmTreeNode>>(2);
  child[0] = std::make_unique<autopas::fmm::FmmTreeNode>(*this->tree, this, boxMin, firstBoxMax);
  child[1] = std::make_unique<autopas::fmm::FmmTreeNode>(*this->tree, this, secondBoxMin, boxMax);
}
