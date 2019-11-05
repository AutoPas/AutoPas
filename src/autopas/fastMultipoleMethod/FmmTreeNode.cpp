/**
 * @file FmmTreeNode.cpp
 * @author Joachim Marin
 * @date 22.10.2019
 */

#include "FmmTreeNode.h"
#include "FmmTree.h"

FmmTreeNode::FmmTreeNode(FmmTree &tree, FmmTreeNode *parent, std::array<double, 3> boxMin, std::array<double, 3> boxMax)
    : tree(&tree), parent(parent), boxMin(boxMin), boxMax(boxMax) {}

void FmmTreeNode::split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin) {
  child = std::vector<std::unique_ptr<FmmTreeNode>>(2);
  child[0] = std::make_unique<FmmTreeNode>(*this->tree, this, boxMin, firstBoxMax);
  child[1] = std::make_unique<FmmTreeNode>(*this->tree, this, secondBoxMin, boxMax);
}
