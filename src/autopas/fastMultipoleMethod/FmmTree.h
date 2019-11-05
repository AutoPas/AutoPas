/**
 * @file FmmTree.h
 * @author Joachim Marin
 * @date 22.10.2019
 */

#pragma once

#include <array>
#include <memory>
#include <vector>
#include "FmmTreeNode.h"

class FmmTree {
 public:
  FmmTree() = default;
  FmmTreeNode &getRoot() { return *root; }
  FmmTreeNode &setRoot(std::array<double, 3> boxMin, std::array<double, 3> boxMax);

 private:
  std::unique_ptr<FmmTreeNode> root;
};


