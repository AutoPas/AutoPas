/**
 * @file FmmTreeNode.h
 * @author Joachim Marin
 * @date 22.10.2019
 */
#pragma once

#include <array>
#include <memory>
#include <vector>

class FmmTree;

class FmmTreeNode {
 public:
  FmmTreeNode(FmmTree &tree, FmmTreeNode *parent, std::array<double, 3> boxMin, std::array<double, 3> boxMax);
  void split(std::array<double, 3> firstBoxMax, std::array<double, 3> secondBoxMin);
  FmmTreeNode &getChild(int index) { return *child[index]; }
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const { return boxMin; }
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const { return boxMax; }

 private:
  FmmTree *tree;
  FmmTreeNode *parent;
  std::array<double, 3> boxMin;
  std::array<double, 3> boxMax;
  std::vector<std::unique_ptr<FmmTreeNode>> child;
};
