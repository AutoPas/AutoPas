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
#include "autopas/particles/FmmParticle.h"

namespace autopas::fmm {
class FmmTree {
 public:
  FmmTree() = default;
  FmmTree(const FmmTree &) = delete;
  FmmTree &operator=(const FmmTree &) = delete;
  FmmTree(FmmTree &&) = default;
  FmmTree &operator=(FmmTree &&) = default;

  FmmTreeNode &getRoot() { return *root; }

  FmmTreeNode &setRoot(std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
    root = std::make_unique<autopas::fmm::FmmTreeNode>(nullptr, boxMin, boxMax);
    return *root;
  }

 private:
  std::unique_ptr<FmmTreeNode> root;
};

}  // namespace autopas::fmm
