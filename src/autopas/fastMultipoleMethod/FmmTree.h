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
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/particles/FmmParticle.h"

namespace autopas::fmm {
template <class ParticleCell>
class FmmTree {
 public:
  explicit FmmTree(ParticleContainerInterface<ParticleCell> *container) : container(container) {}
  FmmTreeNode<ParticleCell> &getRoot() { return *root; }

  [[nodiscard]] ParticleContainerInterface<ParticleCell> *getContainer() const { return container; }
  FmmTreeNode<ParticleCell> &setRoot(std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
    root = std::make_unique<autopas::fmm::FmmTreeNode<ParticleCell>>(*this, nullptr, boxMin, boxMax);
    return *root;
  }

 private:
  std::unique_ptr<FmmTreeNode<ParticleCell>> root;
  ParticleContainerInterface<ParticleCell> *container;
};

}  // namespace autopas::fmm
