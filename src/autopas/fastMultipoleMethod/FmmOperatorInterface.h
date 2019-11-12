#pragma once

#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"

namespace autopas::fmm {

template <class ParticleCell>
class FmmOperatorInterface {
 public:
  virtual void P2M(FmmTreeNode<ParticleCell> &leaf) = 0;
  virtual void M2M(FmmTreeNode<ParticleCell> &parent) = 0;
  virtual void M2L(FmmTreeNode<ParticleCell> &node) = 0;
  virtual void L2L(FmmTreeNode<ParticleCell> &child) = 0;
  virtual void L2P(FmmTreeNode<ParticleCell> &leaf) = 0;
  void RunFmm(FmmTree<ParticleCell> &fmmTree) {
    // Todo
  }
};
}  // namespace autopas::fmm
