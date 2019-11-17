#pragma once

#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"

namespace autopas::fmm {

template <class Particle, class ParticleCell>
class FmmOperatorInterface {
 public:
  virtual void P2M(FmmTreeNode &leaf, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void M2M(FmmTreeNode &parent) = 0;
  virtual void M2L(FmmTreeNode &node) = 0;
  virtual void L2L(FmmTreeNode &child) = 0;
  virtual void L2P(FmmTreeNode &leaf, AutoPas<Particle, ParticleCell> &container) = 0;
  void RunFmm(FmmTree &fmmTree, AutoPas<Particle, ParticleCell> &container) {
    // Todo
  }
};
}  // namespace autopas::fmm
