#pragma once

#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"

namespace autopas::fmm {

template <class Particle, class ParticleCell>
class FmmOperatorInterface {
 private:
  virtual void init(long orderOfExpansion) = 0;
  virtual void P2M(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void M2M(FmmTreeNode &parent, long orderOfExpansion) = 0;
  virtual void M2L(FmmTreeNode &node, long orderOfExpansion) = 0;
  virtual void L2L(FmmTreeNode &child, long orderOfExpansion) = 0;
  virtual void L2P(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void NearField(FmmParticle &p1, FmmParticle &p2) = 0;

  template <typename Function>
  void upwardPass(FmmTreeNode &node, Function function) {
    if (not node.isLeaf()) {
      upwardPass(node.getChild(0), function);
      upwardPass(node.getChild(1), function);
    }
    function(node);
  }
  template <typename Function>
  void downwardPass(FmmTreeNode &node, Function function) {
    function(node);
    if (not node.isLeaf()) {
      downwardPass(node.getChild(0), function);
      downwardPass(node.getChild(1), function);
    }
  }

 public:
  void RunFmm(FmmTree &fmmTree, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    init(orderOfExpansion);
    auto &root = fmmTree.getRoot();
    upwardPass(root, [&](FmmTreeNode &node) {
      node.initCoefficients(orderOfExpansion);
      if (node.isLeaf()) {
        P2M(node, orderOfExpansion, container);
      }
    });
    upwardPass(root, [&](FmmTreeNode &node) {
      if (not node.isLeaf()) {
        M2M(node, orderOfExpansion);
      }
    });
    downwardPass(root, [&](FmmTreeNode &node) { M2L(node, orderOfExpansion); });
    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.getDepth() > 0) {
        L2L(node, orderOfExpansion);
      }
    });
    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.isLeaf()) {
        L2P(node, orderOfExpansion, container);
      }
    });
    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.isLeaf()) {
        for (auto p1 = container.getRegionIterator(node.getBoxMin(), node.getBoxMax()); p1.isValid(); ++p1) {
          for (auto nearCell : node.getNearFieldList()) {
            for (auto p2 = container.getRegionIterator(nearCell->getBoxMin(), nearCell->getBoxMax()); p2.isValid();
                 ++p2) {
              if (p1->getID() != p2->getID()) {
                NearField(*p1, *p2);
              }
            }
          }
        }
      }
    });
  }
};
}  // namespace autopas::fmm
