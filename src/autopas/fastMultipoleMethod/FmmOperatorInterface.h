#pragma once

#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"

namespace autopas::fmm {

template <class Particle, class ParticleCell>
class FmmOperatorInterface {
 private:
  virtual void P2M(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void M2M(FmmTreeNode &parent, long orderOfExpansion) = 0;
  virtual void M2L(FmmTreeNode &node, long orderOfExpansion) = 0;
  virtual void L2L(FmmTreeNode &child, long orderOfExpansion) = 0;
  virtual void L2P(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void NearField(FmmParticle &p1, FmmParticle &p2) = 0;

  void P2MRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    node.initCoefficients(orderOfExpansion);
    if (node.isLeaf()) {
      P2M(node, orderOfExpansion, container);
    } else {
      P2MRec(node.getChild(0), orderOfExpansion, container);
      P2MRec(node.getChild(1), orderOfExpansion, container);
    }
  }

  void M2MRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    if (not node.isLeaf()) {
      M2MRec(node.getChild(0), orderOfExpansion, container);
      M2MRec(node.getChild(1), orderOfExpansion, container);
      M2M(node, orderOfExpansion);
    }
  };

  void M2LRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    M2L(node, orderOfExpansion);
    if (not node.isLeaf()) {
      M2LRec(node.getChild(0), orderOfExpansion, container);
      M2LRec(node.getChild(1), orderOfExpansion, container);
    }
  };

  void L2LRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    if (node.getDepth() > 0) {
      L2L(node, orderOfExpansion);
    }
    if (not node.isLeaf()) {
      L2LRec(node.getChild(0), orderOfExpansion, container);
      L2LRec(node.getChild(1), orderOfExpansion, container);
    }
  };

  void L2PRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    if (node.isLeaf()) {
      L2P(node, orderOfExpansion, container);
    } else {
      L2PRec(node.getChild(0), orderOfExpansion, container);
      L2PRec(node.getChild(1), orderOfExpansion, container);
    }
  }

  void NearFieldRec(FmmTreeNode &node, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
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
    } else {
      NearFieldRec(node.getChild(0), orderOfExpansion, container);
      NearFieldRec(node.getChild(1), orderOfExpansion, container);
    }
  }

 public:
  void RunFmm(FmmTree &fmmTree, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    std::cout << "RunFmm 1" << std::endl;
    P2MRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 2" << std::endl;
    M2MRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 3" << std::endl;
    M2LRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 4" << std::endl;
    L2LRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 5" << std::endl;
    L2PRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 6" << std::endl;
    NearFieldRec(fmmTree.getRoot(), orderOfExpansion, container);
    std::cout << "RunFmm 6" << std::endl;
  }
};
}  // namespace autopas::fmm
