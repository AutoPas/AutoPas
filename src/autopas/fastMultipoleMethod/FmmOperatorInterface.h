#pragma once

#include <iomanip>
#include <iostream>
#include "autopas/fastMultipoleMethod/FmmTree.h"
#include "autopas/fastMultipoleMethod/FmmTreeNode.h"
#include "autopas/utils/Timer.h"

namespace autopas::fmm {

template <class Particle, class ParticleCell>
class FmmOperatorInterface {
 private:
  virtual void init(long orderOfExpansion) = 0;
  virtual void p2m(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void m2m(FmmTreeNode &parent, long orderOfExpansion) = 0;
  virtual void m2l(FmmTreeNode &node, long orderOfExpansion) = 0;
  virtual void l2l(FmmTreeNode &child, long orderOfExpansion) = 0;
  virtual void l2p(FmmTreeNode &leaf, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) = 0;
  virtual void nearField(FmmParticle &p1, FmmParticle &p2) = 0;

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
  void printTime(const std::string &name, long time, long total) {
    std::cout << std::setfill(' ') << std::setw(4) << name << ":";
    std::cout << std::setfill(' ') << std::setw(9) << time;
    std::cout << "us | ";
    std::cout << std::setfill(' ') << std::setw(9) << std::fixed
              << (100 * static_cast<double>(time) / static_cast<double>(total));
    std::cout << "%" << std::endl;
  }

  void runFmm(FmmTree &fmmTree, long orderOfExpansion, AutoPas<Particle, ParticleCell> &container) {
    autopas::utils::Timer timer;
    long initTime, p2mTime, m2mTime, m2lTime, l2lTime, l2pTime, nearFieldTime;
    long totalTime;
    timer.start();

    init(orderOfExpansion);
    initTime = timer.stop();
    timer.start();

    auto &root = fmmTree.getRoot();

    upwardPass(root, [&](FmmTreeNode &node) {
      node.initCoefficients(orderOfExpansion);
      if (node.isLeaf()) {
        p2m(node, orderOfExpansion, container);
      }
    });
    p2mTime = timer.stop();
    timer.start();

    upwardPass(root, [&](FmmTreeNode &node) {
      if (not node.isLeaf()) {
        m2m(node, orderOfExpansion);
      }
    });
    m2mTime = timer.stop();
    timer.start();

    downwardPass(root, [&](FmmTreeNode &node) { m2l(node, orderOfExpansion); });
    m2lTime = timer.stop();
    timer.start();

    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.getDepth() > 0) {
        l2l(node, orderOfExpansion);
      }
    });
    l2lTime = timer.stop();
    timer.start();

    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.isLeaf()) {
        l2p(node, orderOfExpansion, container);
      }
    });
    l2pTime = timer.stop();
    timer.start();

    downwardPass(root, [&](FmmTreeNode &node) {
      if (node.isLeaf()) {
        for (auto p1 = container.getRegionIterator(node.getBoxMin(), node.getBoxMax()); p1.isValid(); ++p1) {
          for (auto nearCell : node.getNearFieldList()) {
            for (auto p2 = container.getRegionIterator(nearCell->getBoxMin(), nearCell->getBoxMax()); p2.isValid();
                 ++p2) {
              if (p1->getID() != p2->getID()) {
                nearField(*p1, *p2);
              }
            }
          }
        }
      }
    });
    nearFieldTime = timer.stop();

    totalTime = initTime + p2mTime + m2mTime + m2lTime + l2lTime + l2pTime + nearFieldTime;
    printTime("init", initTime, totalTime);
    printTime("p2m", p2mTime, totalTime);
    printTime("m2m", m2mTime, totalTime);
    printTime("m2l", m2lTime, totalTime);
    printTime("l2l", l2lTime, totalTime);
    printTime("l2p", l2pTime, totalTime);
    printTime("near", nearFieldTime, totalTime);
  }
};
}  // namespace autopas::fmm
