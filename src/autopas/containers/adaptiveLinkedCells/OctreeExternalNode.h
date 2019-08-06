/**
 * @file OctreeExternalNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <bitset>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include "autopas/containers/adaptiveLinkedCells/OctreeNode.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/inBox.h"

namespace autopas {
namespace internal {

template <class Particle, class ParticleCell>
class OctreeInternalNode;

/**
 * Class representing an external node (leaf) in an octree.
 */
template <class Particle, class ParticleCell>
class OctreeExternalNode : public OctreeNode<Particle, ParticleCell> {
 public:
  OctreeExternalNode() = default;

  /**
   * Create external node by combining the children of an internal node.
   * @param parent Parent node (nullptr for root).
   * @param cells
   * @param index Start index of this node.
   * @param nodes
   * @param boxMin
   * @param boxMax
   * @param level
   */
  OctreeExternalNode(OctreeNode<Particle, ParticleCell> *parent, std::vector<ParticleCell> &cells, unsigned long index,
                     std::array<OctreeNode<Particle, ParticleCell> *, 8> nodes, const std::array<double, 3> &boxMin,
                     const std::array<double, 3> &boxMax)
      : OctreeNode<Particle, ParticleCell>(parent, index, ArrayMath::mulScalar(ArrayMath::sub(boxMax, boxMin), 0.5)),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _cell(&cells.at(index)) {
    for (auto *node : nodes) {
      const size_t index = node->getIndex();
      // Move all particles to base cell
      for (auto it = cells[index].begin(); it.isValid(); ++it) {
        cells[index].addParticle(*it);
      }
      cells[index].clear();
    }
    cells[index].setCellLength(ArrayMath::sub(boxMax, boxMin));
  }

  /**
   * Create external node from particles in cell at position index.
   * @param parent Parent node (nullptr for root).
   * @param cells
   * @param index
   * @param boxMin
   * @param boxMax
   * @param level
   */
  OctreeExternalNode(OctreeNode<Particle, ParticleCell> *parent, std::vector<ParticleCell> &cells, unsigned long index,
                     const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax)
      : OctreeNode<Particle, ParticleCell>(
            parent, index, ArrayMath::add(boxMin, ArrayMath::mulScalar(ArrayMath::sub(boxMax, boxMin), 0.5))),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _cell(&cells.at(index)) {
    cells.at(index).setCellLength(ArrayMath::sub(boxMax, boxMin));
    // add all elements from content to underlying container
  };

  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const override;

  ParticleCell &getCell() const;

  size_t getSize() const override;

  OctreeNode<Particle, ParticleCell> *update(std::vector<ParticleCell> &cells) override;

  void apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func, ExecutionPolicy policy) override;

  operator std::string() const override {
    std::set<int> ids;
    for (auto e : _neighbors) ids.insert(e.first);
    return std::string(OctreeNode<Particle, ParticleCell>::_maxExp - this->_level, '\t') +
           std::to_string(this->_level) + ": [" + ArrayUtils::to_string(this->_center) + "] " +
           std::to_string(this->_index) + " NumPar: " + std::to_string(_cell->numParticles()) + " [" +
           ArrayUtils::to_string(ids) + "]\n";
  }

  bool isUpdateNeeded() const override;

  /**
   * Returns all neighbors of this Octree node.
   * @return Neighbors of this Octree node.
   */
  const std::set<std::pair<unsigned long, std::array<double, 3>>> &getNeighbors() const { return _neighbors; }

  std::set<OctreeExternalNode<Particle, ParticleCell> *> findSmallerNeighborsDirection(
      int xDir, int yDir, int zDir, OctreeInternalNode<Particle, ParticleCell> *gEq) {
    std::set<OctreeExternalNode<Particle, ParticleCell> *> neighbors;
    std::set<OctreeNode<Particle, ParticleCell> *> candidates;

    candidates.insert(gEq);
    while (not candidates.empty()) {
      auto begin = candidates.begin();
      if ((*begin)->isLeaf()) {
        neighbors.insert(dynamic_cast<OctreeExternalNode<Particle, ParticleCell> *>(*begin));
      } else {
        auto relPos{
            dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(this->_parent)->relPosOfChild(this->_center)};
        std::array<int, 3> relPosArray({abs(relPos[0] + xDir), abs(relPos[1] + yDir), abs(relPos[2] + zDir)});

        std::bitset<3> neighborPos;
        for (int i = 0; i < 3; i++) {
          if (relPosArray[i] % 2 == 1) {
            neighborPos.set(i);
          }
        }
        std::array<int, 3> xyz({xDir, yDir, zDir});
        auto *internalNode = dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(*begin);
        candidates.insert(internalNode->_children[neighborPos.to_ulong()]);
        for (int i = 0; i < 3; i++) {
          if (xyz[i] == 0) {
            neighborPos.flip(i);
            candidates.insert(internalNode->_children[neighborPos.to_ulong()]);
            for (int j = i + 1; j < 3; j++) {
              if (xyz[j] == 0) {
                neighborPos.flip(j);
                candidates.insert(internalNode->_children[neighborPos.to_ulong()]);
                neighborPos.flip(i);
                candidates.insert(internalNode->_children[neighborPos.to_ulong()]);
                break;
              }
            }
            break;
          }
        }
      }
      candidates.erase(begin);
    }

    return neighbors;
  }

  std::optional<OctreeNode<Particle, ParticleCell> *> findGreaterEqualNeighborDirection(int xDir, int yDir, int zDir) {
    // code similar to https://geidav.wordpress.com/2017/12/02/advanced-octrees-4-finding-neighbor-nodes/
    // root has neighbors
    if (this->_parent == nullptr) return {};

    // get relative position of this node to parent center
    auto relPos{
        dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(this->_parent)->relPosOfChild(this->_center)};

    std::array<int, 3> relPosArray({relPos[0] + xDir, relPos[1] + yDir, relPos[2] + zDir});

    if (std::any_of(relPosArray.cbegin(), relPosArray.cend(), [](auto e) { return e < 0 or e > 1; })) {
      if (xDir != 0 and relPos[0] >= 0 and relPos[0] < 2) xDir = 0;
      if (yDir != 0 and relPos[1] >= 0 and relPos[1] < 2) yDir = 0;
      if (zDir != 0 and relPos[2] >= 0 and relPos[2] < 2) zDir = 0;
      auto node = dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(this->_parent)
                      ->findGreaterEqualNeighborDirection(xDir, yDir, zDir);
      if (!node) {
        return {};
      }
      if ((*node)->isLeaf()) {
        return *node;
      }
      std::bitset<3> neighborPos;
      for (int i = 0; i < 3; i++) {
        if (relPosArray[i] % 2 == 1) {
          neighborPos.set(i);
        }
      }
      return std::make_optional(
          dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(*node)->_children[neighborPos.to_ulong()]);
    } else {
      std::bitset<3> neighborPos;
      for (int i = 0; i < 3; i++) {
        if (relPosArray[i] == 1) {
          neighborPos.set(i);
        }
      }

      return std::make_optional(
          dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(this->_parent)->_children[neighborPos.to_ulong()]);
    }
  }

  bool isLeaf() const override { return true; }
  /**
   * Updates all Neighbors. This is necessary after splitting or combination of nodes.
   */
  void updateNeigbors() override {
    _neighbors.clear();
    std::set<OctreeExternalNode<Particle, ParticleCell> *> result;
    std::cout << "Index:" << this->getIndex() << " Center: " << ArrayUtils::to_string(this->_center) << std::endl;
    for (int x = -1; x <= 1; x++) {
      for (int y = -1; y <= 1; y++) {
        for (int z = -1; z <= 1; z++) {
          if (x == 0 and y == 0 and z == 0) {
            continue;
          }
          auto n = findGreaterEqualNeighborDirection(x, y, z);
          if (n) {
            std::cout << "NI:" << (*n)->getIndex() << " Center: " << ArrayUtils::to_string((*n)->getCenter())
                      << std::endl;
            if ((*n)->isLeaf()) {
              result.insert(dynamic_cast<OctreeExternalNode<Particle, ParticleCell> *>(*n));
            } else {
              std::cout << "Not a leaf";
              auto neighbors = findSmallerNeighborsDirection(
                  x, y, z, dynamic_cast<OctreeInternalNode<Particle, ParticleCell> *>(*n));
              std::merge(result.begin(), result.end(), neighbors.begin(), neighbors.end(),
                         std::inserter(result, result.begin()));
            }
          }
        }
      }
    }

    for (auto e : result) {
      auto r = ArrayMath::normalize(ArrayMath::sub(this->_center, e->_center));
      _neighbors.insert(std::make_pair(e->getIndex(), r));
      // std::cout << "NI:" << e->getIndex() << std::endl;
    }
  }

  std::vector<Particle> getOutliers() override {
    std::vector<Particle> myInvalidParticles;
    // if empty
    if (not _cell->isNotEmpty()) {
      return myInvalidParticles;
    }

    for (auto &&pIter = _cell->begin(); pIter.isValid(); ++pIter) {
      // if not in cell
      if (utils::notInBox(pIter->getR(), _boxMin, _boxMax)) {
        myInvalidParticles.push_back(*pIter);
        pIter.deleteCurrentParticle();
      }
    }

    return myInvalidParticles;
  }
  /**
   * Sets max. number of elements inside of each node.
   * @param maxElements Max. number of elements.
   */
  static void setMaxElements(const unsigned long maxElements) { _maxElements = maxElements; }

  // private:
  const std::array<double, 3> _boxMin;
  const std::array<double, 3> _boxMax;
  ParticleCell *_cell;
  std::set<std::pair<unsigned long, std::array<double, 3>>> _neighbors;
  static unsigned long _maxElements;
};

template <class Particle, class ParticleCell>
unsigned long OctreeExternalNode<Particle, ParticleCell>::_maxElements = 0ul;

template <class Particle, class ParticleCell>
ParticleCell &OctreeExternalNode<Particle, ParticleCell>::getCell() const {
  return *_cell;
}

template <class Particle, class ParticleCell>
ParticleCell &OctreeExternalNode<Particle, ParticleCell>::getContainingCell(const std::array<double, 3> &) const {
  return getCell();
}

template <class Particle, class ParticleCell>
size_t OctreeExternalNode<Particle, ParticleCell>::getSize() const {
  return _cell->numParticles();
}

template <class Particle, class ParticleCell>
OctreeNode<Particle, ParticleCell> *OctreeExternalNode<Particle, ParticleCell>::update(
    std::vector<ParticleCell> &cells) {
  if (getSize() > _maxElements and this->_level != 0) {
    // split
    return new OctreeInternalNode<Particle, ParticleCell>(this->_parent, cells, this->getIndex(), _boxMin, _boxMax);
  } else {
    return this;
  }
}

template <class Particle, class ParticleCell>
bool OctreeExternalNode<Particle, ParticleCell>::isUpdateNeeded() const {
  // Check for outliers
  /*for (auto iter = this->_cells[this->_index].begin(); iter.isValid(); ++iter) {
    if (not utils::inBox(iter->getR(), _boxMin, _boxMax)) {
      return true;
    }
  }*/

  // Check whether node must be split
  return getSize() > _maxElements;
}

template <class Particle, class ParticleCell>
void OctreeExternalNode<Particle, ParticleCell>::apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func,
                                                       ExecutionPolicy) {
  func(*this);
}
}  // namespace internal
}  // namespace autopas