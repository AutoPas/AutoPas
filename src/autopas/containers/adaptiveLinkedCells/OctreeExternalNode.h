/**
 * @file OctreeExternalNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>
#include <optional>
#include <set>
#include "autopas/containers/adaptiveLinkedCells/OctreeNode.h"
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
      : OctreeNode<Particle, ParticleCell>(parent, index),
        _center(ArrayMath::mulScalar(ArrayMath::sub(boxMax, boxMin), 0.5)),
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
      : OctreeNode<Particle, ParticleCell>(parent, index),
        _center(ArrayMath::mulScalar(ArrayMath::sub(boxMax, boxMin), 0.5)),
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
    /// @todo replace by ArrayUtils::toString()
    return std::to_string(this->_level) + ":" + std::to_string(_center[0]) + std::to_string(_center[1]) +
           std::to_string(_center[2]);
  }

  bool isUpdateNeeded() const override;

  /**
   * Returns all neighbors of this Octree node.
   * @return Neighbors of this Octree node.
   */
  const std::set<std::pair<unsigned long, std::array<double, 3>>> &getNeighbors() const { return _neighbors; }

  /*std::set<OctreeNode<Particle, ParticleCell>> findSmallerNeighborsDirection(
      int xDir, int yDir, int zDir, std::set<OctreeNode<Particle, ParticleCell>> gEq) {
    std::set<OctreeNode<Particle, ParticleCell>> neighbors;

    if (gEq.isEmpty()) {
      return neighbors;
    }

    while (not gEq.isEmpty()) {
      auto begin = gEq.begin();
      if (begin->isLeaf()) {
        neighbors.insert(*begin);
      } else {
        auto relPos{this->_parent.relPosOfChild(this->_center)};
        std::array<int, 3> relPosArray({relPos[0] + xDir, relPos[1] + yDir, relPos[2] + zDir});

        std::bitset<3> neighborPos;
        for (int i = 0; i < 3; i++) {
          if (relPosArray[i] % 2 == 1) {
            neighborPos.set(i);
          }
        }
        begin->_children[neighborPos.to_ulong()];
      }
      gEq.erase(begin);
    }

    return neighbors;
  }

  std::optional<OctreeNode<Particle, ParticleCell>> findGreaterEqualNeighborDirection(int xDir, int yDir, int zDir) {
    // code similar to https://geidav.wordpress.com/2017/12/02/advanced-octrees-4-finding-neighbor-nodes/
    if (this->_parent == nullptr) return {};
    auto relPos{this->_parent.relPosOfChild(this->_center)};
    std::array<int, 3> relPosArray({relPos[0] + xDir, relPos[1] + yDir, relPos[2] + zDir});
    if (std::any_of(relPosArray.cbegin(), relPosArray.cend(), [](auto e) { return e < 0 or e > 1; })) {
      auto node = this->_parent.findGreaterEqualNeighborDirection(xDir, yDir, zDir);
      if (node or node->isLeaf()) {
        return {};
      }
      std::bitset<3> neighborPos;
      for (int i = 0; i < 3; i++) {
        if (relPosArray[i] % 2 == 1) {
          neighborPos.set(i);
        }
      }
      return std::make_optional(node->_children[neighborPos.to_ulong()]);
    } else {
      std::bitset<3> neighborPos;
      for (int i = 0; i < 3; i++) {
        if (relPosArray[i] == 1) {
          neighborPos.set(i);
        }
      }

      return std::make_optional(this->parent._children[neighborPos.to_ulong()]);
    }
  }
*/
  /**
   * Updates all Neighbors. This is necessary after splitting or combination of nodes.
   */
  /*void updateNeigbors() {
    _neighbors.clear();
    std::set<OctreeNode<Particle, ParticleCell>> result;
    for (int x = -1; x <= 1; x++) {
      for (int y = -1; y <= 1; y++) {
        for (int z = -1; z <= 1; z++) {
          if (x == 0 and y == 0 and z == 0) {
            continue;
          }
          auto n = findGreaterEqualNeighborDirection(x, y, z, result);
          if (n) {
            auto neighbors = findSmallerNeighborsDirection(x, y, z, *n);
            std::merge(result.begin(), result.end(), neighbors.begin(), neighbors.end(),
                       std::inserter(_neighbors, _neighbors.begin()))
          }
        }
      }
    }
  }*/
  /**
   * Sets max. number of elements inside of each node.
   * @param maxElements Max. number of elements.
   */
  static void setMaxElements(const unsigned long maxElements) { _maxElements = maxElements; }

 private:
  const std::array<double, 3> _center;
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
  if (getSize() > _maxElements and this->_level > 0) {
    // split
    return new OctreeInternalNode<Particle, ParticleCell>(this, cells, this->getIndex(), _boxMin, _boxMax);
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