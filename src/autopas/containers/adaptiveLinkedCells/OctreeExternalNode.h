/**
 * @file OctreeExternalNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>
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
  return getSize() > _maxElements;  //@Todo outliers
}

template <class Particle, class ParticleCell>
void OctreeExternalNode<Particle, ParticleCell>::apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func,
                                                       ExecutionPolicy) {
  func(*this);
}
}  // namespace internal
}  // namespace autopas