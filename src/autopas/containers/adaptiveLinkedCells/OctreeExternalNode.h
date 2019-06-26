/**
 * @file OctreeExternalNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>
#include "autopas/containers/adaptiveLinkedCells/OctreeInternalNode.h"
#include "autopas/containers/adaptiveLinkedCells/OctreeNode.h"

namespace autopas {

/**
 * Class representing an external node (leaf) in an octree.
 */
template <class Particle, class ParticleCell>
class OctreeExternalNode : public OctreeNode<Particle, ParticleCell> {
 public:
  OctreeExternalNode() = default;

  OctreeExternalNode(std::array<std::unique_ptr<OctreeNode<Particle, ParticleCell>>, 8> content,
                     const std::array<double, 3> &center, const unsigned int level)
      : OctreeNode<Particle, ParticleCell>(level), _center(center) {
    // add all elements from content to underlying container
  }

  OctreeExternalNode(std::vector<ParticleCell> &cells, unsigned long index, const std::array<double, 3> &center,
                     const unsigned int level)
      : OctreeNode<Particle, ParticleCell>(level),
        _center(center),
        _cellView(cells.at(index), std::array<double, 3>{1, 0, 0}){
            // add all elements from content to underlying container
        };

  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const override;

  size_t getSize() const override;

  // std::unique_ptr<OctreeNode<ParticleCell>> update() override;

  bool isUpdateNeeded() const override;

  /**
   * Sets max. number of elements inside of each node.
   * @param maxElements Max. number of elements.
   */
  static void setMaxElements(const unsigned long maxElements) { _maxElements = maxElements; }

 private:
  const std::array<double, 3> _center;
  SortedCellView<Particle, ParticleCell> _cellView;
  static unsigned long _maxElements;
};

template <class Particle, class ParticleCell>
unsigned long OctreeExternalNode<Particle, ParticleCell>::_maxElements = 0ul;

template <class Particle, class ParticleCell>
ParticleCell &OctreeExternalNode<Particle, ParticleCell>::getContainingCell(const std::array<double, 3> &pos) const {
  return *(_cellView._cell);
}

template <class Particle, class ParticleCell>
size_t OctreeExternalNode<Particle, ParticleCell>::getSize() const {
  return _cellView._cell->numParticles();
}
/*
template <class Particle, class ParticleCell>
std::unique_ptr<OctreeNode<ParticleCell>> OctreeExternalNode<ParticleCell>::update() {
  if (getSize() > _maxElements) {
    // split
    return OctreeInternalNode<ParticleCell>(); // _cell, _center, this->_level + 1ul
}  // namespace autopas
else {
  return *this;
}
}
*/

template <class Particle, class ParticleCell>
bool OctreeExternalNode<Particle, ParticleCell>::isUpdateNeeded() const {
  return getSize() > _maxElements;
}

}  // namespace autopas