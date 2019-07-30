/**
 * @file OctreeInternalNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include "autopas/containers/adaptiveLinkedCells/OctreeNode.h"

#include <array>
#include <bitset>
#include <memory>

namespace autopas {
namespace internal {

template <class Particle, class ParticleCell>
class OctreeExternalNode;

/**
 * Class representing an internal node in an octree.
 */
template <class Particle, class ParticleCell>
class OctreeInternalNode : public OctreeNode<Particle, ParticleCell> {
 public:
  OctreeInternalNode() = default;

  OctreeInternalNode(std::vector<ParticleCell> &cells, unsigned long index, const std::array<double, 3> &center,
                     const unsigned int level);

  ~OctreeInternalNode() { delete[] _children.data(); }

  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const override;

  size_t getSize() const override;

  OctreeNode<Particle, ParticleCell> *update(std::vector<ParticleCell> &cells) override;

  void apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func, ExecutionPolicy policy) override;

  operator std::string() const override {
    return std::accumulate(_children.cbegin(), _children.cend(), std::string(),
                           [](auto &acc, auto elem) { return acc + static_cast<std::string>(*elem); });
  }

  bool isUpdateNeeded() const override;

  /**
   * Sets min. number of elements inside of each node.
   * @param minElements Min. number of elements.
   */
  static void setMinElements(const unsigned long minElements) { _minElements = minElements; }

 private:
  const std::array<double, 3> _center;
  std::array<OctreeNode<Particle, ParticleCell> *, 8> _children;
  static unsigned long _minElements;
};

template <class Particle, class ParticleCell>
unsigned long OctreeInternalNode<Particle, ParticleCell>::_minElements = 0ul;

template <class Particle, class ParticleCell>
OctreeInternalNode<Particle, ParticleCell>::OctreeInternalNode(std::vector<ParticleCell> &cells, unsigned long index,
                                                               const std::array<double, 3> &center,
                                                               const unsigned int level)
    : OctreeNode<Particle, ParticleCell>(level, index), _center(center) {
  for (unsigned int i = 0; i < 8; ++i) {
    std::bitset<3> index(i);
    std::array<double, 3> childCenter{_center};
    for (unsigned int d = 0; d < 3; ++d) {
      if (index.test(d)) {
        childCenter[d] += cells[index.to_ulong()].getCellLength()[d] / 2.0;
      } else {
        childCenter[d] -= cells[index.to_ulong()].getCellLength()[d] / 2.0;
      }
    }
    //_children = std::make_unique<OctreeExternalNode>(, childCenter, this->_level + 1)
  }
}

template <class Particle, class ParticleCell>
ParticleCell &OctreeInternalNode<Particle, ParticleCell>::getContainingCell(const std::array<double, 3> &pos) const {
  std::bitset<3> index;
  for (unsigned int d = 0; d < 3; ++d) {
    if (pos[d] > _center[d]) {
      index.set(d);
    }
  }

  return _children[index.to_ulong()]->getContainingCell(pos);
}

template <class Particle, class ParticleCell>
size_t OctreeInternalNode<Particle, ParticleCell>::getSize() const {
  return std::accumulate(_children.cbegin(), _children.cend(), 0ul,
                         [](auto acc, const auto e) { return acc + e->getSize(); });
}

template <class Particle, class ParticleCell>
OctreeNode<Particle, ParticleCell> *OctreeInternalNode<Particle, ParticleCell>::update(
    std::vector<ParticleCell> &cells) {
  if (getSize() < _minElements) {
    // combine
    return new OctreeExternalNode<Particle, ParticleCell>(cells, this->getIndex(), _children, _center, this->_level);
  } else {
    // call update on all children
    std::transform(_children.begin(), _children.end(), _children.begin(), [&](auto e) { return e->update(cells); });
    return this;
  }
}

template <class Particle, class ParticleCell>
bool OctreeInternalNode<Particle, ParticleCell>::isUpdateNeeded() const {
  if (getSize() < _minElements) {
    return true;
  } else {
    // call isUpdateNeeded on all children
    return std::any_of(_children.cbegin(), _children.cend(), [](const auto e) { return e->isUpdateNeeded(); });
  }
}

template <class Particle, class ParticleCell>
void OctreeInternalNode<Particle, ParticleCell>::apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func,
                                                       ExecutionPolicy policy) {
  switch (policy) {
    case ExecutionPolicy::seq: {
      for (int i = 0; i < 8; ++i) {
        _children[i]->apply(func, policy);
      }
    }
    case ExecutionPolicy::par: {
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#pragma omp taskloop if (this->_level > 1)
#endif
      for (int i = 0; i < 8; ++i) {
        _children[i]->apply(func, policy);
      }
    }
  }
}

}  // namespace internal
}  // namespace autopas