/**
 * @file RegionParticleIterator.h specifies the RegionParticleIterator class
 * @author seckler
 * @date 03.04.2018
 */

#pragma once

#include <autopas/utils/inBox.h>
#include <array>
#include <vector>
#include "autopas/iterators/ParticleIterator.h"

namespace autopas {
namespace internal {
/**
 * RegionParticleIterator to iterate over all particles within a specific region
 * @todo optimize the region particle iterater. Currently we iterate over all
 * particles
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 */
template <class Particle, class ParticleCell>
class RegionParticleIterator : public ParticleIterator<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the RegionParticleIterator.
   *
   * @param cont Container of particle cells.
   * @param cellsPerDimension Total number of cells along each dimension.
   * @param startRegion Lower corner of the region to iterate over.
   * @param endRegion Top corner of the region to iterate over.
   * @param startIndex Index of the first cell in the region of interest.
   * @param endIndex Index of the last cell in the region of interest.
   * @param flagManager The CellBorderAndFlagManager that shall be used to
   * query the cell types. Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells
   * shall be iterated through.
   */
  explicit RegionParticleIterator(std::vector<ParticleCell> *cont, const std::array<size_t, 3> cellsPerDimension,
                                  std::array<double, 3> startRegion, std::array<double, 3> endRegion, size_t startIndex,
                                  size_t endIndex, CellBorderAndFlagManager *flagManager = nullptr,
                                  IteratorBehavior behavior = haloAndOwned)
      : ParticleIterator<Particle, ParticleCell>(cont, startIndex, flagManager, behavior),
        _cellsPerDim(cellsPerDimension),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _startIndex(startIndex),
        _endIndex(endIndex) {
    // ParticleIterator's constructor will initialize the Iterator, such that it
    // points to the first particle if one is found, otherwise the pointer is
    // not valid
    if (this->isValid()) {  // if there is NO particle, we can not dereference it, so we need a check.
      if (notInBox(this->operator*().getR(), _startRegion, _endRegion)) {
        operator++();
      }
    }
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++
   * @todo optimize! this version is currently very slow
   */
  inline RegionParticleIterator<Particle, ParticleCell> &operator++() override {
    do {
      ParticleIterator<Particle, ParticleCell>::operator++();
    } while (ParticleIterator<Particle, ParticleCell>::isValid() &&
             notInBox(this->operator*().getR(), _startRegion, _endRegion));
    return *this;
  }

  // @todo add test of clone
  inline ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new RegionParticleIterator<Particle, ParticleCell>(*this);
  }

 private:
  std::array<size_t, 3> _cellsPerDim;
  std::array<double, 3> _startRegion;
  std::array<double, 3> _endRegion;
  size_t _startIndex;
  size_t _endIndex;
};
}  // namespace internal
}  // namespace autopas