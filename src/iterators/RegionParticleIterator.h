/**
 * @file RegionParticleIterator.h specifies the RegionParticleIterator class
 * @author seckler
 * @date 03.04.2018
 */
#include <array>
#include <vector>
#include "ParticleIterator.h"

#pragma once

namespace autopas {
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
   * @param cont container of particle cells
   * @param startRegion lower corner of the region to iterate over
   * @param endRegion top corner of the region to iterate over
   */
  explicit RegionParticleIterator(std::vector<ParticleCell>* cont,
                                  std::array<double, 3> startRegion,
                                  std::array<double, 3> endRegion)
      : ParticleIterator<Particle, ParticleCell>(cont),
        _startRegion(startRegion),
        _endRegion(endRegion) {
    // ParticleIterator's constructor will initialize the Iterator, such that it
    // points to the first particle if one is found, otherwise the pointer is
    // not valid
    if (this->isValid()) {  // if there is NO particle, we can not dereference
                            // it, so we need a check.
      if (not(this->operator*()).inBox(_startRegion, _endRegion)) {
        operator++();
      }
    }
  }

  // todo: optimize! this currently just loops over all particles.
  // Better: implement smart version of next_non_empty_cell.
  /**
   * Increment operator of the iterator
   * @todo optimize! this version is currently very slow
   * @return
   */
  inline RegionParticleIterator<Particle, ParticleCell>& operator++() {
    do {
      ParticleIterator<Particle, ParticleCell>::operator++();
    } while (ParticleIterator<Particle, ParticleCell>::isValid() &&
             !(this->operator*()).inBox(_startRegion, _endRegion));
    return *this;
  }

 private:
  std::array<double, 3> _startRegion;
  std::array<double, 3> _endRegion;
};

}  // namespace autopas