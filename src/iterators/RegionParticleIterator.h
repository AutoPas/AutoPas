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
template <class Particle, class ParticleCell>
class RegionParticleIterator : public ParticleIterator<Particle, ParticleCell> {
 public:
  explicit RegionParticleIterator(std::vector<ParticleCell>* cont,
                                  std::array<double, 3> startRegion,
                                  std::array<double, 3> endRegion)
      : ParticleIterator<Particle, ParticleCell>(cont),
        _startRegion(startRegion),
        _endRegion(endRegion) {
    if(not (this->operator*()).inBox(_startRegion, _endRegion)){
      operator++();
    }
  }

  // TODO: optimize! this currently just loops over all particles.
  // Better: implement smart version of next_non_empty_cell.
  inline RegionParticleIterator<Particle, ParticleCell>& operator++() {
    do {
      ParticleIterator<Particle,ParticleCell>::operator++();
    } while (ParticleIterator<Particle,ParticleCell>::isValid() &&
             !(this->operator*()).inBox(_startRegion, _endRegion));
    return *this;
  }

 private:
  std::array<double, 3> _startRegion;
  std::array<double, 3> _endRegion;
};

}  // namespace autopas