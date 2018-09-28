/**
 * @file RegionParticleIterator.h specifies the RegionParticleIterator class
 * @author seckler
 * @date 03.04.2018
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/inBox.h"

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
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
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
    _startIndex3D = utils::ThreeDimensionalMapping::oneToThreeD(_startIndex, _cellsPerDim);
    _endIndex3D = utils::ThreeDimensionalMapping::oneToThreeD(_endIndex, _cellsPerDim);

    _shortJump = _cellsPerDim[0] - (_endIndex3D[0] - _startIndex3D[0]) - 1;
    _longJump = _cellsPerDim[0] * (_cellsPerDim[1] - (_endIndex3D[1] - _startIndex3D[1]) - 1);

    // ParticleIterator's constructor will initialize the Iterator, such that it
    // points to the first particle if one is found, otherwise the pointer is
    // not valid
    if (ParticleIterator<Particle, ParticleCell>::isValid()) {  // if there is NO particle, we can not dereference it,
                                                                // so we need a check.
      if (utils::notInBox(this->operator*().getR(), _startRegion, _endRegion)) {
        operator++();
      }
    }
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++
   */
  inline RegionParticleIterator<Particle, ParticleCell> &operator++() override {
    do {
      ParticleIterator<Particle, ParticleCell>::operator++();
    } while (ParticleIterator<Particle, ParticleCell>::isValid() &&
             utils::notInBox(this->operator*().getR(), _startRegion, _endRegion) && this->getCurrentCellId() <= _endIndex);
    return *this;
  }

  bool isValid() const override {
    return ParticleIterator<Particle, ParticleCell>::isValid() &&
           utils::inBox(this->operator*().getR(), _startRegion, _endRegion);
  }

  // @todo add test of clone
  inline ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new RegionParticleIterator<Particle, ParticleCell>(*this);
  }

 private:
  /**
   * @copydoc
   * @note overrides call in ParticleIterator<Particle, ParticleCell>::operator++()
   */
  void next_non_empty_cell() override {
    // find the next non-empty cell
    const int stride = autopas_get_num_threads();  // num threads
    for (this->_iteratorAcrossCells += stride; this->getCurrentCellId() <= _endIndex;
         this->_iteratorAcrossCells += stride) {
      auto posDim0 = this->getCurrentCellId() % _cellsPerDim[0];
      if (posDim0 > _endIndex3D[0] || posDim0 < _startIndex3D[0]) {
        this->_iteratorAcrossCells += _shortJump;
      }
      auto posDim1 = (this->getCurrentCellId() / _cellsPerDim[0]) % _cellsPerDim[1];
      if (posDim1 > _endIndex3D[1] || posDim1 < _startIndex3D[1]) {
        this->_iteratorAcrossCells += _longJump;
      }

      if (this->_iteratorAcrossCells < this->_vectorOfCells->end() and this->_iteratorAcrossCells->isNotEmpty() and
          this->isCellTypeBehaviorCorrect()) {
        this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
        break;
      }
    }
  }

  std::array<size_t, 3> _cellsPerDim;
  std::array<double, 3> _startRegion;
  std::array<double, 3> _endRegion;
  std::array<size_t, 3> _startIndex3D;
  std::array<size_t, 3> _endIndex3D;
  size_t _startIndex;
  size_t _endIndex;
  /**
   * When the iterator passes _endIndex3D[0] this is the jump distance to increase y by 1 and setting x to
   * _startIndex3D[0] (or behind this depending on the stride).
   */
  size_t _shortJump;
  /**
   * When the iterator passes _endIndex3D[1] this is the jump distance to increase z by 1 and setting y to
   * _startIndex3D[1] (or behind this depending on the stride).
   */
  size_t _longJump;
};
}  // namespace internal
}  // namespace autopas