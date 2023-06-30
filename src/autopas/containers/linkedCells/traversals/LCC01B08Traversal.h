/**
* @file LCC01B08Traversal.h
* @author Luis Gall
* @date 28.06.2023
*/

#pragma once

#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          bool combineSoA = false>

class LCC01B08Traversal : public LCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA> {

 public:

  explicit LCC01B08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                             const double interactionLength, const std::array<double, 3> &cellLength)
      : LCC01Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3, combineSoA>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _ownCellOffsets{}
  {
    computeOffsets();
  }

  // TODO : And only applicable for verlet list cells pair neighbor list rebuilding
  [[nodiscard]] bool isApplicable() const override {
    return dataLayout == DataLayoutOption::aos and not combineSoA;
  }

  // Offsets should be computed like in a c18 based traversal -> copied from LCC18Traversal.h
  void computeOffsets() override {
    _ownCellOffsets.resize(2 * this->_overlap[1] + 1, std::vector<offsetArray_t>(2 * this->_overlap[0] + 1));
    const std::array<long, 3> _overlap_s = utils::ArrayUtils::static_cast_array<long>(this->_overlap);

    const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);

    for (long z = 0l; z <= _overlap_s[2]; ++z) {
      for (long y = -_overlap_s[1]; y <= _overlap_s[1]; ++y) {
        for (long x = -_overlap_s[0]; x <= _overlap_s[0]; ++x) {
          const long offset = utils::ThreeDimensionalMapping::threeToOneD(
              x, y, z, utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension));

          if (offset < 0l) {
            continue;
          }
          // add to each applicable special case
          for (long yArray = -_overlap_s[1]; yArray <= _overlap_s[1]; ++yArray) {
            if (std::abs(yArray + y) <= _overlap_s[1]) {
              for (long xArray = -_overlap_s[0]; xArray <= _overlap_s[0]; ++xArray) {
                if (std::abs(xArray + x) <= _overlap_s[0]) {
                  std::array<double, 3> pos = {};
                  pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
                  pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
                  pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
                  // calculate distance between base cell and other cell
                  const double distSquare = utils::ArrayMath::dot(pos, pos);
                  // only add cell offset if cell is within cutoff radius
                  if (distSquare <= interactionLengthSquare) {
                    _ownCellOffsets[yArray + _overlap_s[1]][xArray + _overlap_s[0]].push_back(
                        std::make_pair(offset, utils::ArrayMath::normalize(pos)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

 private:
  using offsetArray_t = std::vector<std::pair<unsigned long, std::array<double, 3>>>;

  std::vector<std::vector<offsetArray_t>> _ownCellOffsets;

  unsigned long getIndex(
      const unsigned long pos, const unsigned int dim) const {
    unsigned long index;
    if (pos < this->_overlap[dim]) {
      index = pos;
    } else if (pos < this->_cellsPerDimension[dim] - this->_overlap[dim]) {
      index = this->_overlap[dim];
    } else {
      index = pos - this->_cellsPerDimension[dim] + 2 * this->_overlap[dim] + 1ul;
    }
    return index;
  }

  // copied from LCC18Traversal
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) override {
    const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

    const unsigned long xArray = getIndex(x, 0);
    const unsigned long yArray = getIndex(y, 1);

    ParticleCell &baseCell = cells[baseIndex];

    this->_cellFunctor.setOnlyDirty(this->_onlyDirty);

    offsetArray_t &offsets = _ownCellOffsets[yArray][xArray];
    for (auto const &[offset, r] : offsets) {
      unsigned long otherIndex = baseIndex + offset;
      ParticleCell &otherCell = cells[otherIndex];

      if (baseIndex == otherIndex) {
        this->_cellFunctor.processCell(baseCell);
      } else {
        this->_cellFunctor.processCellPair(baseCell, otherCell, r);
      }
    }
  }
};

}