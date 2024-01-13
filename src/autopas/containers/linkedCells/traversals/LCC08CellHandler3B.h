/**
* @file LCC08CellHandler3B.h
* @author N. Deng
* @date 20.10.2023
 */

#pragma once

#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

#include <chrono>
#include <array>
namespace autopas {

/**
* This class provides the base for traversals using the c08 base step.
*
* The base step processBaseCell() computes one set of triwise interactions
* between three cells for each spatial direction based on the baseIndex.
* After executing the base step on all cells all triwise interactions for
* all cells are done.
*
* @tparam ParticleCell the type of cells
* @tparam PairwiseFunctor The functor that defines the interaction of three particles.
* @tparam useSoA
* @tparam useNewton3
 */
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC08CellHandler3B {
 public:
  /**
  * Constructor of the LCC08CellHandler3B.
  * @param functor The functor that defines the interaction of two particles.
  * @param cellsPerDimension The number of cells per dimension.
  * @param interactionLength Interaction length (cutoff + skin).
  * @param cellLength cell length.
  * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
  * in that case the interactionLength is needed!
   */
  explicit LCC08CellHandler3B(Functor *functor, const std::array<unsigned long, 3> &cellsPerDimension,
                              const double interactionLength, const std::array<double, 3> &cellLength,
                              const std::array<unsigned long, 3> &overlap)
      : _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
        _cellOffsets{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes all interactions
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);


 /**
  * @copydoc autopas::CellTraversal::setSortingThreshold()
  */
  void setSortingThreshold(size_t sortingThreshold) { _cellFunctor.setSortingThreshold(sortingThreshold); }

protected:

  /**
    * Combination of triplets for processBaseCell().
   */
  std::vector<std::tuple<long, long, long, std::array<double, 3>>> _cellOffsets;

  /**
   * Computes triplets used in processBaseCell()
   */
  void computeOffsets(const std::array<unsigned long, 3> &cellsPerDimension);


  /**
  * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<unsigned long, 3> _overlap;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between three cells.
   */
  internal::CellFunctor3B<typename ParticleCell::ParticleType, ParticleCell, Functor, dataLayout, useNewton3, true>
      _cellFunctor;

  Functor *_functor;

  /**
  * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
  * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;
};


template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long baseIndex) {
  for (auto const &[offset1, offset2, offset3, r] : _cellOffsets) {
    const unsigned long index1 = baseIndex + offset1;
    const unsigned long index2 = baseIndex + offset2;
    const unsigned long index3 = baseIndex + offset3;

    ParticleCell &cell1 = cells[index1];
    ParticleCell &cell2 = cells[index2];
    ParticleCell &cell3 = cells[index3];

    if (index1 == index2 && index1 == index3 && index2 == index3) {
      this->_cellFunctor.processCell(cell1);
    } else if (index1 == index2 && index1 != index3) {
      this->_cellFunctor.processCellPair(cell1, cell3);
    } else if (index1 != index2 && index1 == index3) {
      this->_cellFunctor.processCellPair(cell1, cell2);
    } else if (index1 != index2 && index2 == index3) {
      this->_cellFunctor.processCellPair(cell1, cell2);
    } else {
      this->_cellFunctor.processCellTriple(cell1, cell2, cell3);
    }
  }
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3>::computeOffsets(
    const std::array<unsigned long, 3> &cellsPerDimension) {
  using namespace utils::ArrayMath::literals;

  static std::chrono::duration<double> accumulatedDuration = std::chrono::duration<double>::zero();
  static long hitCounter = 0l;
  // Start the timer
  auto startTime = std::chrono::high_resolution_clock::now();

  // Helper function to get minimal distance between two cells
  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * _cellLength[0],
                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * _cellLength[1],
                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * _cellLength[2]};
  };

  auto is_valid_distance = [&] (long x1, long y1, long z1, long x2, long y2, long z2, auto interactionlengthsq) {
    auto dist = cellDistance(x1, y1, z1, x2, y2, z2);
    return utils::ArrayMath::dot(dist, dist) > interactionlengthsq;
  };
  std::cout << "Interaction Length : " << _interactionLength << ", Cell Length : " << _cellLength[0] << "  " <<
      _cellLength[1] << "  " << _cellLength[2];

  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);
  _cellOffsets.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});
  // offsets for the first cell
  for (long x1 = 0; x1 <= static_cast<long>(this->_overlap[0]); ++x1) {
    for (long y1 = 0; y1 <= static_cast<long>(this->_overlap[1]); ++y1) {
      for (long z1 = 0; z1 <= static_cast<long>(this->_overlap[2]); ++z1) {

        // offsets for the second cell
        for (long x2 = 0; x2 <= static_cast<long>(this->_overlap[0]); ++x2) {
          for (long y2 = 0; y2 <= static_cast<long>(this->_overlap[1]); ++y2) {
            for (long z2 = 0; z2 <= static_cast<long>(this->_overlap[2]); ++z2) {
              // check distance between cell 1 and cell 2
              const auto dist12 = cellDistance(x1, y1, z1, x2, y2, z2);
              const double dist12Square = utils::ArrayMath::dot(dist12, dist12);
              if (dist12Square > interactionLengthSquare) continue;

              // offsets for the third cell
              for (long x3 = 0; x3 <= static_cast<long>(this->_overlap[0]); ++x3) {
                for (long y3 = 0; y3 <= static_cast<long>(this->_overlap[1]); ++y3) {
                  for (long z3 = 0; z3 <= static_cast<long>(this->_overlap[2]); ++z3) {
                    // check distance between cell 1 and cell 3
                    const auto dist13 = cellDistance(x1, y1, z1, x3, y3, z3);
                    const double dist13Square = utils::ArrayMath::dot(dist13, dist13);
                    if (dist13Square > interactionLengthSquare) continue;

                    // check distance between cell 2 and cell 3
                    const auto dist23 = cellDistance(x2, y2, z2, x3, y3, z3);
                    const double dist23Squared = utils::ArrayMath::dot(dist23, dist23);
                    if (dist23Squared > interactionLengthSquare) continue;
                    const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
                        x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
                        x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    const long offset3 = utils::ThreeDimensionalMapping::threeToOneD(
                        x3, y3, z3, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));

                    if (offset1 > offset2 or offset2 >= offset3) continue;

                    if ((x1 == 0 or x2 == 0 or x3 == 0) and(y1 == 0 or y2 == 0 or y3 == 0) and (z1 == 0 or z2 == 0 or z3 == 0)) {
                      hitCounter++;

                      const std::array<double, 3> sortDirection = {(x1 + x2) * this->_cellLength[0],
                                                                   (y1 + y2) * this->_cellLength[1],
                                                                   (z1 + z2) * this->_cellLength[2]};
                      _cellOffsets.emplace_back(offset1, offset2, offset3, utils::ArrayMath::normalize(sortDirection));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "Hits for overlap " << _overlap[0] << " : " << hitCounter << std::endl;
  hitCounter = 1;
}

}  // namespace autopas