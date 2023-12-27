///**
//* @file LCC08CellHandler3B.h
//* @author N. Deng
//* @date 20.10.2023
// */
//
//#pragma once
//
//#include "autopas/containers/cellTraversals/CellTraversal.h"
//#include "autopas/baseFunctors/CellFunctor3B.h"
//#include "autopas/utils/ThreeDimensionalMapping.h"
//
//#include <chrono>
//#include <array>
//namespace autopas {
//
///**
//* This class provides the base for traversals using the c08 base step.
//*
//* The base step processBaseCell() computes one set of triwise interactions
//* between three cells for each spatial direction based on the baseIndex.
//* After executing the base step on all cells all triwise interactions for
//* all cells are done.
//*
//* @tparam ParticleCell the type of cells
//* @tparam PairwiseFunctor The functor that defines the interaction of three particles.
//* @tparam useSoA
//* @tparam useNewton3
// */
//template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
//class LCC08CellHandler3B {
// public:
//  /**
//  * Constructor of the LCC08CellHandler3B.
//  * @param functor The functor that defines the interaction of two particles.
//  * @param cellsPerDimension The number of cells per dimension.
//  * @param interactionLength Interaction length (cutoff + skin).
//  * @param cellLength cell length.
//  * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
//  * in that case the interactionLength is needed!
//   */
//  explicit LCC08CellHandler3B(Functor *functor, const std::array<unsigned long, 3> &cellsPerDimension,
//                              const double interactionLength, const std::array<double, 3> &cellLength,
//                              const std::array<unsigned long, 3> &overlap)
//      : _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
//        _cellOffsets{},
//        _interactionLength(interactionLength),
//        _cellLength(cellLength),
//        _overlap(overlap) {
//    computeOffsets(cellsPerDimension);
//  }
//
//  /**
//   * Computes all interactions
//   * @param cells vector of all cells.
//   * @param x X-index of base cell.
//   * @param y Y-index of base cell.
//   * @param z Z-index of base cell.
//   */
//  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);
//
//
// /**
//  * @copydoc autopas::CellTraversal::setSortingThreshold()
//  */
//  void setSortingThreshold(size_t sortingThreshold) { _cellFunctor.setSortingThreshold(sortingThreshold); }
//
//protected:
//
//  /**
//    * Combination of triplets for processBaseCell().
//   */
//  std::vector<std::tuple<long, long, long, std::array<double, 3>>> _cellOffsets;
//
//  /**
//   * Computes triplets used in processBaseCell()
//   */
//  void computeOffsets(const std::array<unsigned long, 3> &cellsPerDimension);
//
//
//  /**
//  * Overlap of interacting cells. Array allows asymmetric cell sizes.
//   */
//  const std::array<unsigned long, 3> _overlap;
//
// private:
//  /**
//   * CellFunctor to be used for the traversal defining the interaction between three cells.
//   */
//  internal::CellFunctor3B<typename ParticleCell::ParticleType, ParticleCell, Functor, dataLayout, useNewton3, true>
//      _cellFunctor;
//
//  Functor *_functor;
//
//  /**
//  * Interaction length (cutoff + skin).
//   */
//  const double _interactionLength;
//
//  /**
//  * Cell length in CellBlock3D.
//   */
//  const std::array<double, 3> _cellLength;
//};
//
//
//template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
//inline void LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3>::processBaseCell(
//    std::vector<ParticleCell> &cells, unsigned long baseIndex) {
//  for (auto const &[offset1, offset2, offset3, r] : _cellOffsets) {
//    const unsigned long index1 = baseIndex + offset1;
//    const unsigned long index2 = baseIndex + offset2;
//    const unsigned long index3 = baseIndex + offset3;
//
//    ParticleCell &cell1 = cells[index1];
//    ParticleCell &cell2 = cells[index2];
//    ParticleCell &cell3 = cells[index3];
//
//    if (index1 == index2 && index1 == index3 && index2 == index3) {
//      this->_cellFunctor.processCell(cell1);
//    } else if (index1 == index2 && index1 != index3) {
//      this->_cellFunctor.processCellPair(cell1, cell3);
//    } else if (index1 != index2 && index1 == index3) {
//      this->_cellFunctor.processCellPair(cell1, cell2);
//    } else if (index1 != index2 && index2 == index3) {
//      this->_cellFunctor.processCellPair(cell1, cell2);
//    } else {
//      this->_cellFunctor.processCellTriple(cell1, cell2, cell3);
//    }
//  }
//}
//
//template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
//inline void LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3>::computeOffsets(
//    const std::array<unsigned long, 3> &cellsPerDimension) {
//  using namespace utils::ArrayMath::literals;
//
//  static std::chrono::duration<double> accumulatedDuration = std::chrono::duration<double>::zero();
//
//  // Start the timer
//  auto startTime = std::chrono::high_resolution_clock::now();
//
//  // Helper function to get minimal distance between two cells
//  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
//    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * this->_cellLength[0],
//                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * this->_cellLength[1],
//                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * this->_cellLength[2]};
//  };
//
//  auto is_valid_distance = [&] (long x1, long y1, long z1, long x2, long y2, long z2, auto interactionlengthsq) {
//    auto dist = cellDistance(x1, y1, z1, x2, y2, z2);
//    return utils::ArrayMath::dot(dist, dist) > interactionlengthsq;
//  };
//
//  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);
//  _cellOffsets.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});
//  // offsets for the first cell
//  for (long x1 = 0; x1 <= static_cast<long>(this->_overlap[0]); ++x1) {
//    for (long y1 = 0; y1 <= static_cast<long>(this->_overlap[1]); ++y1) {
//      for (long z1 = 0; z1 <= static_cast<long>(this->_overlap[2]); ++z1) {
//
//        // offsets for the second cell
//        for (long x2 = 0; x2 <= static_cast<long>(this->_overlap[0]); ++x2) {
//          for (long y2 = 0; y2 <= static_cast<long>(this->_overlap[1]); ++y2) {
//            for (long z2 = 0; z2 <= static_cast<long>(this->_overlap[2]); ++z2) {
//              // check distance between cell 1 and cell 2
//              const auto dist12 = cellDistance(x1, y1, z1, x2, y2, z2);
//              const double dist12Square = utils::ArrayMath::dot(dist12, dist12);
//              if (dist12Square > interactionLengthSquare) continue;
//
//              // offsets for the third cell
//              for (long x3 = 0; x3 <= static_cast<long>(this->_overlap[0]); ++x3) {
//                for (long y3 = 0; y3 <= static_cast<long>(this->_overlap[1]); ++y3) {
//                  for (long z3 = 0; z3 <= static_cast<long>(this->_overlap[2]); ++z3) {
//                    // check distance between cell 1 and cell 3
//                    const auto dist13 = cellDistance(x1, y1, z1, x3, y3, z3);
//                    const double dist13Square = utils::ArrayMath::dot(dist13, dist13);
//                    if (dist13Square > interactionLengthSquare) continue;
//
//                    // check distance between cell 2 and cell 3
//                    const auto dist23 = cellDistance(x2, y2, z2, x3, y3, z3);
//                    const double dist23Squared = utils::ArrayMath::dot(dist23, dist23);
//                    if (dist23Squared > interactionLengthSquare) continue;
//
//                    const long offset1 = utils::ThreeDimensionalMapping::threeToOneD(
//                        x1, y1, z1, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
//
//                    const long offset2 = utils::ThreeDimensionalMapping::threeToOneD(
//                        x2, y2, z2, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
//
//                    const long offset3 = utils::ThreeDimensionalMapping::threeToOneD(
//                        x3, y3, z3, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
//
//                    if (offset1 > offset2 or offset2 >= offset3) continue;
//
//                    if ((x1 == 0 or x2 == 0 or x3 == 0) and(y1 == 0 or y2 == 0 or y3 == 0) and (z1 == 0 or z2 == 0 or z3 == 0)) {
//
//                      const std::array<double, 3> sortDirection = {(x1 + x2) * this->_cellLength[0],
//                                                                   (y1 + y2) * this->_cellLength[1],
//                                                                   (z1 + z2) * this->_cellLength[2]};
//                      _cellOffsets.emplace_back(offset1, offset2, offset3, utils::ArrayMath::normalize(sortDirection));
//                    }
//                  }
//                }
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//  accumulatedDuration += std::chrono::high_resolution_clock::now() - startTime;
//
//  // If needed, you can print the accumulated time after each call
//  std::cout << "Accumulated execution time in computeOffsets: " << accumulatedDuration.count() << " seconds." << std::endl;
//
//
//}
//
//}  // namespace autopas
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
  long ovX = static_cast<long>(this->_overlap[0]);
  long ovY = static_cast<long>(this->_overlap[1]);
  long ovZ = static_cast<long>(this->_overlap[2]);

  static std::chrono::duration<double> accumulatedDuration = std::chrono::duration<double>::zero();

  // Start the timer
  auto startTime = std::chrono::high_resolution_clock::now();

  // Helper function to get minimal distance between two cells
  auto cellDistance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    return std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * this->_cellLength[0],
                                 std::max(0l, (std::abs(y1 - y2) - 1l)) * this->_cellLength[1],
                                 std::max(0l, (std::abs(z1 - z2) - 1l)) * this->_cellLength[2]};
  };

  auto pairToTriplet = [](const std::vector<long>& pair, std::vector<std::tuple<long, long, long>>& offsets) {
    long firstOffset = pair[0];
    long secondOffset = pair[1];

    if (firstOffset > secondOffset) {
      std::swap(firstOffset, secondOffset);
    }
    offsets.push_back({firstOffset, firstOffset, secondOffset});
  };


  auto appendPlaneOffsets = [](long offset1, long offset2, const std::vector<std::vector<long>>& plane, std::vector<std::tuple<long, long, long>>& offsets) {
    for (long x = 1; x < plane.size(); ++x) {
      for (long y = 1; y < plane[x].size(); ++y) {
        offsets.push_back({offset1, offset2, plane[x][y]});
      }
    }
  };

  auto offsets3Edges2Same = [](const std::vector<long>& edgeFirst, const std::vector<long>& edgeSecond, std::vector<std::tuple<long, long, long>>& offsets) {
    for (long indexFirst = 1; indexFirst < edgeFirst.size(); ++indexFirst) {
      for (long indexSecond = 1; indexSecond < edgeSecond.size(); ++indexSecond) {
        for (long indexFirstNext = indexFirst + 1; indexFirstNext < edgeFirst.size(); ++indexFirstNext) {
          offsets.push_back({edgeFirst[indexFirst], edgeFirst[indexFirstNext], edgeSecond[indexSecond]});
        }
      }
    }
  };

  auto appendInnerOffsets = [](long offset1, long offset2, const std::vector<std::vector<std::vector<long>>>& allCells, std::vector<std::tuple<long, long, long>>& offsets) {
    for (long x = 1; x < allCells.size(); ++x) {
      for (long y = 1; y < allCells[x].size(); ++y) {
        for (long z = 1; z < allCells[x][y].size(); ++z) {
          offsets.push_back({offset1, offset2, allCells[x][y][z]});
        }
      }
    }
  };
  auto offsetsTwoEdgesOneWildCard = [&pairToTriplet, &appendPlaneOffsets, &appendInnerOffsets, &offsets3Edges2Same](
                                        const std::vector<long>& edgeFirst,
                                        const std::vector<long>& edgeSecond,
                                        const std::vector<std::vector<std::vector<long>>>& planes,
                                        const std::vector<std::vector<std::vector<long>>>& allCells,
                                        std::vector<std::tuple<long, long, long>>& offsets
                                    ) {
    std::size_t edgeFirstLength = edgeFirst.size();
    std::size_t edgeSecondLength = edgeSecond.size();

    for (long indexFirst = 1; indexFirst < edgeFirstLength; ++indexFirst) {
      for (long indexSecond = 1; indexSecond < edgeSecondLength; ++indexSecond) {
        pairToTriplet({edgeFirst[indexFirst], edgeSecond[indexSecond]}, offsets);

        for (const auto & plane : planes) {
          appendPlaneOffsets(edgeFirst[indexFirst], edgeSecond[indexSecond], plane, offsets);
        }

        appendInnerOffsets(edgeFirst[indexFirst], edgeSecond[indexSecond], allCells, offsets);
      }
    }

    offsets3Edges2Same(edgeFirst, edgeSecond, offsets);
    offsets3Edges2Same(edgeSecond, edgeFirst, offsets);
  };

  auto offsets1Edge2Wildcards = [&pairToTriplet, &appendPlaneOffsets, &appendInnerOffsets](
                                    const std::vector<long>& edge,
                                    const std::vector<std::vector<long>>& plane1,
                                    const std::vector<std::vector<long>>& plane2,
                                    const std::vector<std::vector<long>>& plane3,
                                    const std::vector<std::vector<std::vector<long>>>& allCells,
                                    std::vector<std::tuple<long, long, long>>& offsets
                                ) {
    for (long indexEdge = 1; indexEdge < edge.size(); ++indexEdge) {
      for (long x3 = 1; x3 < plane3.size(); ++x3) {
        for (long y3 = 1; y3 < plane3[x3].size(); ++y3) {
          pairToTriplet({edge[indexEdge], plane3[x3][y3]}, offsets);
          appendPlaneOffsets(edge[indexEdge], plane3[x3][y3], plane1, offsets);
          appendPlaneOffsets(edge[indexEdge], plane3[x3][y3], plane2, offsets);

          for (long x32 = 1; x32 < plane3.size(); ++x32) {
            for (long y32 = 1; y32 < plane3[0].size(); ++y32) {
              if (plane3[x3][y3] < plane3[x32][y32]) {
                offsets.push_back({edge[indexEdge], plane3[x3][y3], plane3[x32][y32]});
              }
            }
          }

          appendInnerOffsets(edge[indexEdge], plane3[x3][y3], allCells, offsets);

          for (long indexEdgeNext = indexEdge + 1; indexEdgeNext < edge.size(); ++indexEdgeNext) {
            offsets.push_back({edge[indexEdge], edge[indexEdgeNext], plane3[x3][y3]});
          }
        }
      }
    }
  };

  std::vector<std::vector<std::vector<long>>> all(ovX, std::vector<std::vector<long>>(ovY, std::vector<long>(ovZ)));
  std::vector<std::vector<long>> planeXy(ovX, std::vector<long>(ovY));
  std::vector<std::vector<long>> planeXz(ovX, std::vector<long>(ovZ));
  std::vector<std::vector<long>> planeYz(ovY, std::vector<long>(ovZ));
  std::vector<long> edgeX(ovX), edgeY(ovY), edgeZ(ovZ);

  // Initializing 3D vector 'all'
  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      for (long z = 0; z < ovZ; ++z) {
        all[x][y][z] = (ovX * ovX * x) + (y * ovY) + z;
      }
    }
  }
  // Initializing 2D vectors for planes
  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      planeXy[x][y] = utils::ThreeDimensionalMapping::threeToOneD(
          x, y, 0L, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
    }
  }
  for (long x = 0; x < ovX; ++x) {
    for (long z = 0; z < ovZ; ++z) {
      planeXz[x][z] = utils::ThreeDimensionalMapping::threeToOneD(
          x, 0L, z, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
    }
  }
  for (long y = 0; y < ovY; ++y) {
    for (long z = 0; z < ovZ; ++z) {
      planeXy[y][z] = utils::ThreeDimensionalMapping::threeToOneD(
          0L, y, z, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
    }
  }

  // Initializing 1D vectors for edges
  for (long x = 0; x < ovX; ++x) {
    edgeX[x] = utils::ThreeDimensionalMapping::threeToOneD(
        x, 0L, 0L, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
  }
  for (long y = 0; y < ovY; ++y) {
    edgeY[y] = utils::ThreeDimensionalMapping::threeToOneD(
        0L, y, 0L, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
  }
  for (long z = 0; z < ovZ; ++z) {
    edgeZ[z] = utils::ThreeDimensionalMapping::threeToOneD(
        0L, 0L, z, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
  }
  std::vector<std::vector<std::vector<long>>> planes = {planeXy, planeXz, planeYz};

  // Initialize offsets with base cell combos
  std::vector<std::tuple<long, long, long>> offsets = {{0, 0, 0}};

  long totalElements = ovX * ovY * ovZ;
  for (long a = 0; a < totalElements; ++a) {
    for (long b = a + 1; b < totalElements; ++b) {
      offsets.push_back({0, a, b});
    }
  }
  // Combining edges and planes
  offsetsTwoEdgesOneWildCard(edgeX, edgeY, planes, all, offsets);
  offsetsTwoEdgesOneWildCard(edgeX, edgeZ, planes, all, offsets);
  offsetsTwoEdgesOneWildCard(edgeY, edgeZ, planes, all, offsets);

  // Combining an edge with two planes
  offsets1Edge2Wildcards(edgeX, planeXy, planeXz, planeYz, all, offsets);
  offsets1Edge2Wildcards(edgeY, planeXy, planeYz, planeXz, all, offsets);
  offsets1Edge2Wildcards(edgeZ, planeYz, planeXz, planeXy, all, offsets);

  for (auto const &[offset1, offset2, offset3] : offsets) {
    _cellOffsets.emplace_back(offset1, offset2, offset3, std::array<double, 3>{1., 1., 1.});
  }

  accumulatedDuration += std::chrono::high_resolution_clock::now() - startTime;

  // If needed, you can print the accumulated time after each call
  std::cout << "Accumulated execution time in computeOffsets: " << accumulatedDuration.count() << " seconds. Size : " << _cellOffsets.size() << std::endl;


}

}  // namespace autopas
