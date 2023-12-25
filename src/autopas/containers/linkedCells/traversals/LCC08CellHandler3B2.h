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
  void offsetsTwoEdgesOneWildCard(
      const std::vector<long>& edge1,
      const std::vector<long>& edge2,
      const std::vector<std::vector<std::vector<long>>>& planes,
      const std::vector<std::vector<std::vector<long>>>& all,
      const std::function<std::vector<long>(const std::vector<long>&)>& pair_to_triple,
      const std::function<void(long, long, const std::vector<std::vector<long>>&, std::vector<std::vector<long>>&)> &append_plane_offsets,
      const std::function<void(long, long, const std::vector<std::vector<std::vector<long>>>&, std::vector<std::vector<long>>&)> &append_inner_offsets,
      const std::function<void(const std::vector<long>&, const std::vector<long>&, std::vector<std::vector<long>>&)> &offsets_3_edges_2_same
  ) {
    long ov_1 = edge1.size();
    long ov_2 = edge2.size();

    for (long a = 1; a < ov_1; ++a) {
      for (long b = 1; b < ov_2; ++b) {
        _cellOffsets.emplace_back(pair_to_triple({edge1[a], edge2[b]}));

        for (const auto& plane : planes) {
          append_plane_offsets(edge1[a], edge2[b], plane, _cellOffsets);
        }

        append_inner_offsets(edge1[a], edge2[b], all, _cellOffsets);
      }
    }

    offsets_3_edges_2_same(edge1, edge2, _cellOffsets);
    offsets_3_edges_2_same(edge2, edge1, _cellOffsets);
  }

  void offsetsOneEdgeTwoUnknown(std::vector<long> vector1, std::vector<long> vector2, std::vector<long> vector3,
                                std::vector<long> vector4, std::vector < std::vector<std::vector<long>>) {}
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
  auto pair_to_offsets = [](const std::vector<long>& pair, std::vector<std::vector<long>>& offsets) {
    long first = pair[0];
    long second = pair[1];
    if (first > second) {
      std::swap(first, second);
    }
    offsets.push_back({first, first, second});
  };

  auto append_inner_offsets = [](long o1, long o2, const std::vector<std::vector<std::vector<long>>>& all, std::vector<std::vector<long>>& offsets) {
    for (long x = 1; x < all.size(); ++x) {
      for (long y = 1; y < all[x].size(); ++y) {
        for (long z = 1; z < all[x][y].size(); ++z) {
          offsets.push_back({o1, o2, all[x][y][z]});
        }
      }
    }
  };

  auto append_plane_offsets = [](long o1, long o2, const std::vector<std::vector<long>>& plane, std::vector<std::vector<long>>& offsets) {
    for (long x = 1; x < plane.size(); ++x) {
      for (long y = 1; y < plane[x].size(); ++y) {
        offsets.push_back({o1, o2, plane[x][y]});
      }
    }
  };
  auto offsets_3_edges_2_same = [](const std::vector<long>& edge1, const std::vector<long>& edge2, std::vector<std::vector<long>>& offsets) {
    for (long a = 1; a < edge1.size(); ++a) {
      for (long b = 1; b < edge2.size(); ++b) {
        for (long a1 = a + 1; a1 < edge1.size(); ++a1) {
          offsets.push_back({edge1[a], edge1[a1], edge2[b]});
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
  std::vector<std::vector<long>> offsets = {{0, 0, 0}};
  // Fill offsets with combinations and other patterns

  // Combining edges and planes
  offsetsTwoEdgesOneWildCard(edgeX, edgeY, planes, all, offsets);
  offsetsTwoEdgesOneWildCard(edgeX, edgeZ, planes, all, offsets);
  offsetsTwoEdgesOneWildCard(edgeY, edgeZ, planes, all, offsets);

  // Combining an edge with two planes
  offsetsOneEdgeTwoUnknown(edgeX, planeXy, planeXz, planeYz, all, offsets);
  offsetsOneEdgeTwoUnknown(edgeY, planeXy, planeYz, planeXz, all, offsets);
  offsetsOneEdgeTwoUnknown(edgeZ, planeYz, planeXz, planeXy, all, offsets);


  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);
  _cellOffsets.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});

  accumulatedDuration += std::chrono::high_resolution_clock::now() - startTime;

  // If needed, you can print the accumulated time after each call
  std::cout << "Accumulated execution time in computeOffsets: " << accumulatedDuration.count() << " seconds." << std::endl;


}

}  // namespace autopas
