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
  long ovX = static_cast<long>(this->_overlap[0]) + 1l;
  long ovY = static_cast<long>(this->_overlap[1]) + 1l;
  long ovZ = static_cast<long>(this->_overlap[2]) + 1l;

  static std::chrono::duration<double> accumulatedDuration = std::chrono::duration<double>::zero();

  // Start the timer
  auto startTime = std::chrono::high_resolution_clock::now();
  const auto interactionLengthSquare(this->_interactionLength * this->_interactionLength);
  const auto cellLength1 = _cellLength[0];
  const auto cellLength2 = _cellLength[1];
  const auto cellLength3 = _cellLength[2];
  std::vector<std::vector<std::vector<long>>> cells(ovX, std::vector<std::vector<long>>(ovY, std::vector<long>(ovZ)));
  auto is_valid_distance = [&](long x1, long y1, long z1, long x2, long y2, long z2) {
    auto dist = std::array<double, 3>{std::max(0l, (std::abs(x1 - x2) - 1l)) * cellLength1,
                                      std::max(0l, (std::abs(y1 - y2) - 1l)) * cellLength2,
                                      std::max(0l, (std::abs(z1 - z2) - 1l)) * cellLength3};
    return utils::ArrayMath::dot(dist, dist) <= interactionLengthSquare;
  };

  auto emplaceOffset = [&](long x1, long y1, long z1, long x2, long y2, long z2, long x3, long y3, long z3) {
    if (is_valid_distance(x2, y2, z2, x3, y3, z3)) {
      if (is_valid_distance(x1, y1, z1, x3, y3, z3)) {
        std::array<double, 3> sortingDirection = {(x1 + x2) * cellLength1, (y1 + y2) * cellLength2,(z1 + z2) * cellLength3};
        _cellOffsets.emplace_back(cells[x1][y1][z1], cells[x2][y2][z2], cells[x3][y3][z3],
                                  utils::ArrayMath::normalize(sortingDirection));
      }
    }
  };


  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      for (long z = 0; z < ovZ; ++z) {
        cells[x][y][z] = utils::ThreeDimensionalMapping::threeToOneD(
            x, y, z, utils::ArrayUtils::static_cast_copy_array<long>(cellsPerDimension));
      }
    }
  }
  // base cell offsets
  _cellOffsets.emplace_back(0, 0, 0, std::array<double, 3>{1., 1., 1.});
  for (long x = 0; x < ovX; ++x) {
    for (long y = 0; y < ovY; ++y) {
      for (long z = 0; z < ovZ; ++z) {
        if (!is_valid_distance(0, 0, 0, x, y, z)) {
          continue;
        }
        for (long z2 = z + 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, 0, 0, x, y, z, x, y, z2);
        }
        for (long y2 = y + 1; y2 < ovZ; ++y2) {
          for (long z2 = 0; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, 0, x, y, z, x, y2, z2);
          }
        }
        for (long x2 = x + 1; x2 < ovX; ++x2) {
          for (long y2 = 0; y2 < ovY; ++y2) {
            for (long z2 = 0; z2 < ovZ; ++z2) {
              emplaceOffset(0, 0, 0, x, y, z, x2, y2, z2);
            }
          }
        }
      }
    }
  }
  // 3 planes
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z = 1; z < ovZ; ++z) {
          if (!is_valid_distance(x, y, 0, x2, 0, z)) {
            continue;
          }
          for (long y2 = 1; y2 < ovY; ++y2) {
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(x, y, 0, x2, 0, z, 0, y2, z2);
            }
          }
        }
      }
    }
  }
  // 3 edges
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long z = 1; z < ovZ; ++z) {
        emplaceOffset(x, 0, 0, 0, y, 0, 0, 0, z);
      }
    }
  }
  // 2 edges, 1 Wildcard
  // edge x, edge y
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      if (!is_valid_distance(x, 0, 0, 0, y, 0)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[x][0][0] < cells[0][y][0]) {
        emplaceOffset(x, 0, 0, x, 0, 0, 0, y, 0);
      } else {
        emplaceOffset(0, y, 0, 0, y, 0, x, 0, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long x2 = x + 1; x2 < ovX; ++x2) {
        emplaceOffset(x, 0, 0, 0, y, 0, x2, 0, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z = 1; z < ovZ; ++z) {
          emplaceOffset(x, 0, 0, 0, y, 0, x2, 0, z);
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          emplaceOffset(x, 0, 0, 0, y, 0, x2, y2, 0);
          // 2 edges, 1 inner
          for(long z = 1; z < ovZ; ++z) {
            emplaceOffset(x,0, 0, 0, y, 0, x2, y2, z);
          }
        }
      }
      for (long y2 = 1; y2 < ovY; ++y2) {
        for(long z = 1; z < ovZ; ++z) {
          emplaceOffset(x, 0, 0, 0, y, 0, 0, y2, z);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long y = 1; y < ovY; ++y) {
    for (long x = 1; x < ovX; ++x) {
      for (long y2 = y + 1; y2 < ovY; ++y2) {
        emplaceOffset(0, y, 0, x, 0, 0, 0, y2, 0);
      }
    }
  }

  // edge x, edge z
  for (long x = 1; x < ovX; ++x) {
    for (long z = 1; z < ovZ; ++z) {
      if (!is_valid_distance(x, 0, 0, 0, 0, z)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[x][0][0] < cells[0][0][z]) {
        emplaceOffset(x, 0, 0, x, 0, 0, 0, 0, z);
      } else {
        emplaceOffset(0, 0, z, 0, 0, z, x, 0, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long x2 = x + 1; x2 < ovX; ++x2) {
        emplaceOffset(x, 0, 0, 0, 0, z, x2, 0, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x2 = 1; x2 < ovX; ++x2) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(x, 0, 0, 0, 0, z, x2, 0, z2);
        }
        for (long y = 1; y < ovY; ++y) {
          emplaceOffset(x, 0, 0, 0, 0, z, x2, y, 0);
          // 2 edges, 1 inner
          for(long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(x,0, 0, 0, 0, z, x2, y, z2);
          }
        }
      }
      for (long y = 1; y < ovY; ++y) {
        for(long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(x, 0, 0, 0, 0, z, 0, y, z2);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long z = 1; z < ovY; ++z) {
    for (long x = 1; x < ovX; ++x) {
      for (long z2 = z + 1; z2 < ovY; ++z2) {
        emplaceOffset(0, 0, z, x, 0, 0, 0, 0, z2);
      }
    }
  }

  // edge y, edge z
  for (long y = 1; y < ovY; ++y) {
    for (long z = 1; z < ovZ; ++z) {
      if (!is_valid_distance(0, y, 0, 0, 0, z)) {
        continue;
      }
      // 2 edges (same cell on one edge counted twice)
      if (cells[0][y][0] < cells[0][0][z]) {
        emplaceOffset(0, y, 0, 0, y, 0, 0, 0, z);
      } else {
        emplaceOffset(0, 0, z, 0, 0, z, 0, y, 0);
      }
      // 2 edges (2 different cells on one edge)
      for (long y2 = y + 1; y2 < ovX; ++y2) {
        emplaceOffset(0, y, 0, 0, 0, z, 0, y2, 0);
      }
      // 2 edges, 1 plane/inner
      for (long x = 1; x < ovX; ++x) {
        for (long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, y, 0, 0, 0, z, x, 0, z2);
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          emplaceOffset(0, y, 0, 0, 0, z, x, y2, 0);
          // 2 edges, 1 inner
          for(long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, y, 0, 0, 0, z, x, y2, z2);
          }
        }
      }
      for (long y2 = 1; y2 < ovY; ++y2) {
        for(long z2 = 1; z2 < ovZ; ++z2) {
          emplaceOffset(0, y, 0, 0, 0, z, 0, y2, z2);
        }
      }
    }
  }
  // 2 edges (2 different cells on one edge)
  for (long z = 1; z < ovY; ++z) {
    for (long y = 1; y < ovY; ++y) {
      for (long z2 = z + 1; z2 < ovY; ++z2) {
        emplaceOffset(0, 0, z, 0, y, 0, 0, 0, z2);
      }
    }
  }

  // edge x, plane yz
  for (long x = 1; x < ovX; ++x) {
    for (long y = 1; y < ovY; ++y) {
      for (long z = 1; z < ovZ; ++z) {
        if (!is_valid_distance(x, 0, 0, 0, y, z)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[x][0][0] < cells[0][y][z]) {
          emplaceOffset(x, 0, 0, x, 0, 0, 0, y, z);
        } else {
          emplaceOffset(0, y, z, 0, y, z, x, 0, 0);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long x2 = x + 1; x2 < ovX; ++x2) {
          emplaceOffset(x, 0, 0, x2, 0, 0, 0, y, z);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            emplaceOffset(x, 0, 0, 0, y, z, x2, y2, 0);
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(x, 0, 0, 0, y, z, x2, y2, z2);
            }
          }
        }
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(x, 0, 0, 0, y, z, x2, 0, z2);
          }
        }
        // 1 edge 1 plane (2 different cells same plane)
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            if (cells[0][y][z] < cells[0][y2][z2]) {
              emplaceOffset(x, 0, 0, 0, y, z, 0, y2, z2);
            }
          }
        }
      }
    }
  }
  // edge y, plane xz
  for (long y = 1; y < ovY; ++y) {
    for (long x = 1; x < ovX; ++x) {
      for (long z = 1; z < ovZ; ++z) {
        if (!is_valid_distance(0, y, 0, x, 0, z)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[0][y][0] < cells[x][0][z]) {
          emplaceOffset(0, y, 0, 0, y, 0, x, 0, z);
        } else {
          emplaceOffset(x, 0, z, x, 0, z, 0, y, 0);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long y2 = y + 1; y2 < ovX; ++y2) {
          emplaceOffset(0, y, 0, 0, y2, 0, x, 0, z);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            emplaceOffset(0, y, 0, x, 0, z, x2, y2, 0);
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(0, y, 0, x, 0, z, x2, y2, z2);
            }
          }
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, y, 0, x, 0, z, 0, y2, z2);
          }
        }
        // 1 edge 1 plane (2 different cells same plane)
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            if (cells[x][0][z] < cells[x2][0][z2]) {
              emplaceOffset(0, y, 0, x, 0, z, x2, 0, z2);
            }
          }
        }
      }
    }
  }
  // edge z, plane xy
  for (long z = 1; z < ovZ; ++z) {
    for (long x = 1; x < ovX; ++x) {
      for (long y = 1; y < ovY; ++y) {
        if (!is_valid_distance(0, 0, z, x, y, 0)) {
          continue;
        }
        // 1 edge, 1 plane (same cell on one edge counted twice)
        if (cells[0][0][z] < cells[x][y][0]) {
          emplaceOffset(0, 0, z, 0, 0, z, x, y, 0);
        } else {
          emplaceOffset(x, y, 0, x, y, 0, 0, 0, z);
        }
        // 1 edge, 1 plane (2 different cells same edge)
        for (long z2 = z + 1; z2 < ovX; ++z2) {
          emplaceOffset(0, 0, z, 0, 0, z2, x, y, 0);
        }
        // 1 edge, 1 plane, 1 plane/inner
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long y2 = 1; y2 < ovY; ++y2) {
            // 1 edge 1 plane (2 different cells same plane)
            if (cells[x][y][0] < cells[x2][y2][0]) {
              emplaceOffset(0, 0, z, x, y, 0, x2, y2, 0);
            }
            // inner
            for (long z2 = 1; z2 < ovZ; ++z2) {
              emplaceOffset(0, 0, z, x, y, 0, x2, y2, z2);
            }
          }
        }
        for (long y2 = 1; y2 < ovY; ++y2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, z, x, y, 0, 0, y2, z2);
          }
        }
        for (long x2 = 1; x2 < ovX; ++x2) {
          for (long z2 = 1; z2 < ovZ; ++z2) {
            emplaceOffset(0, 0, z, x, y, 0, x2, 0, z2);
          }
        }
      }
    }
  }

  accumulatedDuration += std::chrono::high_resolution_clock::now() - startTime;
  std::cout << "Accumulated execution time in computeOffsets: " << accumulatedDuration.count()
              << " seconds. Size : " << _cellOffsets.size() << std::endl;

}
}  // namespace autopas
