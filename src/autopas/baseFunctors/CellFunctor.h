/**
 * @file CellFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor_T. It is an internal object to handle interactions between
 * two cells of particles. Also supports PsVL.
 * @tparam ParticleCell_T
 * @tparam ParticleFunctor_T the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The particle functor, which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit CellFunctor(ParticleFunctor_T &f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell_T &cell);

  /**
   * Process the interactions inside one cell for PseudoVerletLists.
   * @param cellIndex All pairwise interactions of particles inside this cell are calculated.
   */
  void processCellPsVL(unsigned long cellIndex);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell_T &cell1, ParticleCell_T &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Process the interactions between the particles of cell1 with particles of cell2 for PseudoVerletLists.
   * @param cell1Index
   * @param cell2Index
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairPsVL(unsigned long cell1Index, unsigned long cell2Index,
                           const std::array<double, 3> &sortingDirection);

  /**
   * Getter
   * @return
   */
  [[nodiscard]] DataLayoutOption::Value getDataLayout() const { return _dataLayout; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getNewton3() const { return _useNewton3; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getBidirectional() const { return bidirectional; }

  /**
   * Set the sorting-threshold
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in two cells from which sorting should be enabled.
   */
  void setSortingThreshold(size_t sortingThreshold);

  /**
   * Sets the orientationList used by PseudoVerletLists.
   * @param list the orientationList
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list);

 private:
  /**
   * Evaluate whether the AoSFunctor should use sorting, depending on the set sorting threshold.
   * @param particleCount Total number of involved particles.
   * @param sortingDirection No sorting when the sorting direction is {0., 0., 0.}.
   * @return whether the AoSFunctor should use the SortedCellView.
   */
  [[nodiscard]] bool shouldUseSorting(size_t particleCount, const std::array<double, 3> &sortingDirection) const {
    return particleCount >= _sortingThreshold and
           (sortingDirection[0] != 0.0 or sortingDirection[1] != 0.0 or sortingDirection[2] != 0.0);
  }

  /**
   * Applies the functor to all particle pairs exploiting Newton's third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside a cell.
   * The value of _useNewton3 defines how to apply the aos functor:
   * - If _useNewton3 is true: The aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - If _useNewton3 is false: The aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
   * @param cell
   */
  void processCellAoSImpl(ParticleCell_T &cell);

  /**
   * Applies the functor to all particle pairs inside one cell for PseudoVerletLists, using the precomputed
   * SortedCellView along the cell diagonal.
   * @param cellIndex
   */
  void processCellAoSPsVLImpl(unsigned long cellIndex);

  /**
   * Applies the AoS functor to all particle pairs between cell1 and cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSImpl(ParticleCell_T &cell1, ParticleCell_T &cell2,
                              const std::array<double, 3> &sortingDirection);

  /**
   * Applies the AoS functor to all particle pairs between cell1 and cell2 for PseudoVerletLists, using the precomputed
   * SortedCellViews for the given direction.
   * @param cell1Index
   * @param cell2Index
   * @param directionIndex Index of the direction connecting the centers of cell1 and cell2.
   */
  void processCellPairAoSPsVLImpl(unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex);

  /**
   * Applies the SoA functor to all particle pairs between cell1 and cell2.
   * @param cell1
   * @param cell2
   */
  void processCellPairSoAImpl(ParticleCell_T &cell1, ParticleCell_T &cell2);

  /**
   * Utilizes the threeToOneD mapping to calculate the directionIndex relative to the center of a 3x3 cube of cells.
   * For PseudoVerletLists.
   * @param sortingDirection
   * @return directionIndex
   */
  static signed long getDirectionIndex(const std::array<double, 3> &sortingDirection);

  /**
   * Flips negative directionIndices for PsVL as only positive directions are stored in _orientationList.
   * For PseudoVerletLists.
   * @param directionIndex
   * @return flipped directionIndex
   */
  static signed long flipDirectionIndex(signed long directionIndex);

  /**
   * The functor used for particle interactions.
   */
  ParticleFunctor_T &_functor;

  /**
   * This parameter indicates the maximal distance the sorted particles are to interact.
   */
  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  /**
   * This parameter determines which data layout is used.
   */
  const DataLayoutOption::Value _dataLayout;

  /**
   * This parameter determines if the newton3 optimization is used.
   */
  const bool _useNewton3;

  /**
   * Orientation List: For each cell, 13 sortedCellViews are stored, each of which sorts in the direction of the
   * neighboring cell. Used by PseudoVerletLists.
   */
  std::vector<std::vector<SortedCellView<ParticleCell_T>>> *_orientationList = nullptr;
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) {
  _orientationList = &list;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
signed long CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::flipDirectionIndex(
    signed long directionIndex) {
  return -(directionIndex + 2);
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
signed long CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::getDirectionIndex(
    const std::array<double, 3> &sortingDirection) {
  const auto &[x, y, z] = sortingDirection;
  // each direction get values in {-1,0,1}, according to the sign of the direction in sortingDirection.
  int xDirection = (x > 0) - (x < 0);
  int yDirection = (y > 0) - (y < 0);
  int zDirection = (z > 0) - (z < 0);

  constexpr std::array<int, 3> dims{3, 3, 3};
  // subtract an offset of 14 so the direction with index 0 is the first vector in positive direction.
  return autopas::utils::ThreeDimensionalMapping::threeToOneD(xDirection + 1, yDirection + 1, zDirection + 1, dims) -
         14;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCell(ParticleCell_T &cell) {
  const bool isAoS = _dataLayout == DataLayoutOption::aos ? true : false;
  const bool isSoA = _dataLayout == DataLayoutOption::soa ? true : false;

  // Return early if the cell is empty.
  if ((isSoA and cell._particleSoABuffer.size() == 0) or (isAoS and cell.isEmpty())) {
    return;
  }
  // Avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  if (not cell.canHaveOwnedParticles()) {
    return;
  }

  if (isAoS) {
    processCellAoSImpl(cell);
  } else if (isSoA) {
    _functor.SoAFunctorSingle(cell._particleSoABuffer, _useNewton3);
  } else {
    utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                       DataLayoutOption(_dataLayout));
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPsVL(unsigned long cellIndex) {
  const bool isAoS = _dataLayout == DataLayoutOption::aos ? true : false;
  const bool isSoA = _dataLayout == DataLayoutOption::soa ? true : false;

  // Empty-cell and ownership guards are handled by the PsVL traversal - no checking here unlike in processCell

  if (isAoS) {
    processCellAoSPsVLImpl(cellIndex);
  } else if (isSoA) {
    utils::ExceptionHandler::exception("PsVL currently don't support SoA!");
  } else {
    utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                       DataLayoutOption(_dataLayout));
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPair(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  const bool isAoS = _dataLayout == DataLayoutOption::aos ? true : false;
  const bool isSoA = _dataLayout == DataLayoutOption::soa ? true : false;

  // Return early if a cell is empty.
  if ((isSoA and (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0)) or
      (isAoS and (cell1.isEmpty() or cell2.isEmpty()))) {
    return;
  }

  if (not cell1.canHaveOwnedParticles()) {
    // Nothing to do if cell1 has no owned particles and we don't write to cell2 particles.
    if constexpr (not bidirectional) {
      if (not _useNewton3) {
        return;
      }
    }
    // Nothing to do if both cells cannot have owned particles.
    if (not cell2.canHaveOwnedParticles()) {
      return;
    }
  }

  if (isAoS) {
    processCellPairAoSImpl(cell1, cell2, sortingDirection);
  } else if (isSoA) {
    processCellPairSoAImpl(cell1, cell2);
  } else {
    utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                       DataLayoutOption(_dataLayout));
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairPsVL(
    unsigned long cell1Index, unsigned long cell2Index, const std::array<double, 3> &sortingDirection) {
  const bool isAoS = _dataLayout == DataLayoutOption::aos ? true : false;
  const bool isSoA = _dataLayout == DataLayoutOption::soa ? true : false;

  // Empty-cell and ownership guards are handled by the PsVL traversal - no checking here unlike in processCellPair

  if (isAoS) {
    const signed long directionIndex = getDirectionIndex(sortingDirection);
    processCellPairAoSPsVLImpl(cell1Index, cell2Index, directionIndex);
  } else if (isSoA) {
    utils::ExceptionHandler::exception("PsVL currently don't support SoA!");
  } else {
    utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                       DataLayoutOption(_dataLayout));
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSImpl(ParticleCell_T &cell) {
  // helper function
  const auto interactParticles = [this](auto &p1, auto &p2) {
    this->_functor.AoSFunctor(p1, p2, this->_useNewton3);
    if (not this->_useNewton3) {
      this->_functor.AoSFunctor(p2, p1, false);
    }
  };

  if (cell.size() >= _sortingThreshold) {
    SortedCellView<ParticleCell_T> cellSorted(cell, utils::ArrayMath::normalize(cell.getCellLength()));

    for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;
      // start inner loop ahead of the outer loop
      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p2Projection - p1Projection) > _sortingCutoff) {
          break;
        }
        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  } else {
    for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
      auto p2Ptr = p1Ptr;
      ++p2Ptr;
      for (; p2Ptr != cell.end(); ++p2Ptr) {
        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSPsVLImpl(unsigned long cellIndex) {
  // helper function
  const auto interactParticles = [this](auto &p1, auto &p2) {
    this->_functor.AoSFunctor(p1, p2, this->_useNewton3);
    if (not this->_useNewton3) {
      this->_functor.AoSFunctor(p2, p1, false);
    }
  };

  auto &cellSorted = (*_orientationList)[cellIndex][6];  // Index 6 for the cell diagonal (direction {1,1,1})
  for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
    auto &[p1Projection, p1Ptr] = *cellIter1;
    // start inner loop ahead of the outer loop
    for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
      auto &[p2Projection, p2Ptr] = *cellIter2;
      if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
        break;
      }
      interactParticles(*p1Ptr, *p2Ptr);
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  const auto interactParticles = [this](auto &p1, auto &p2) {
    this->_functor.AoSFunctor(p1, p2, this->_useNewton3);
    if constexpr (bidirectional) {
      if (not this->_useNewton3) {
        this->_functor.AoSFunctor(p2, p1, false);
      }
    }
  };

  if (shouldUseSorting(cell1.size() + cell2.size(), sortingDirection)) {
    // Use sorted cell views
    SortedCellView<ParticleCell_T> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell_T> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p2Projection - p1Projection) > _sortingCutoff) {
          break;
        }

        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  } else {
    // Without sorting
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        interactParticles(p1, p2);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSPsVLImpl(
    unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex) {
  const auto interactParticles = [this](auto &p1, auto &p2) {
    this->_functor.AoSFunctor(p1, p2, this->_useNewton3);
    if constexpr (bidirectional) {
      if (not this->_useNewton3) {
        this->_functor.AoSFunctor(p2, p1, false);
      }
    }
  };

  bool reverseCell2View = false;
  if (directionIndex < 0) {
    reverseCell2View = true;
    directionIndex = flipDirectionIndex(directionIndex);
  }
  auto &cell1Sorted = (*_orientationList)[cell1Index][directionIndex];
  auto &cell2Sorted = (*_orientationList)[cell2Index][directionIndex];

  if (not reverseCell2View) {
    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  } else {
    for (auto &[p1Proj, p1Ptr] : cell1Sorted._particles) {
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAImpl(ParticleCell_T &cell1,
                                                                                           ParticleCell_T &cell2) {
  _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, _useNewton3);
  if constexpr (bidirectional) {
    if (not _useNewton3) {
      _functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
    }
  }
}
}  // namespace autopas::internal