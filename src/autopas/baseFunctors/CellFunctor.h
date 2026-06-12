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

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor_T. It is an internal object to handle interactions between
 * two cells of particles.
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
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell_T &cell1, ParticleCell_T &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

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

 private:
  /**
   *
   * @param particleCount
   * @param sortingDirection
   * @return
   */
  [[nodiscard]] bool shouldUseSorting(size_t particleCount, const std::array<double, 3> &sortingDirection) const {
    return particleCount >= _sortingThreshold and
           (sortingDirection[0] != 0.0 or sortingDirection[1] != 0.0 or sortingDirection[2] != 0.0);
  }

  /**
   * Applies the functor to all particle pairs exploiting Newton's third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside a cell.
   * The value of newton3 defines how to apply the aos functor:
   * - If newton3 is true: The aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - If newton3 is false: The aos functor will be applied twice for each pair (i,j and j,i), passing newton3=false.
   * @param cell
   */
  void processCellAoSImpl(ParticleCell_T &cell);

  /**
   * Applies the AoS functor to all particle pairs between cell1 and cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSImpl(ParticleCell_T &cell1, ParticleCell_T &cell2,
                              const std::array<double, 3> &sortingDirection);

  /**
   * Applies the SoA functor to all particle pairs between cell1 and cell2.
   * @param cell1
   * @param cell2
   */
  void processCellPairSoAImpl(ParticleCell_T &cell1, ParticleCell_T &cell2);

  ParticleFunctor_T &_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  const DataLayoutOption::Value _dataLayout;

  const bool _useNewton3;
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCell(ParticleCell_T &cell) {
  // avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  const bool cellHasOwnedParticles = toInt64(cell.getPossibleParticleOwnerships() & OwnershipState::owned);
  if (not cellHasOwnedParticles) {
    return;
  }

  switch (_dataLayout) {
    case DataLayoutOption::aos:
      if (cell.isEmpty()) {
        return;
      }
      processCellAoSImpl(cell);
      break;
    case DataLayoutOption::soa:
      if (cell._particleSoABuffer.size() == 0) {
        return;
      }
      _functor.SoAFunctorSingle(cell._particleSoABuffer, _useNewton3);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPair(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  const bool cell1HasOwnedParticles = toInt64(cell1.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell2HasOwnedParticles = toInt64(cell2.getPossibleParticleOwnerships() & OwnershipState::owned);

  // Avoid force calculations if both cells cannot contain owned particles.
  if (not cell1HasOwnedParticles and not cell2HasOwnedParticles) {
    return;
  }

  // In Newton3==false and functor is not bidirectional, only interactions from cell1 to cell2 are computed. Therefore,
  // if cell1 cannot contain owned particles, there is nothing useful to do.
  if constexpr (not bidirectional) {
    if (not _useNewton3 and not cell1HasOwnedParticles) {
      return;
    }
  }

  switch (_dataLayout) {
    case DataLayoutOption::aos:
      if (cell1.isEmpty() or cell2.isEmpty()) {
        return;
      }
      processCellPairAoSImpl(cell1, cell2, sortingDirection);
      break;
    case DataLayoutOption::soa:
      if (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0) {
        return;
      }
      processCellPairSoAImpl(cell1, cell2);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSImpl(ParticleCell_T &cell) {
  // helper function
  const auto interactParticles = [&f = this->_functor, n3 = this->_useNewton3](auto &p1, auto &p2) {
    if (n3) {
      f.AoSFunctor(p1, p2, true);
    } else {
      if (not p1.isHalo()) {
        f.AoSFunctor(p1, p2, false);
      }
      if (not p2.isHalo()) {
        f.AoSFunctor(p2, p1, false);
      }
    }
  };

  if (cell.size() >= _sortingThreshold) {
    SortedCellView<ParticleCell_T> cellSorted(cell, utils::ArrayMath::normalize(cell.getCellLength()));

    for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;
      // start inner loop ahead of the outer loop
      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (p2Projection - p1Projection > _sortingCutoff) {
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
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  const auto interactParticles = [&f = this->_functor, n3 = this->_useNewton3](auto &p1, auto &p2) {
    if (n3) {
      f.AoSFunctor(p1, p2, true);
    } else {
      f.AoSFunctor(p1, p2, false);

      if constexpr (bidirectional) {
        f.AoSFunctor(p2, p1, false);
      }
    }
  };

  if (shouldUseSorting(cell1.size() + cell2.size(), sortingDirection)) {
    // Use sorted cell views
    SortedCellView<ParticleCell_T> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell_T> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        // p2Projection > p1Projection guaranteed by sorting direction
        if (p2Projection - p1Projection > _sortingCutoff) {
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
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAImpl(ParticleCell_T &cell1,
                                                                                           ParticleCell_T &cell2) {
  if (_useNewton3) {
    _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
  } else {
    _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);

    if constexpr (bidirectional) {
      _functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
    }
  }
}
}  // namespace autopas::internal
