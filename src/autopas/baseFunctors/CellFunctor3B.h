/**
 * @file CellFunctor3B.h
 * @author M. Muehlhaeusser
 * @date 25/07/23
 */

#pragma once

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor_T. It is an internal object to handle interactions between
 * up to three cells of particles.
 * @tparam ParticleCell_T
 * @tparam ParticleFunctor_T the functor which is used for particle interactions
 * @tparam bidirectional if no newton3 is used, processCellPair(..) and processCellPair(..) should handle interactions
 * for particles from all cells.
 */
template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional = true>
class CellFunctor3B {
 public:
  /**
   * The constructor of CellFunctor3B.
   * @param f The particle functor, which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit CellFunctor3B(ParticleFunctor_T &f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All triwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell_T &cell);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2. This includes both triwise
   * interactions with 2 particles from cell1 and 1 particle from cell2 as well as 1 particle from cell1 and 2 particles
   * from cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell_T &cell1, ParticleCell_T &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Process the interactions between 3 particles, all located in a different cell.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles will be sorted.If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellTriple(ParticleCell_T &cell1, ParticleCell_T &cell2, ParticleCell_T &cell3,
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
   * If the sum of the number of particles in three cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in three cells from which sorting should be enabled.
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
   * Applies the functor to all particle triplets exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside a cell.
   * The value of newton3 defines how to apply the aos functor:
   * - If _useNewton3 is true: The aos functor will be applied once for each triplet (only i,j,k), passing newton3=true.
   * - If _useNewton3 is false: The aos functor will be applied three times for each triplet (i,j,k and j,i,k and
   * k,i,j), passing newton3=false.
   * @tparam newton3
   * @param cell
   */
  template <bool newton3>
  void processCellAoSImpl(ParticleCell_T &cell);

  /**
   * Applies the functor to all particle triplets between cell1 and cell2
   * exploiting newtons third law of motion.
   * @tparam newton3
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  template <bool newton3>
  void processCellPairAoSImpl(ParticleCell_T &cell1, ParticleCell_T &cell2,
                              const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle triples between cell1, cell2 and cell3.
   * @tparam newton3 Whether to use Newton's third law of motion.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles from cell1 and cell2 are sorted.
   *
   * @warning If sorting is used, sortingDirection should correspond to the pairwise sorting direction between cell1 and
   * cell2.
   */
  template <bool newton3>
  void processCellTripleAoSImpl(ParticleCell_T &cell1, ParticleCell_T &cell2, ParticleCell_T &cell3,
                                const std::array<double, 3> &sortingDirection);

  /**
   * Applies the SoA functor to all particle triplets between cell1 and cell2.
   * @tparam newton3 Whether to use Newton's third law of motion.
   * @param cell1
   * @param cell2
   */
  template <bool newton3>
  void processCellPairSoAImpl(ParticleCell_T &cell1, ParticleCell_T &cell2);

  /**
   * Applies the SoA functor to all particle triplets between cell1, cell2 and cell3.
   * @tparam newton3 Whether to use Newton's third law of motion.
   * @param cell1
   * @param cell2
   * @param cell3
   */
  template <bool newton3>
  void processCellTripleSoAImpl(ParticleCell_T &cell1, ParticleCell_T &cell2, ParticleCell_T &cell3);

  ParticleFunctor_T &_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles of all involved cells.
   * The default threshold is (blindly) taken from CellFunctor.h. For some more details see:
   * https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  const DataLayoutOption::Value _dataLayout;

  const bool _useNewton3;
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCell(ParticleCell_T &cell) {
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
      if (_useNewton3) {
        processCellAoSImpl<true>(cell);
      } else {
        processCellAoSImpl<false>(cell);
      }
      break;
    case DataLayoutOption::soa:
      if (cell._particleSoABuffer.size() == 0) {
        return;
      }
      if (_useNewton3) {
        _functor.SoAFunctorSingle(cell._particleSoABuffer, true);
      } else {
        _functor.SoAFunctorSingle(cell._particleSoABuffer, false);
      }
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPair(

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
      if (_useNewton3) {
        processCellPairAoSImpl<true>(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSImpl<false>(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0) {
        return;
      }
      if (_useNewton3) {
        processCellPairSoAImpl<true>(cell1, cell2);
      } else {
        processCellPairSoAImpl<false>(cell1, cell2);
      }
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellTriple(

    ParticleCell_T &cell1, ParticleCell_T &cell2, ParticleCell_T &cell3,
    const std::array<double, 3> &sortingDirection) {
  const bool cell1HasOwnedParticles = toInt64(cell1.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell2HasOwnedParticles = toInt64(cell2.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell3HasOwnedParticles = toInt64(cell3.getPossibleParticleOwnerships() & OwnershipState::owned);

  // Avoid force calculations if all three cells cannot contain owned particles.
  if (not cell1HasOwnedParticles and not cell2HasOwnedParticles and not cell3HasOwnedParticles) {
    return;
  }

  // In Newton3==false and functor is not bidirectional, only interactions from cell1 are computed. Therefore,
  // if cell1 cannot contain owned particles, there is nothing useful to do.
  if constexpr (not bidirectional) {
    if (not _useNewton3 and not cell1HasOwnedParticles) {
      return;
    }
  }

  switch (_dataLayout) {
    case DataLayoutOption::aos:
      if (cell1.isEmpty() or cell2.isEmpty() or cell3.isEmpty()) {
        return;
      }
      if (_useNewton3) {
        processCellTripleAoSImpl<true>(cell1, cell2, cell3, sortingDirection);
      } else {
        processCellTripleAoSImpl<false>(cell1, cell2, cell3, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0 or
          cell3._particleSoABuffer.size() == 0) {
        return;
      }
      if (_useNewton3) {
        processCellTripleSoAImpl<true>(cell1, cell2, cell3);
      } else {
        processCellTripleSoAImpl<false>(cell1, cell2, cell3);
      }
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSImpl(ParticleCell_T &cell) {
  // helper function
  const auto interactParticles = [&f = this->_functor](auto &p1, auto &p2, auto &p3) {
    if constexpr (newton3) {
      f.AoSFunctor(p1, p2, p3, true);
    } else {
      if (not p1.isHalo()) {
        f.AoSFunctor(p1, p2, p3, false);
      }
      if (not p2.isHalo()) {
        f.AoSFunctor(p2, p1, p3, false);
      }
      if (not p3.isHalo()) {
        f.AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if (cell.size() > _sortingThreshold) {
    SortedCellView<ParticleCell_T> cellSorted(cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;

      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (p2Projection - p1Projection > _sortingCutoff) {
          break;
        }

        for (auto cellIter3 = std::next(cellIter2); cellIter3 != cellSorted._particles.end(); ++cellIter3) {
          auto &[p3Projection, p3Ptr] = *cellIter3;
          if (p3Projection - p1Projection > _sortingCutoff or p3Projection - p2Projection > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr);
        }
      }
    }
  } else {
    for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
      auto p2Ptr = p1Ptr;
      ++p2Ptr;
      for (; p2Ptr != cell.end(); ++p2Ptr) {
        auto p3Ptr = p2Ptr;
        ++p3Ptr;
        for (; p3Ptr != cell.end(); ++p3Ptr) {
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr);
        }
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  const auto interactParticles = [&f = this->_functor](auto &p1, auto &p2, auto &p3, const bool p2FromCell1) {
    if constexpr (newton3) {
      f.AoSFunctor(p1, p2, p3, true);
    } else {
      f.AoSFunctor(p1, p2, p3, false);
      if (p2FromCell1) {
        f.AoSFunctor(p2, p1, p3, false);  // because of no newton3 and p2 is still in cell1
      } else {
        if constexpr (bidirectional) {
          f.AoSFunctor(p2, p1, p3, false);
        }
      }
      if constexpr (bidirectional) {
        f.AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if (shouldUseSorting(cell1.size() + cell2.size(), sortingDirection)) {
    SortedCellView<ParticleCell_T> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell_T> cell2Sorted(cell2, sortingDirection);

    // Particle 1 from cell1
    for (auto cellIter1 = cell1Sorted._particles.begin(); cellIter1 != cell1Sorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;

      // Particle 2 in cell1, particle 3 in cell2
      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cell1Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (p2Projection - p1Projection > _sortingCutoff) {
          break;
        }
        for (auto &[p3Projection, p3Ptr] : cell2Sorted._particles) {
          if (p3Projection - p2Projection > _sortingCutoff or p3Projection - p1Projection > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }

      // Particle 2 and 3 in cell 2
      for (auto cellIter2 = cell2Sorted._particles.begin(); cellIter2 != cell2Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (p2Projection - p1Projection > _sortingCutoff) {
          break;
        }
        for (auto cellIter3 = std::next(cellIter2); cellIter3 != cell2Sorted._particles.end(); ++cellIter3) {
          auto &[p3Projection, p3Ptr] = *cellIter3;
          if (p3Projection - p2Projection > _sortingCutoff or p3Projection - p1Projection > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, false);
        }
      }
    }
  } else {  // no sorting
    // Particle 1 always from cell1
    for (auto p1Ptr = cell1.begin(); p1Ptr != cell1.end(); ++p1Ptr) {
      // Particle 2 still in cell1, particle 3 in cell2
      auto p2Ptr = p1Ptr;
      ++p2Ptr;
      for (; p2Ptr != cell1.end(); ++p2Ptr) {
        for (auto &p3 : cell2) {
          interactParticles(*p1Ptr, *p2Ptr, p3, true);
        }
      }

      // Particles 2 and 3 in cell2
      for (auto p2Ptr = cell2.begin(); p2Ptr != cell2.end(); ++p2Ptr) {
        auto p3Ptr = p2Ptr;
        ++p3Ptr;
        for (; p3Ptr != cell2.end(); ++p3Ptr) {
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, false);
        }
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellTripleAoSImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, ParticleCell_T &cell3,
    const std::array<double, 3> &sortingDirection) {
  const auto interactParticles = [&f = this->_functor](auto &p1, auto &p2, auto &p3) {
    if constexpr (newton3) {
      f.AoSFunctor(p1, p2, p3, true);
    } else {
      f.AoSFunctor(p1, p2, p3, false);

      if constexpr (bidirectional) {
        f.AoSFunctor(p2, p1, p3, false);
        f.AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if (shouldUseSorting(cell1.size() + cell2.size() + cell3.size(), sortingDirection)) {
    SortedCellView<ParticleCell_T> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell_T> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        // p2Projection > p1Projection guaranteed by sorting direction
        if (p2Projection - p1Projection > _sortingCutoff) {
          break;
        }
        for (auto &p3 : cell3) {
          interactParticles(*p1Ptr, *p2Ptr, p3);
        }
      }
    }
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        for (auto &p3 : cell3) {
          interactParticles(p1, p2, p3);
        }
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAImpl(ParticleCell_T &cell1,
                                                                                             ParticleCell_T &cell2) {
  if constexpr (newton3) {
    _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
  } else {
    _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
    if constexpr (bidirectional) {
      _functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor3B<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellTripleSoAImpl(ParticleCell_T &cell1,
                                                                                               ParticleCell_T &cell2,
                                                                                               ParticleCell_T &cell3) {
  if constexpr (newton3) {
    _functor.SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, true);
  } else {
    _functor.SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, false);
    if constexpr (bidirectional) {
      _functor.SoAFunctorTriple(cell2._particleSoABuffer, cell1._particleSoABuffer, cell3._particleSoABuffer, false);
      _functor.SoAFunctorTriple(cell3._particleSoABuffer, cell1._particleSoABuffer, cell2._particleSoABuffer, false);
    }
  }
}
}  // namespace autopas::internal
