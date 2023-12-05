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
 * type ParticleFunctor. It is an internal object to handle interactions between
 * up to three cells of particles.
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which is used for particle interactions
 * @tparam dataLayout the DataLayout to be used
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used, processCellPair(..) and processCellPair(..) should handle interactions
 * for particles from all cells.
 */
template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3 = false, bool bidirectional = true>
class CellFunctor3B {
 public:
  /**
   * The constructor of CellFunctor3B.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   */
  explicit CellFunctor3B(ParticleFunctor *f, const double sortingCutoff) : _functor(f), _sortingCutoff(sortingCutoff) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All triwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell &cell);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2. This includes both 3-body
   * interactions with 2 particles from cell1 and 1 particle from cell2 as well as 1 particle from cell1 and 2 particles
   * from cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Process the interactions between 3 particles, all located in a different cell.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles will be sorted.If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellTriple(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                         const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Getter
   * @return
   */
  [[nodiscard]] DataLayoutOption::Value getDataLayout() const { return dataLayout; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getNewton3() const { return useNewton3; }

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
   * Applies the functor to all particle triplets exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of newton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or not:
   * - if newton3 is true: the aos functor will be applied once for each triplet (only i,j,k), passing newton3=true.
   * - if newton3 is false: the aos functor will be applied three times for each triplet (i,j,k and j,i,k and k,i,j),
   * passing newton3=false.
   * @tparam newton3 defines whether or not to use newton3
   * @param cell
   */
  template <bool newton3>
  void processCellAoS(ParticleCell &cell);

  /**
   * Applies the functor to all particle triplets between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle triplets between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle triples between cell1, cell2 and cell3
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles are sorted.
   */
  void processCellTripleAoSN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                              const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle triples between cell1, cell2 and cell3
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles are sorted.
   */
  void processCellTripleAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                                const std::array<double, 3> &sortingDirection);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellTripleSoAN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3);

  void processCellTripleSoANoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles of all involved cells.
   * The default threshold is taken from from CellFunctor.h. For more details see:
   * https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};
};

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::setSortingThreshold(
    size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((DataLayout == DataLayoutOption::soa && cell._particleSoABuffer.size() == 0) ||
      (DataLayout == DataLayoutOption::aos && cell.size() == 0)) {
    return;
  }

  // avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  const bool cellHasOwnedParticles = toInt64(cell.getPossibleParticleOwnerships() & OwnershipState::owned);
  if (not cellHasOwnedParticles) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      processCellAoS<useNewton3>(cell);
      break;
    case DataLayoutOption::soa:
      if constexpr (useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPair(

    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if ((DataLayout == DataLayoutOption::soa and
       (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0)) or
      (DataLayout == DataLayoutOption::aos and (cell1.size() == 0 or cell2.size() == 0))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if constexpr (useNewton3) {
        processCellPairAoSN3(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSNoN3(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if constexpr (useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellTriple(

    ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3, const std::array<double, 3> &sortingDirection) {
  if ((DataLayout == DataLayoutOption::soa and
       (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0 or
        cell3._particleSoABuffer.size() == 0)) or
      (DataLayout == DataLayoutOption::aos and (cell1.size() == 0 or cell2.size() == 0 or cell3.size() == 0))) {
    return;
  }

  // avoid force calculations if all three cells can not contain owned particles or if newton3==false and cell1 does not
  // contain owned particles
  const bool cell1HasOwnedParticles = toInt64(cell1.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell2HasOwnedParticles = toInt64(cell2.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell3HasOwnedParticles = toInt64(cell3.getPossibleParticleOwnerships() & OwnershipState::owned);

  if (((not cell1HasOwnedParticles) and (not useNewton3) and (not bidirectional)) or
      ((not cell1HasOwnedParticles) and (not cell2HasOwnedParticles) and (not cell3HasOwnedParticles))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if constexpr (useNewton3) {
        processCellTripleAoSN3(cell1, cell2, cell3, sortingDirection);
      } else {
        processCellTripleAoSNoN3(cell1, cell2, cell3, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if constexpr (useNewton3) {
        processCellTripleSoAN3(cell1, cell2, cell3);
      } else {
        processCellTripleSoANoN3(cell1, cell2, cell3);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
template <bool newton3>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellAoS(
    ParticleCell &cell) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2, auto &p3) {
    if constexpr (newton3) {
      _functor->AoSFunctor(p1, p2, p3, true);
    } else {
      if (p1.isOwned()) {
        _functor->AoSFunctor(p1, p2, p3, false);
      }
      if (p2.isOwned()) {
        _functor->AoSFunctor(p2, p1, p3, false);
      }
      if (p3.isOwned()) {
        _functor->AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if (cell.size() > _sortingThreshold) {
    SortedCellView<Particle, ParticleCell> cellSorted(
        cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;

      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }

        for (auto cellIter3 = std::next(cellIter2); cellIter3 != cellSorted._particles.end(); ++cellIter3) {
          auto &[p3Projection, p3Ptr] = *cellIter3;
          if (std::abs(p1Projection - p3Projection) > _sortingCutoff or
              std::abs(p2Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr);
        }
      }
    }
  } else {
    for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
      for (auto p2Ptr = std::next(p1Ptr); p2Ptr != cell.end(); ++p2Ptr) {
        for (auto p3Ptr = std::next(p2Ptr); p3Ptr != cell.end(); ++p3Ptr) {
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2,
                                                        const std::array<double, 3> &sortingDirection) {
  if ((cell1.size() + cell2.size() > _sortingThreshold) and sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    // Particle 1 from cell1
    for (auto cellIter1 = cell1Sorted._particles.begin(); cellIter1 != cell1Sorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;

      // Particle 2 in cell1, particle 3 in cell2
      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cell1Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        for (auto &[p3Projection, p3Ptr] : cell2Sorted._particles) {
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          _functor->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }

      // Particle 2 and 3 in cell 2
      for (auto cellIter2 = cell2Sorted._particles.begin(); cellIter2 != cell2Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        for (auto cellIter3 = std::next(cellIter2); cellIter3 != cell2Sorted._particles.end(); ++cellIter3) {
          auto &[p3Projection, p3Ptr] = *cellIter3;
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          _functor->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }
    }
  } else {  // no sorting
    // Particle 1 always from cell1
    for (auto p1Ptr = cell1.begin(); p1Ptr != cell1.end(); ++p1Ptr) {
      // Particle 2 still in cell1, particle 3 in cell2
      for (auto p2Ptr = std::next(p1Ptr); p2Ptr != cell1.end(); ++p2Ptr) {
        for (auto &p3 : cell2) {
          _functor->AoSFunctor(*p1Ptr, *p2Ptr, p3, true);
        }
      }

      // Particles 2 and 3 in cell2
      for (auto p2Ptr = cell2.begin(); p2Ptr != cell2.end(); ++p2Ptr) {
        for (auto p3Ptr = std::next(p2Ptr); p3Ptr != cell2.end(); ++p3Ptr) {
          _functor->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                          const std::array<double, 3> &sortingDirection) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2, auto &p3, const bool p2FromCell1) {
    _functor->AoSFunctor(p1, p2, p3, false);
    if (p2FromCell1) {
      _functor->AoSFunctor(p2, p1, p3, false);  // because of no newton3 and p2 is still in cell1
    } else {
      if constexpr (bidirectional) {
        if (p2.isOwned()) {
          _functor->AoSFunctor(p2, p1, p3, false);
        }
      }
    }
    if constexpr (bidirectional) {
      if (p3.isOwned()) {
        _functor->AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if ((cell1.size() + cell2.size() > _sortingThreshold) and sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    // Particle 1 always from cell1
    for (auto cellIter1 = cell1Sorted._particles.begin(); cellIter1 != cell1Sorted._particles.end(); ++cellIter1) {
      auto &[p1Projection, p1Ptr] = *cellIter1;

      // Particle 2 in cell1, particle 3 in cell2
      for (auto cellIter2 = std::next(cellIter1); cellIter2 != cell1Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        for (auto &[p3Projection, p3Ptr] : cell2Sorted._particles) {
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }

      // Particles 2 and 3 both in cell2
      for (auto cellIter2 = cell2Sorted._particles.begin(); cellIter2 != cell2Sorted._particles.end(); ++cellIter2) {
        auto &[p2Projection, p2Ptr] = *cellIter2;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        for (auto cellIter3 = std::next(cellIter2); cellIter3 != cell2Sorted._particles.end(); ++cellIter3) {
          auto &[p3Projection, p3Ptr] = *cellIter3;
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, false);
        }
      }
    }
  } else {  // no sorting
    // Particle 1 from cell1
    for (auto p1Ptr = cell1.begin(); p1Ptr != cell1.end(); ++p1Ptr) {
      // Particle 2 in cell1, particle 3 in cell2
      for (auto p2Ptr = std::next(p1Ptr); p2Ptr != cell1.end(); ++p2Ptr) {
        for (auto &p3 : cell2) {
          interactParticles(*p1Ptr, *p2Ptr, p3, true);
        }
      }

      // Particles 2 and 3 both in cell2
      for (auto p2Ptr = cell2.begin(); p2Ptr != cell2.end(); ++p2Ptr) {
        for (auto p3Ptr = std::next(p2Ptr); p3Ptr != cell2.end(); ++p3Ptr) {
          interactParticles(*p1Ptr, *p2Ptr, *p3Ptr, false);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleAoSN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                                                          const std::array<double, 3> &sortingDirection) {
  if ((cell1.size() + cell2.size() + cell3.size() > _sortingThreshold) and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection.first) > _sortingCutoff) {
          break;
        }
        for (auto &[p3Projection, p3Ptr] : cell3Sorted._particles) {
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          _functor->AoSFunctor(*p1Ptr, *p2Ptr, *p3Ptr, true);
        }
      }
    }
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        for (auto &p3 : cell3) {
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                            ParticleCell &cell3,
                                                            const std::array<double, 3> &sortingDirection) {
  // helper function
  const auto interactParticlesNoN3 = [&](auto &p1, auto &p2, auto &p3) {
    _functor->AoSFunctor(p1, p2, p3, false);
    if constexpr (bidirectional) {
      if (p2.isOwned()) {
        _functor->AoSFunctor(p2, p1, p3, false);
      }
      if (p3.isOwned()) {
        _functor->AoSFunctor(p3, p1, p2, false);
      }
    }
  };

  if (cell1.size() + cell2.size() + cell3.size() > _sortingThreshold and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }

        for (auto &[p3Projection, p3Ptr] : cell3Sorted._particles) {
          if (std::abs(p2Projection - p3Projection) > _sortingCutoff or
              std::abs(p1Projection - p3Projection) > _sortingCutoff) {
            break;
          }
          interactParticlesNoN3(*p1Ptr, *p2Ptr, *p3Ptr);
        }
      }
    }
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        for (auto &p3 : cell3) {
          interactParticlesNoN3(p1, p2, p3);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoAN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoANoN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if (bidirectional) _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleSoAN3(ParticleCell &cell1, ParticleCell &cell2,
                                                          ParticleCell &cell3) {
  _functor->SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleSoANoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                            ParticleCell &cell3) {
  _functor->SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, false);
  if (bidirectional) {
    _functor->SoAFunctorTriple(cell2._particleSoABuffer, cell1._particleSoABuffer, cell3._particleSoABuffer, false);
    _functor->SoAFunctorTriple(cell3._particleSoABuffer, cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  }
}

}  // namespace autopas::internal
