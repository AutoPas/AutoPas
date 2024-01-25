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
 * two cells of particles.
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which is used for particle interactions
 * @tparam DataLayout the DataLayout to be used
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
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
   * Process the interactions between the particles of cell1 with particles of cell2.
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
   * @param sortingDirection Normalized vector along which particles will be sorted (currently disabled)
   */
  void processCellTriple(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                         const std::array<double, 3> &sortingDirection = {0., 0., 0.});

  /**
   * Set the sorting-threshold
   * If the sum of the number of particles in three cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in three cells for which sorting should be enabled
   */
  void setSortingThreshold(size_t sortingThreshold);

 private:
  /**
   * Applies the functor to all particle triplets exploiting newtons third law of
   * motion.
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
   * Min. number of particles to start sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
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

  switch (DataLayout) {
    case DataLayoutOption::aos:
      processCellAoS<useNewton3>(cell);
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
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
  if ((DataLayout == DataLayoutOption::soa &&
       (cell1._particleSoABuffer.size() == 0 || cell2._particleSoABuffer.size() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.size() == 0 || cell2.size() == 0))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if (useNewton3) {
        processCellPairAoSN3(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSNoN3(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
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
  if ((DataLayout == DataLayoutOption::soa &&
       (cell1._particleSoABuffer.size() == 0 || cell2._particleSoABuffer.size() == 0 ||
        cell3._particleSoABuffer.size() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.size() == 0 || cell2.size() == 0 || cell3.size() == 0))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if (useNewton3) {
        processCellTripleAoSN3(cell1, cell2, cell3, sortingDirection);
      } else {
        processCellTripleAoSNoN3(cell1, cell2, cell3, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
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
  if (cell.size() > _sortingThreshold) {
    SortedCellView<Particle, ParticleCell> cellSorted(
        cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    auto p1Iter = cellSorted._particles.begin();
    for (; p1Iter != cellSorted._particles.end(); ++p1Iter) {
      Particle &p1 = *p1Iter->second;
      auto p2Iter = p1Iter;
      ++p2Iter;

      for (; p2Iter != cellSorted._particles.end(); ++p2Iter) {
        if (std::abs(p1Iter->first - p2Iter->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter->second;
        auto p3Iter = p2Iter;
        ++p3Iter;

        for (; p3Iter != cellSorted._particles.end(); ++p3Iter) {
          if (std::abs(p2Iter->first - p3Iter->first) > _sortingCutoff ||
              std::abs(p1Iter->first - p3Iter->first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter->second;
          if constexpr (newton3) {
            _functor->AoSFunctor(p1, p2, p3, true);
          } else {
            _functor->AoSFunctor(p1, p2, p3, false);
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  } else {
    for (auto p1Iter = cell.begin(); p1Iter != cell.end(); ++p1Iter) {
      Particle &p1 = *p1Iter;
      auto p2Iter = p1Iter;
      ++p2Iter;

      for (; p2Iter != cell.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;
        auto p3Iter = p2Iter;
        ++p3Iter;

        for (; p3Iter != cell.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;

          if constexpr (newton3) {
            _functor->AoSFunctor(p1, p2, p3, true);
          } else {
            _functor->AoSFunctor(p1, p2, p3, false);
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
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
    auto p1Iter = cell1Sorted._particles.begin();
    for (; p1Iter != cell1Sorted._particles.end(); ++p1Iter) {
      Particle &p1 = *p1Iter->second;

      // Particle 2 in cell1, particle 3 in cell2
      auto p2Iter = p1Iter;
      ++p2Iter;
      for (; p2Iter != cell1Sorted._particles.end(); ++p2Iter) {
        if (std::abs(p1Iter->first - p2Iter->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter->second;
        for (auto &p3Iter : cell2Sorted._particles) {
          if (std::abs(p2Iter->first - p3Iter.first) > _sortingCutoff ||
              std::abs(p1Iter->first - p3Iter.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter.second;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }

      // Particle 2 and 3 in cell 2
      p2Iter = cell2Sorted._particles.begin();
      for (; p2Iter != cell2Sorted._particles.end(); ++p2Iter) {
        if (std::abs(p1Iter->first - p2Iter->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter->second;
        auto p3Iter = p2Iter;
        ++p3Iter;

        for (; p3Iter != cell2Sorted._particles.end(); ++p3Iter) {
          if (std::abs(p2Iter->first - p3Iter->first) > _sortingCutoff ||
              std::abs(p1Iter->first - p3Iter->first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter->second;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  } else {
    // Particle 1 always from cell1
    for (auto p1Iter = cell1.begin(); p1Iter != cell1.end(); ++p1Iter) {
      Particle &p1 = *p1Iter;

      // Particle 2 still in cell1, particle 3 in cell2
      auto p2Iter = p1Iter;
      ++p2Iter;
      for (; p2Iter != cell1.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;

        for (auto p3Iter = cell2.begin(); p3Iter != cell2.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }

      // Particles 2 and 3 in cell2
      for (p2Iter = cell2.begin(); p2Iter != cell2.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;

        auto p3Iter = p2Iter;
        ++p3Iter;
        for (; p3Iter != cell2.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;
          _functor->AoSFunctor(p1, p2, p3, true);
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
  if ((cell1.size() + cell2.size() > _sortingThreshold) and sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    // Particle 1 always from cell1
    for (auto p1Iter = cell1Sorted._particles.begin(); p1Iter != cell1Sorted._particles.end(); ++p1Iter) {
      Particle &p1 = *p1Iter->second;

      // Particle 2 in cell1, particle 3 in cell2
      auto p2Iter = p1Iter;
      ++p2Iter;
      for (; p2Iter != cell1Sorted._particles.end(); ++p2Iter) {
        if (std::abs(p1Iter->first - p2Iter->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter->second;

        for (auto &p3Iter : cell2Sorted._particles) {
          if (std::abs(p2Iter->first - p3Iter.first) > _sortingCutoff ||
              std::abs(p1Iter->first - p3Iter.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter.second;
          _functor->AoSFunctor(p1, p2, p3, false);
          _functor->AoSFunctor(p2, p1, p3, false);  // because of no newton3 and p2 is still in cell1
          if (bidirectional) {
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }

      // Particles 2 and 3 both in cell2
      p2Iter = cell2Sorted._particles.begin();
      for (; p2Iter != cell2Sorted._particles.end(); ++p2Iter) {
        if (std::abs(p1Iter->first - p2Iter->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter->second;

        auto p3Iter = p2Iter;
        ++p3Iter;
        for (; p3Iter != cell2Sorted._particles.end(); ++p3Iter) {
          if (std::abs(p2Iter->first - p3Iter->first) > _sortingCutoff ||
              std::abs(p1Iter->first - p3Iter->first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter->second;
          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  } else {
    // Particle 1 from cell1
    for (auto p1Iter = cell1.begin(); p1Iter != cell1.end(); ++p1Iter) {
      Particle &p1 = *p1Iter;

      // Particle 2 in cell1, particle 3 in cell2
      auto p2Iter = p1Iter;
      ++p2Iter;
      for (; p2Iter != cell1.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;

        for (auto p3Iter = cell2.begin(); p3Iter != cell2.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;
          _functor->AoSFunctor(p1, p2, p3, false);
          _functor->AoSFunctor(p2, p1, p3, false);  // because of no newton3 and p2 is still in cell1
          if (bidirectional) {
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }

      // Particles 2 and 3 both in cell2
      for (p2Iter = cell2.begin(); p2Iter != cell2.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;

        auto p3Iter = p2Iter;
        ++p3Iter;
        for (; p3Iter != cell2.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;
          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
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
  auto cell1Count = 0;
  auto cell2Count = 0;
  auto cell3Count = 0;
  if ((cell1.size() + cell2.size() + cell3.size() > _sortingThreshold) and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &p1Iter : cell1Sorted._particles) {
      Particle &p1 = *p1Iter.second;

      for (auto &p2Iter : cell2Sorted._particles) {
        if (std::abs(p1Iter.first - p2Iter.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter.second;

        for (auto &p3Iter : cell3Sorted._particles) {
          if (std::abs(p2Iter.first - p3Iter.first) > _sortingCutoff ||
              std::abs(p1Iter.first - p3Iter.first) > _sortingCutoff) {
            break;
          }

          Particle &p3 = *p3Iter.second;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  } else {
    for (auto p1Iter = cell1.begin(); p1Iter != cell1.end(); ++p1Iter) {
      Particle &p1 = *p1Iter;
      cell1Count++;

      for (auto p2Iter = cell2.begin(); p2Iter != cell2.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;
        cell2Count++;

        for (auto p3Iter = cell3.begin(); p3Iter != cell3.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;
          cell3Count++;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  }
  std::cout << "Particles for cells 1 2 3 : " << cell1Count << "  " << cell2Count << "  " << cell3Count << std::endl;
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                            ParticleCell &cell3,
                                                            const std::array<double, 3> &sortingDirection) {
  if (cell1.size() + cell2.size() + cell3.size() > _sortingThreshold and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &p1Iter : cell1Sorted._particles) {
      Particle &p1 = *p1Iter.second;

      for (auto &p2Iter : cell2Sorted._particles) {
        if (std::abs(p1Iter.first - p2Iter.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *p2Iter.second;

        for (auto &p3Iter : cell3Sorted._particles) {
          if (std::abs(p2Iter.first - p3Iter.first) > _sortingCutoff ||
              std::abs(p1Iter.first - p3Iter.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *p3Iter.second;

          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  } else {
    for (auto p1Iter = cell1.begin(); p1Iter != cell1.end(); ++p1Iter) {
      Particle &p1 = *p1Iter;

      for (auto p2Iter = cell2.begin(); p2Iter != cell2.end(); ++p2Iter) {
        Particle &p2 = *p2Iter;

        for (auto p3Iter = cell3.begin(); p3Iter != cell3.end(); ++p3Iter) {
          Particle &p3 = *p3Iter;

          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
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
