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

 private:
  /**
   * Applies the functor to all particle pairs exploiting newtons third law of
   * motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of newton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or not:
   * - if newton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if newton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing newton3=false.
   * @tparam newton3 defines whether or not to use newton3
   * @param cell
   */
  template <bool newton3>
  void processCellAoS(ParticleCell &cell);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
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
  void processCellTripleAoSN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle triples between cell1, cell2 and cell3
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param cell3
   * @param sortingDirection Normalized vector along which particles are sorted.
   */
  void processCellTripleAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3, const std::array<double, 3> &sortingDirection);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellTripleSoAN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3);

  void processCellTripleSoANoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting.
   * @todo Currently, this is disabled because of https://github.com/AutoPas/AutoPas/issues/418
   */
  constexpr static unsigned long _startSorting = std::numeric_limits<unsigned long>::max();
};

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((DataLayout == DataLayoutOption::soa && cell._particleSoABuffer.getNumberOfParticles() == 0) ||
      (DataLayout == DataLayoutOption::aos && cell.numParticles() == 0)) {
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
  if ((DataLayout == DataLayoutOption::soa && (cell1._particleSoABuffer.getNumberOfParticles() == 0 ||
                                               cell2._particleSoABuffer.getNumberOfParticles() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.numParticles() == 0 || cell2.numParticles() == 0))) {
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
  if ((DataLayout == DataLayoutOption::soa && (cell1._particleSoABuffer.getNumberOfParticles() == 0 ||
                                               cell2._particleSoABuffer.getNumberOfParticles() == 0 ||
                                               cell3._particleSoABuffer.getNumberOfParticles() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.numParticles() == 0 || cell2.numParticles() == 0 || cell3.numParticles() == 0))) {
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
  if (cell.numParticles() > _startSorting) {
    SortedCellView<Particle, ParticleCell> cellSorted(
        cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    auto outerOuter = cellSorted._particles.begin();
    for (; outerOuter != cellSorted._particles.end(); ++outerOuter) {
      Particle &p1 = *outerOuter->second;
      auto outer = outerOuter;
      ++outer;

      for (; outer != cellSorted._particles.end(); ++outer) {
        if (std::abs(outerOuter->first - outer->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer->second;
        auto inner = outer;
        ++inner;

        for (; inner != cellSorted._particles.end(); ++inner) {
          if (std::abs(outer->first - inner->first) > _sortingCutoff || std::abs(outerOuter->first - inner->first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner->second;
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
    for (auto outerOuter = cell.begin(); outerOuter != cell.end(); ++outerOuter) {
      Particle &p1 = *outerOuter;
      auto outer = outerOuter;
      ++outer;

      for (; outer != cell.end(); ++outer) {
        Particle &p2 = *outer;
        auto inner = outer;
        ++inner;

        for (; inner != cell.end(); ++inner) {
          Particle &p3 = *inner;

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
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if (cell1.numParticles() + cell2.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    auto outerOuter = cell1Sorted._particles.begin();
    for (; outerOuter != cell1Sorted._particles.end(); ++outerOuter) {
      Particle &p1 = *outerOuter.second;

      // Particle 1 and 2 in cell 1, particle 3 in cell 2
      auto outer = outerOuter;
      ++outer;
      for (; outer != cell1Sorted._particles.end(); ++outer) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;
        for (auto &inner : cell2Sorted._particles) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff ||
              std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }

      // Particle 1 in cell 1, particle 2 and 3 in cell 2
      outer = cell2Sorted._particles.begin();
      for (; outer != cell2Sorted._particles.end(); ++outer) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;
        auto inner = outer;
        ++inner;

        for (; inner != cell2Sorted._particles.end(); ++inner) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff ||
              std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;
          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  } else {

    auto outerOuter = cell1.begin();
    auto innerStart = cell2.begin();
    for (; outerOuter != cell1.end(); ++outerOuter) {
      Particle &p1 = *outerOuter;

      // Particle 2 still in cell 1, particle 3 in cell 2
      auto outer = outerOuter;
      ++outer;
      for (; outer != cell1.end(); ++outer) {
        Particle &p2 = *outer;

        for (auto inner = innerStart; inner != cell2.end(); ++inner) {
          Particle &p3 = *inner;

          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }

      outer = innerStart;
      // Particles 2 and 3 in cell 2
      for (; outer != cell2.end(); ++outer) {
        Particle &p2 = *outer;
        auto inner = outer;
        ++inner;

        for (; inner != cell2.end(); ++inner) {
          Particle &p3 = *inner;

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
  if (cell1.numParticles() + cell2.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    for (auto outerOuter = cell1Sorted._particles.begin(); outerOuter != cell1._particles.end(); ++outerOuter) {
      Particle &p1 = *outerOuter.second;

      // Particles 1 and 2 in cell 1, particle 3 in cell 2
      auto outer = outerOuter;
      ++outer;
      for (; outer != cell1Sorted._particles.end(); ++outer) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;

        for (auto &inner : cell2Sorted._particles) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff || std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;
          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }

      // Particle 1 in cell 1, particles 2 and 3 in cell 2
      outer = cell2Sorted._particles.begin();
      for (; outer != cell2Sorted._particles.end(); ++outer) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;
        auto inner = outer;
        ++inner;

        for (; inner != cell2Sorted._particles.end(); ++inner) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff || std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;
          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  } else {
    auto innerStart = cell2.begin();

    for (auto outerOuter = cell1.begin(); outerOuter != cell1.end(); ++outerOuter) {
      Particle &p1 = *outerOuter;

      // Particles 1 and 2 in cell 1, particle 3 in cell 2
      auto outer = outerOuter;
      ++outer;
      for (; outer != cell1.end(); ++outer) {
        Particle &p2 = *outer;
        for (auto inner = innerStart; inner != cell2.end(); ++inner) {
          Particle &p3 = *inner;

          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }

      // Particle 1 in cell 1, particles 2 and 3 in cell 2
      outer = innerStart;
      for (; outer != cell2.end(); ++outer) {
        Particle &p2 = *outer;
        auto inner = outer;
        ++inner;
        for (; inner != cell2.end(); ++inner) {
          Particle &p3 = *inner;

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
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellTripleAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3, const std::array<double, 3> &sortingDirection) {
  if (cell1.numParticles() + cell2.numParticles() + cell3.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &outerOuter : cell1Sorted._particles) {
      Particle &p1 = *outerOuter.second;

      for (auto &outer : cell2Sorted._particles) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;

        for (auto &inner : cell3Sorted._particles) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff ||
              std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;

          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  } else {
    auto outerOuter = cell1.begin();
    auto outer = cell2.begin();
    auto inner = cell3.begin();

    for (; outerOuter != cell1.end(); ++outerOuter) {
      Particle &p1 = *outerOuter;

      for (; outer != cell2.end(); ++outer) {
        Particle &p2 = *outer;

        for (; inner != cell3.end(); ++inner) {
          Particle &p3 = *inner;

          _functor->AoSFunctor(p1, p2, p3, true);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3,
                                                          const std::array<double, 3> &sortingDirection) {
  if (cell1.numParticles() + cell2.numParticles() + cell3.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell3Sorted(cell3, sortingDirection);

    for (auto &outerOuter : cell1._particles) {
      Particle &p1 = *outerOuter.second;

      for (auto &outer : cell2Sorted._particles) {
        if (std::abs(outerOuter.first - outer.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *outer.second;

        for (auto &inner : cell3Sorted._particles) {
          if (std::abs(outer.first - inner.first) > _sortingCutoff || std::abs(outerOuter.first - inner.first) > _sortingCutoff) {
            break;
          }
          Particle &p3 = *inner.second;

          _functor->AoSFunctor(p1, p2, p3, false);
          if (bidirectional) {
            _functor->AoSFunctor(p2, p1, p3, false);
            _functor->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  } else {
    auto outerOuter = cell1.begin();
    auto outer = cell2.begin();
    auto inner = cell3.begin();

    for (; outerOuter != cell1.end(); ++outerOuter) {
      Particle &p1 = *outerOuter;

      for (; outer != cell2.end(); ++outer) {
        Particle &p2 = *outer;

        for (; inner != cell3.end(); ++inner) {
          Particle &p3 = *inner;

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
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairSoAN3(
    ParticleCell &cell1, ParticleCell &cell2) {
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
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellTripleSoAN3(
    ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3) {
  _functor->SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor3B<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                   bidirectional>::processCellTripleSoANoN3(ParticleCell &cell1, ParticleCell &cell2, ParticleCell &cell3) {
  _functor->SoAFunctorTriple(cell1._particleSoABuffer, cell2._particleSoABuffer, cell3._particleSoABuffer, false);
  if (bidirectional) {
    _functor->SoAFunctorTriple(cell2._particleSoABuffer, cell1._particleSoABuffer, cell3._particleSoABuffer, false);
    _functor->SoAFunctorTriple(cell3._particleSoABuffer, cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  }
}

}  // namespace autopas::internal
