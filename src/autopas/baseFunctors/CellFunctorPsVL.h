/**
 * @file CellFunctorPsVL.h
 *
 * @date 06.12.2025
 * @author Lars Doll
 */

#pragma once

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles using the PseudoVerletList algorithm.
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell, class ParticleFunctor, bool bidirectional = true>
class CellFunctorPsVL {
 public:
  /**
   * The constructor of CellFunctorPsVL.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This parameter normally should be the cutoff + skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit CellFunctorPsVL(ParticleFunctor *f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cellIndex
   */
  void processCell(unsigned long cellIndex);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1Index
   * @param cell2Index
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(unsigned long cell1Index, unsigned long cell2Index,
    const std::array<double, 3> &sortingDirection);

  /**
   * Getter
   * @return
   */
  [[nodiscard]] DataLayoutOption getDataLayout() const { return _dataLayout; }

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

  void setOrientationLists(std::vector<std::vector<SortedCellView<ParticleCell>>> &lists);

 private:

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of _useNewton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or
   * not:
   * - if _useNewton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if _useNewton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
   * @param cellIndex
   */
  void processCellAoS(unsigned long cellIndex);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1Index
   * @param cell2Index
   * @param directionIndex Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1Index
   * @param cell2Index
   * @param directionIndex Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  signed long getDirectionIndex(const std::array<double, 3> &sortingDirection);

  signed long flipDirectionIndex(signed long directionIndex);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  DataLayoutOption _dataLayout;

  bool _useNewton3;

  std::vector<std::vector<SortedCellView<ParticleCell>>>* _orientationList = nullptr;
};

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
signed long CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::flipDirectionIndex(signed long directionIndex) {
  return -(directionIndex+2);
}


template <class ParticleCell, class ParticleFunctor, bool bidirectional>
signed long CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::getDirectionIndex(
    const std::array<double, 3> &sortingDirection) {
  const auto& [x, y, z] = sortingDirection;
  int xInt = (x > 0) - (x < 0);
  int yInt = (y > 0) - (y < 0);
  int zInt = (z > 0) - (z < 0);

  constexpr std::array<int, 3> dims{3, 3, 3};
  return autopas::utils::ThreeDimensionalMapping::threeToOneD(xInt+1,yInt+1,zInt+1,dims) - 14;
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::setOrientationLists(
    std::vector<std::vector<SortedCellView<ParticleCell>>> &lists) {
  _orientationList = &lists;
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCell(unsigned long cellIndex) {

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      processCellAoS(cellIndex);
      break;
    case DataLayoutOption::soa:

      break;
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPair(
  unsigned long cell1Index, unsigned long cell2Index, const std::array<double, 3> &sortingDirection) {

  signed long directionIndex = getDirectionIndex(sortingDirection);

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      if (_useNewton3) {
        processCellPairAoSN3(cell1Index, cell2Index, directionIndex);
      } else {
        processCellPairAoSNoN3(cell1Index, cell2Index, directionIndex);
      }
      break;
    case DataLayoutOption::soa:

      break;
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellAoS(unsigned long cellIndex) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2) {
    if (_useNewton3) {
      _functor->AoSFunctor(p1, p2, true);
    } else {
      if (not p1.isHalo()) {
        _functor->AoSFunctor(p1, p2, false);
      }
      if (not p2.isHalo()) {
        _functor->AoSFunctor(p2, p1, false);
      }
    }
  };

  auto &cellSorted = (*_orientationList)[cellIndex][6]; //Index 6 for the cell diagonal (direction {1,1,1})

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

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSN3(
  unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex) {
  bool revSecondView = false;
  if (directionIndex<0) {
    revSecondView = true;
    directionIndex = flipDirectionIndex(directionIndex);
  }
  auto &cell1Sorted = (*_orientationList)[cell1Index][directionIndex];
  auto &cell2Sorted = (*_orientationList)[cell2Index][directionIndex];

  for (auto &[p1Proj, p1Ptr] : cell1Sorted._particles) {
    if (revSecondView == false) {
      //forward
      for (auto it = cell2Sorted._particles.begin(); it != cell2Sorted._particles.end(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        _functor->AoSFunctor(*p1Ptr, *p2Ptr, true);
           }
    } else {
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        _functor->AoSFunctor(*p1Ptr, *p2Ptr, true);
           }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSNoN3(unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex) {
  const auto interactParticlesNoN3 = [&](auto &p1, auto &p2) {
    _functor->AoSFunctor(p1, p2, false);
    if constexpr (bidirectional) {
      _functor->AoSFunctor(p2, p1, false);
    }
  };

  bool revSecondView = false;
  if (directionIndex<0) {
    revSecondView = true;
    directionIndex = flipDirectionIndex(directionIndex);
  }

  auto &cell1Sorted = (*_orientationList)[cell1Index][directionIndex];
  auto &cell2Sorted = (*_orientationList)[cell2Index][directionIndex];

  for (auto &[p1Proj, p1Ptr] : cell1Sorted._particles) {
    if (revSecondView == false) {
      //forward
      for (auto it = cell2Sorted._particles.begin(); it != cell2Sorted._particles.end(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1Ptr, *p2Ptr);
      }
    } else {
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoAN3(ParticleCell &cell) {

}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoANoN3(ParticleCell &cell) {

}
}  // namespace autopas::internal
