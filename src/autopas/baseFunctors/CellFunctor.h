/**
 * @file CellFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <algorithm>

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoASortedView.h"

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

  /**
   * Set the SoA sorting-threshold.
   * If the sum of the number of particles in two SoA buffers is greater than this value, the SoA path uses
   * SoAFunctorPairSorted instead of SoAFunctorPair.
   * @param soaSortingThreshold Sum of the number of particles in two cells from which SoA sorting should be enabled.
   */
  void setSoASortingThreshold(size_t soaSortingThreshold);

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

  [[nodiscard]] bool shouldUseSoASorting(size_t particleCount, const std::array<double, 3> &sortingDirection) const {
    return particleCount >= _soaSortingThreshold and
           (sortingDirection[0] != 0.0 or sortingDirection[1] != 0.0 or sortingDirection[2] != 0.0);
  }

  /**
   * Applies the functor to all particle pairs exploiting Newton's third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside a cell.
   * The value of newton3 defines how to apply the aos functor:
   * - If _useNewton3 is true: The aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - If _useNewton3 is false: The aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
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
   * Uses SoAFunctorPairSorted when shouldUseSoASorting() is true; otherwise SoAFunctorPair.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairSoAImpl(ParticleCell_T &cell1, ParticleCell_T &cell2,
                              const std::array<double, 3> &sortingDirection);

  ParticleFunctor_T &_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start AoS sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  /**
   * Min. number of particles to start SoA sorting. This is the sum of the SoA buffer sizes of two cells.
   */
  size_t _soaSortingThreshold{8};

  const DataLayoutOption::Value _dataLayout;

  const bool _useNewton3;
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSoASortingThreshold(size_t soaSortingThreshold) {
  _soaSortingThreshold = soaSortingThreshold;
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
    processCellPairSoAImpl(cell1, cell2, sortingDirection);
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
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  if constexpr (requires {
                  ParticleFunctor_T::getNeededAttr();
                  ParticleFunctor_T::getComputedAttr();
                }) {
    if (shouldUseSoASorting(cell1._particleSoABuffer.size() + cell2._particleSoABuffer.size(), sortingDirection)) {
      using Particle_T = typename ParticleCell_T::ParticleType;
      SoASortedView<Particle_T, ParticleFunctor_T> view1(cell1._particleSoABuffer, sortingDirection);
      SoASortedView<Particle_T, ParticleFunctor_T> view2(cell2._particleSoABuffer, sortingDirection);

      const auto &projIdx1 = view1.projIdx;
      const auto &projIdx2 = view2.projIdx;
      const size_t n1 = view1.getView().size();
      const size_t n2 = view2.getView().size();

      {
        const double threshold = projIdx2[0].first - _sortingCutoff;
        auto it = std::upper_bound(projIdx1.begin(), projIdx1.end(), threshold,
                                   [](double val, const auto &elem) { return val < elem.first; });
        const size_t start_i = static_cast<size_t>(it - projIdx1.begin());

        std::vector<size_t> maxIndex(n1, 0);
        std::vector<size_t> minIndex(n1, 0);
        size_t jUpper = 0, jLower = 0;
        for (size_t i = start_i; i < n1; ++i) {
          while (jUpper < n2 and projIdx2[jUpper].first <= projIdx1[i].first + _sortingCutoff) ++jUpper;
          maxIndex[i] = jUpper;
          while (jLower < n2 and projIdx2[jLower].first < projIdx1[i].first - _sortingCutoff) ++jLower;
          minIndex[i] = jLower;
        }
        _functor.SoAFunctorPairSorted(view1.getView(), view2.getView(), {start_i, maxIndex, minIndex}, _useNewton3);
      }

      if constexpr (bidirectional) {
        if (not _useNewton3) {
          const double threshold = projIdx1[0].first - _sortingCutoff;
          auto it = std::upper_bound(projIdx2.begin(), projIdx2.end(), threshold,
                                     [](double val, const auto &elem) { return val < elem.first; });
          const size_t start_i = static_cast<size_t>(it - projIdx2.begin());

          std::vector<size_t> maxIndex(n2, 0);
          std::vector<size_t> minIndex(n2, 0);
          size_t jUpper = 0, jLower = 0;
          for (size_t i = start_i; i < n2; ++i) {
            while (jUpper < n1 and projIdx1[jUpper].first <= projIdx2[i].first + _sortingCutoff) ++jUpper;
            maxIndex[i] = jUpper;
            while (jLower < n1 and projIdx1[jLower].first < projIdx2[i].first - _sortingCutoff) ++jLower;
            minIndex[i] = jLower;
          }
          _functor.SoAFunctorPairSorted(view2.getView(), view1.getView(), {start_i, maxIndex, minIndex}, false);
        }
      }

      view1.scatterBack();
      view2.scatterBack();
      return;
    }
  }
  _functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, _useNewton3);
  if constexpr (bidirectional) {
    if (not _useNewton3) {
      _functor.SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
    }
  }
}
}  // namespace autopas::internal
