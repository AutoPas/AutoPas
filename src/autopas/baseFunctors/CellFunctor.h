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
#include "autopas/utils/SortedSoAView.h"
#include "autopas/utils/WrapOpenMP.h"

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
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {
    if (dataLayout == DataLayoutOption::soa) {
      _soaThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

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
   * Set the aos-sorting-threshold
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param aosSortingThreshold Sum of the number of particles in two cells from which sorting should be enabled.
   */
  void setAoSSortingThreshold(size_t aosSortingThreshold);

  /**
   * Set the SoA sorting-threshold.
   * If the sum of the number of particles in two SoA buffers is greater than this value, the SoA path uses
   * SoAFunctorPairSorted instead of SoAFunctorPair.
   * @param soaSortingThreshold Sum of the number of particles in two cells from which SoA sorting should be enabled.
   */
  void setSoASortingThreshold(size_t soaSortingThreshold);

  /**
   * Computes the per-particle index bounds into projIdxJ needed by SoAFunctorPairSorted.
   * @param projIdxI Sorted projections of the outer-loop cell.
   * @param projIdxJ Sorted projections of the inner-loop cell.
   * @param maxIndexCache Cache to store the computed maxIndex.
   * @param minIndexCache Cache to store the computed minIndex.
   * @return SoASortedPairMeta with startI and per-i upper/lower bounds into j.
   */
  [[nodiscard]] SoASortingData computeSortingData(const std::vector<std::pair<double, size_t>> &projIdxI,
                                                  const std::vector<std::pair<double, size_t>> &projIdxJ,
                                                  std::vector<size_t> &maxIndexCache,
                                                  std::vector<size_t> &minIndexCache) const;

 private:
  /**
   * Evaluate whether the AoSFunctor should use sorting, depending on the set sorting threshold.
   * @param particleCount Total number of involved particles.
   * @param sortingDirection No sorting when the sorting direction is {0., 0., 0.}.
   * @return whether the AoSFunctor should use the SortedCellView.
   */
  [[nodiscard]] bool shouldUseAoSSorting(size_t particleCount, const std::array<double, 3> &sortingDirection) const {
    return particleCount >= _aosSortingThreshold and
           (sortingDirection[0] != 0.0 or sortingDirection[1] != 0.0 or sortingDirection[2] != 0.0);
  }

  /**
   * Evaluate whether the SoA path should use sorted pair iteration, depending on the set SoA sorting threshold.
   * @param particleCount Total number of involved particles (sum of both SoA buffer sizes).
   * @param sortingDirection No sorting when the sorting direction is {0., 0., 0.}.
   * @return Whether the SoA path should use SoAFunctorPairSorted.
   */
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
  size_t _aosSortingThreshold{8};

  /**
   * Min. number of particles to start SoA sorting. This is the sum of the SoA buffer sizes of two cells.
   */
  size_t _soaSortingThreshold{25};

  const DataLayoutOption::Value _dataLayout;

  const bool _useNewton3;

  /**
   * Per-thread persistent buffers for the SoA sorting path, allocated once and reused across calls.
   */
  struct alignas(64) SoAThreadData {
    SoA<typename ParticleCell_T::ParticleType::SoAArraysType> sortedSoa1;
    SoA<typename ParticleCell_T::ParticleType::SoAArraysType> sortedSoa2;
    std::vector<std::pair<double, size_t>> projIdx1;
    std::vector<std::pair<double, size_t>> projIdx2;
    std::vector<size_t> maxIndex;
    std::vector<size_t> minIndex;
  };

  std::vector<SoAThreadData> _soaThreadData{};
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setAoSSortingThreshold(size_t aosSortingThreshold) {
  _aosSortingThreshold = aosSortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSoASortingThreshold(size_t soaSortingThreshold) {
  _soaSortingThreshold = soaSortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCell(ParticleCell_T &cell) {
  const bool isAoS = _dataLayout == DataLayoutOption::aos;
  const bool isSoA = _dataLayout == DataLayoutOption::soa;

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
  const bool isAoS = _dataLayout == DataLayoutOption::aos;
  const bool isSoA = _dataLayout == DataLayoutOption::soa;

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

  if (cell.size() >= _aosSortingThreshold) {
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

  if (shouldUseAoSSorting(cell1.size() + cell2.size(), sortingDirection)) {
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
SoASortingData CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::computeSortingData(
    const std::vector<std::pair<double, size_t>> &projIdxI, const std::vector<std::pair<double, size_t>> &projIdxJ,
    std::vector<size_t> &maxIndexCache, std::vector<size_t> &minIndexCache) const {
  const size_t nI = projIdxI.size();
  const size_t nJ = projIdxJ.size();

  // Compute startI: the first i-particle that can interact with any j-particle.
  // Any i with projI[i] <= projJ[0] - cutoff is further than cutoff from every j along the sorting axis,
  // so it cannot contribute an interaction and is skipped.
  const double threshold = projIdxJ[0].first - _sortingCutoff;
  auto start_iter = std::upper_bound(projIdxI.begin(), projIdxI.end(), threshold,
                                     [](double val, const auto &elem) { return val < elem.first; });
  const size_t startI = static_cast<size_t>(start_iter - projIdxI.begin());

  // Compute maxIndexCache and minIndexCache in a single O(nI + nJ) sweep.
  // Both bounds are monotonically non-decreasing with i because projIdxI is sorted, so each pointer only
  // advances, this is approximate but avoids a O(n^2) scan.
  // For each i, the j-particles in [minIndex, maxIndex) are the candidates satisfying the 1-D cutoff check.
  maxIndexCache.assign(nI, 0);
  minIndexCache.assign(nI, 0);
  size_t jUpper = 0, jLower = 0;
  for (size_t i = startI; i < nI; ++i) {
    // jUpper: first j where projJ > projI[i] + cutoff (exclusive upper bound).
    while (jUpper < nJ and projIdxJ[jUpper].first <= projIdxI[i].first + _sortingCutoff) ++jUpper;
    maxIndexCache[i] = jUpper;
    // jLower: first j where projJ >= projI[i] - cutoff (lower bound).
    while (jLower < nJ and projIdxJ[jLower].first < projIdxI[i].first - _sortingCutoff) ++jLower;
    minIndexCache[i] = jLower;
  }

  return {startI, maxIndexCache, minIndexCache};
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAImpl(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  if constexpr (ParticleFunctor_T::supportsSoASorting) {
    if (shouldUseSoASorting(cell1._particleSoABuffer.size() + cell2._particleSoABuffer.size(), sortingDirection)) {
      using Particle_T = ParticleCell_T::ParticleType;
      auto &thread_data = _soaThreadData[autopas::autopas_get_thread_num()];

      SortedSoAView<Particle_T, ParticleFunctor_T> view1(cell1._particleSoABuffer, sortingDirection,
                                                         thread_data.sortedSoa1, thread_data.projIdx1);
      SortedSoAView<Particle_T, ParticleFunctor_T> view2(cell2._particleSoABuffer, sortingDirection,
                                                         thread_data.sortedSoa2, thread_data.projIdx2);

      _functor.SoAFunctorPairSorted(
          view1.getView(), view2.getView(),
          computeSortingData(view1.projIdx, view2.projIdx, thread_data.maxIndex, thread_data.minIndex), _useNewton3);

      if constexpr (bidirectional) {
        if (not _useNewton3) {
          _functor.SoAFunctorPairSorted(
              view2.getView(), view1.getView(),
              computeSortingData(view2.projIdx, view1.projIdx, thread_data.maxIndex, thread_data.minIndex), false);
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
