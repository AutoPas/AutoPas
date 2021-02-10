/**
 * @file DSKokkosParallelTraversal.h
 * @author F. Gratl
 * @copier L. Gaertner
 * @date 11/23/18
 */

#pragma once

#include <vector>

#include "DSTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/DataLayoutConverter.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3 = true, bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This paramater indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   */
  explicit CellFunctor(ParticleFunctor *f, const double sortingCutoff) : _functor(f), _sortingCutoff(sortingCutoff) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
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

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairCudaNoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairCudaN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  void processCellCudaNoN3(ParticleCell &cell);

  void processCellCudaN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting.
   * @todo Currently, this is disabled because of https://github.com/AutoPas/AutoPas/issues/418
   */
  constexpr static unsigned long _startSorting = std::numeric_limits<unsigned long>::max();
};

/**
 * This sum defines the traversal typically used by the DirectSum container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class DSKokkosParallelTraversal : public CellPairTraversal<ParticleCell>, public DSTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor for the DirectSum traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the directsum traversal, please don't use the interaction length here.)
   */
  explicit DSKokkosParallelTraversal(PairwiseFunctor *pairwiseFunctor, double cutoff)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ds_sequential; }

  [[nodiscard]] bool isApplicable() const override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    if (dataLayout == DataLayoutOption::cuda)
      return nDevices > 0;
    else
      return true;
  }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    auto &cells = *(this->_cells);
    for (auto &cell : cells) {
      _dataLayoutConverter.loadDataLayout(cell);
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    for (auto &cell : cells) {
      _dataLayoutConverter.storeDataLayout(cell);
    }
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, dataLayout, useNewton3,
                        true>
      _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void DSKokkosParallelTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  // Assume cell[0] is the main domain and cell[1] is the halo
  _cellFunctor.processCell(cells[0]);
  _cellFunctor.processCellPair(cells[0], cells[1]);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((DataLayout == DataLayoutOption::soa && cell._particleSoABuffer.getNumParticles() == 0) ||
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
    case DataLayoutOption::cuda:
      if (useNewton3) {
        processCellCudaN3(cell);
      } else {
        processCellCudaNoN3(cell);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPair(

    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if ((DataLayout == DataLayoutOption::soa &&
       (cell1._particleSoABuffer.getNumParticles() == 0 || cell2._particleSoABuffer.getNumParticles() == 0)) ||
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
    case DataLayoutOption::cuda:
      if (useNewton3) {
        processCellPairCudaN3(cell1, cell2);
      } else {
        processCellPairCudaNoN3(cell1, cell2);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
template <bool newton3>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellAoS(
    ParticleCell &cell) {
  if (cell.numParticles() > _startSorting) {
    SortedCellView<Particle, ParticleCell> cellSorted(
        cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    auto outer = cellSorted._particles.begin();
    for (; outer != cellSorted._particles.end(); ++outer) {
      Particle &p1 = *outer->second;

      auto inner = outer;
      ++inner;
      for (; inner != cellSorted._particles.end(); ++inner) {
        if (std::abs(outer->first - inner->first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *inner->second;
        if constexpr (newton3) {
          _functor->AoSFunctor(p1, p2, true);
        } else {
          _functor->AoSFunctor(p1, p2, false);
          _functor->AoSFunctor(p2, p1, false);
        }
      }
    }
  } else {
    auto outer = getStaticCellIter(cell);
    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      auto inner = outer;
      ++inner;
      for (; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        if constexpr (newton3) {
          _functor->AoSFunctor(p1, p2, true);
        } else {
          _functor->AoSFunctor(p1, p2, false);
          _functor->AoSFunctor(p2, p1, false);
        }
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if (cell1.numParticles() + cell2.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> baseSorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> outerSorted(cell2, sortingDirection);

    for (auto &outer : baseSorted._particles) {
      Particle &p1 = *outer.second;

      for (auto &inner : outerSorted._particles) {
        if (std::abs(outer.first - inner.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *inner.second;
        _functor->AoSFunctor(p1, p2, true);
      }
    }
  } else {
    auto outer = getStaticCellIter(cell1);
    auto innerStart = getStaticCellIter(cell2);

    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      for (auto inner = innerStart; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, true);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
    bidirectional>::processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                           const std::array<double, 3> &sortingDirection) {
  if (cell1.numParticles() + cell2.numParticles() > _startSorting and
      sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> baseSorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> outerSorted(cell2, sortingDirection);

    for (auto &outer : baseSorted._particles) {
      Particle &p1 = *outer.second;

      for (auto &inner : outerSorted._particles) {
        if (std::abs(outer.first - inner.first) > _sortingCutoff) {
          break;
        }
        Particle &p2 = *inner.second;
        _functor->AoSFunctor(p1, p2, false);
        if (bidirectional) _functor->AoSFunctor(p2, p1, false);
      }
    }
  } else {
    auto innerStart = getStaticCellIter(cell2);

    for (auto outer = cell1.begin(); outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      for (auto inner = innerStart; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, false);
        if (bidirectional) _functor->AoSFunctor(p2, p1, false);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairSoAN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
    bidirectional>::processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if (bidirectional) _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoAN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoANoN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
    bidirectional>::processCellPairCudaNoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->CudaFunctor(cell1._particleSoABufferDevice, cell2._particleSoABufferDevice, false);
  if (bidirectional) _functor->CudaFunctor(cell2._particleSoABufferDevice, cell1._particleSoABufferDevice, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellCudaNoN3(
    ParticleCell &cell) {
  _functor->CudaFunctor(cell._particleSoABufferDevice, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairCudaN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->CudaFunctor(cell1._particleSoABufferDevice, cell2._particleSoABufferDevice, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout,
    bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellCudaN3(
    ParticleCell &cell) {
  _functor->CudaFunctor(cell._particleSoABufferDevice, true);
}

}  // namespace autopas
