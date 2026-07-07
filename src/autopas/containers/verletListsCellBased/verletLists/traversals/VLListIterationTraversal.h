/**
 * @file VLListIterationTraversal.h
 *
 * @date 7.4.2019
 * @author jspahl
 */
#pragma once

#include "VLTraversalInterface.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class VLListIterationTraversal : public TraversalInterface, public VLTraversalInterface<ParticleCell> {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for Verlet Traversal
   * @param pairwiseFunctor Functor to be used with this Traversal
   * @param dataLayout
   * @param useNewton3
   */
  explicit VLListIterationTraversal(PairwiseFunctor &pairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _functor(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vl_list_iteration; }

  /**
   * VL List iteration is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  void initTraversal() override {
    auto &cells = *(this->_cells);

    // Build a flat index→pointer array so the AoS path can resolve CRS indices
    // back to particle references without a hash-map lookup per neighbor.
    // The iteration order must match buildParticleIndex() in VerletLists exactly:
    // iterate over all cells in order, all particles within each cell in order.
    _particlePtrs.clear();
    _particlePtrs.reserve(this->_neighborList->size());
    for (auto &cell : cells) {
      for (auto &p : cell) {
        _particlePtrs.push_back(&p);
      }
    }

    if (_dataLayout == DataLayoutOption::soa) {
      // Pre-compute per-cell offsets so each SoALoader call can be parallelized
      // independently without knowing the preceding cells' sizes.
      std::vector<size_t> offsets(cells.size() + 1);
      std::inclusive_scan(
          cells.begin(), cells.end(), offsets.begin() + 1,
          [](const size_t &partialSum, const auto &cell) { return partialSum + cell.size(); }, 0);

      _soa.resizeArrays(offsets.back());

      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        _functor.SoALoader(cells[i], _soa, offsets[i], /*skipSoAResize*/ true);
      }
    }
  }

  void endTraversal() override {
    auto &cells = *(this->_cells);
    if (_dataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor.SoAExtractor(cell, _soa, offset);
        offset += cell.size();
      }
    }
  }

  void traverseParticles() override {
    auto &neighborList = *(this->_neighborList);
    const size_t numParticles = neighborList.size();

    switch (this->_dataLayout) {
      case DataLayoutOption::aos: {
        if (not _useNewton3) {
          // Each particle i owns its own list slice — no write conflict between iterations.
          AUTOPAS_OPENMP(parallel for schedule(static))
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *_particlePtrs[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *_particlePtrs[neighborsIPtr[j]], false);
            }
          }
        } else {
          // @todo: VL parallelization with N3 should be implemented. Requires a coloring scheme or similar.
          for (size_t i = 0; i < numParticles; ++i) {
            ParticleType &particleI = *_particlePtrs[i];
            const size_t numNeighbors = neighborList.count(i);
            const size_t *neighborsIPtr = neighborList.begin(i);
            for (size_t j = 0; j < numNeighbors; ++j) {
              _functor.AoSFunctor(particleI, *_particlePtrs[neighborsIPtr[j]], true);
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        // Pass the CRS slice directly as a raw pointer + count.
        // LJFunctor (and any functor overriding the pointer overload) receives
        // nl.begin(i) directly — zero allocation, zero copy per particle.
        // Functors that only override the vector overload fall back to the default
        // shim in PairwiseFunctor which builds a temporary vector.
        if (not _useNewton3) {
          AUTOPAS_OPENMP(parallel for schedule(dynamic, std::max(numParticles / (autopas::autopas_get_max_threads() * 10), 1ul)))
          for (size_t i = 0; i < numParticles; ++i) {
            _functor.SoAFunctorVerlet(_soa, i, neighborList.begin(i), neighborList.count(i), false);
          }
        } else {
          for (size_t i = 0; i < numParticles; ++i) {
            _functor.SoAFunctorVerlet(_soa, i, neighborList.begin(i), neighborList.count(i), true);
          }
        }
        return;
      }
      default: {
        utils::ExceptionHandler::exception("VerletList dataLayout {} not available", _dataLayout);
      }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  PairwiseFunctor &_functor;

  /**
   * SoA buffer of verlet lists.
   */
  SoA<typename ParticleType::SoAArraysType> _soa;

  /**
   * Flat array mapping SoA index i -> pointer to particle i.
   * Built in initTraversal() in the same cell/particle iteration order as VerletLists::buildParticleIndex().
   * Used during AoS force computations to resolve CRS neighbor indices to particles.
   */
  std::vector<ParticleType *> _particlePtrs;
};

}  // namespace autopas