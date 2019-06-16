/**
 * @file TraversalVerlet.h
 *
 * @date 7.4.2019
 * @author jspahl
 */

#pragma once

#include "VerletTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides a Traversal for the verlet lists container.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class TraversalVerlet
    : public CellPairTraversal<ParticleCell>,
      public VerletTraversalInterface<
          typename VerletListHelpers<typename ParticleCell::ParticleType>::VerletListParticleCellType> {
  using Particle = typename ParticleCell::ParticleType;
  // using LinkedParticleCell = typename VerletListHelpers<typename Particle>::VerletListParticleCellType ;
  typedef
      typename VerletListHelpers<typename ParticleCell::ParticleType>::VerletListParticleCellType LinkedParticleCell;

 public:
  DataLayoutOption getDataLayout() override { return DataLayout; }

  /**
   * Constructor for Verlet Traversal
   * @param dims dimensions of the underlying container
   * @param pairwiseFunctor Functor to be used with this Traversal
   */
  TraversalVerlet(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims), _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() override { return TraversalOption::verletTraversal; }

  bool isApplicable() override { return DataLayout == DataLayoutOption::aos || DataLayout == DataLayoutOption::soa; }

  void initTraversal(std::vector<ParticleCell> &cells) override {
    if (DataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoALoader(cell, _soa, offset);
        offset += cell.numParticles();
      }
    }
  }

  void endTraversal(std::vector<ParticleCell> &cells) override {
    if (DataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoAExtractor(cell, _soa, offset);
        offset += cell.numParticles();
      }
    }
  }

  /**
   * Initializes Traversal and copies data to this traversal soa storage
   * @param cells content of the container the Traversal is to be called on
   */
  void initTraversal(std::vector<LinkedParticleCell> &cells) override {
    if (DataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoALoader(cell, _soa, offset);
        offset += cell.numParticles();
      }
    }
  }

  /**
   * Ends Traversal and writes data from this Traversals soa back to the cells
   * @param cells content of the container the Traversal is to be called on
   */
  void endTraversal(std::vector<LinkedParticleCell> &cells) override {
    if (DataLayout == DataLayoutOption::soa) {
      size_t offset = 0;
      for (auto &cell : cells) {
        _functor->SoAExtractor(cell, _soa, offset);
        offset += cell.numParticles();
      }
    }
  }

  /**
   * Iterates over the Particles as specified in the Neighbor lists
   * @param aosNeighborLists neighbor lists in aos format
   * @param soaNeighborLists neighbor lists as index list for the soa format
   */
  void iterateVerletLists(
      std::unordered_map<Particle *, std::vector<Particle *>> aosNeighborLists,
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> soaNeighborLists) override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
#if defined(AUTOPAS_OPENMP)
        if (not useNewton3) {
          size_t buckets = aosNeighborLists.bucket_count();
          // @todo find a sensible chunk size
#pragma omp parallel for schedule(dynamic)
          for (size_t b = 0; b < buckets; b++) {
            auto endIter = aosNeighborLists.end(b);
            for (auto it = aosNeighborLists.begin(b); it != endIter; ++it) {
              Particle &i = *(it->first);
              for (auto j_ptr : it->second) {
                Particle &j = *j_ptr;
                _functor->AoSFunctor(i, j, false);
              }
            }
          }
        } else
#endif
        {
          for (auto &list : aosNeighborLists) {
            Particle &i = *list.first;
            for (auto j_ptr : list.second) {
              Particle &j = *j_ptr;
              _functor->AoSFunctor(i, j, useNewton3);
            }
          }
        }
        return;
      }

      case DataLayoutOption::soa: {
        const size_t iFrom = 0;
        const size_t iTo = soaNeighborLists.size();

#if defined(AUTOPAS_OPENMP)
        if (not useNewton3) {
          // @todo find a sensible chunk size
          const size_t chunkSize = std::max((iTo - iFrom) / (omp_get_max_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
          for (size_t i = iFrom; i < iTo; i++) {
            _functor->SoAFunctor(_soa, soaNeighborLists, i, i + 1, useNewton3);
          }
        } else
#endif
        {
          // iterate over SoA
          _functor->SoAFunctor(_soa, soaNeighborLists, iFrom, iTo, useNewton3);
        }
        return;
      }
      default: { utils::ExceptionHandler::exception("VerletList DataLayout {} not available", DataLayout); }
    }
  }

 private:
  /**
   * Functor for Traversal
   */
  PairwiseFunctor *_functor;

  /**
   *global SoA of verlet lists
   */
  SoA<typename Particle::SoAArraysType> _soa;
};

}  // namespace autopas