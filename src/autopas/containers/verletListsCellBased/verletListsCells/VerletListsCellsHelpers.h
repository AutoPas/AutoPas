/**
 * @file VerletListsCellsHelpers.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * class of helpers for verlet lists
 * @tparam Particle
 */
template <class ParticleCell>
class VerletListsCellsHelpers {
 public:
  using Particle = typename ParticleCell::ParticleType;
  /// Verlet list storage
  using VerletList_storage_type = std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>;

  /// using declaration for verlet-list particle cell type
  using VerletListParticleCellType = ParticleCell;

  /**
   * This functor can generate verlet lists using the typical pairwise traversal.
   */
  class VerletListGeneratorFunctor : public Functor<Particle> {
   public:
    /**
     * Constructor
     * @param verletLists a verletlist for each cell
     * @param cellMap used to get the verletlist of a particle
     * @param cutoffskin cutoff + skin
     */
    VerletListGeneratorFunctor(VerletList_storage_type &verletLists,
                               std::unordered_map<Particle *, std::pair<size_t, size_t>> &cellMap, double cutoffskin)
        : Functor<Particle>(0.),
          _verletLists(verletLists),
          _cellMap(cellMap),
          _cutoffskinsquared(cutoffskin * cutoffskin) {}

    bool isRelevantForTuning() override { return false; }

    bool allowsNewton3() override {
      utils::ExceptionHandler::exception(
          "VerletListGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool allowsNonNewton3() override {
      utils::ExceptionHandler::exception(
          "VerletListGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
      return true;
    }

    bool isAppropriateClusterSize(unsigned int clusterSize, DataLayoutOption::Value dataLayout) const override {
      return false;  // this functor shouldn't be called with clusters!
    }

    void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
      if (i.isDummy() or j.isDummy()) {
        return;
      }
      auto dist = utils::ArrayMath::sub(i.getR(), j.getR());
      double distsquare = utils::ArrayMath::dot(dist, dist);
      if (distsquare < _cutoffskinsquared) {
        // this is thread safe, only if particle i is accessed by only one
        // thread at a time. which is ensured, as particle i resides in a
        // specific cell and each cell is only accessed by one thread at a time
        // (ensured by traversals)
        // also the list is not allowed to be resized!

        auto indices = _cellMap[&i];
        _verletLists[indices.first][indices.second].second.push_back(&j);
      }
    }

   private:
    VerletList_storage_type &_verletLists;
    std::unordered_map<Particle *, std::pair<size_t, size_t>> &_cellMap;
    double _cutoffskinsquared;
  };

};  // class VerletListsCellsHelpers
}  // namespace autopas
