/**
 * @file VerletListsCellsHelpers.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include <atomic>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/SoA.h"
namespace autopas {

/**
 * class of helpers for verlet lists
 * @tparam Particle
 */
template <class Particle>
class VerletListsCellsHelpers {
 public:
  /// Verlet list storage
  typedef std::vector<std::vector<std::pair<Particle*, std::vector<Particle*>>>> VerletList_storage_type;

  /// typedef for verlet-list particle cell type
  typedef FullParticleCell<Particle> VerletListParticleCellType;

  /**
   * This functor can generate verlet lists using the typical pairwise
   * traversal.
   */
  class VerletListGeneratorFunctor : public Functor<Particle, VerletListParticleCellType> {
    typedef VerletListParticleCellType ParticleCell;

   public:
    /**
     * Constructor
     * @param verletLists a verletlist for each cell
     * @param cellMap used to get the verletlist of a particle
     * @param cutoffskin
     */
    VerletListGeneratorFunctor(VerletList_storage_type& verletLists,
                               std::unordered_map<Particle*, std::pair<size_t, size_t>>& cellMap, double cutoffskin)
        : _verletLists(verletLists), _cellMap(cellMap), _cutoffskinsquared(cutoffskin * cutoffskin) {}

    void AoSFunctor(Particle& i, Particle& j, bool newton3) override {
      auto dist = ArrayMath::sub(i.getR(), j.getR());
      double distsquare = ArrayMath::dot(dist, dist);
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
    VerletList_storage_type& _verletLists;
    std::unordered_map<Particle*, std::pair<size_t, size_t>>& _cellMap;
    double _cutoffskinsquared;
  };

};  // class VerletListsCellsHelpers
}  // namespace autopas
