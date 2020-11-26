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
template <class Particle>
class VerletListsCellsHelpers {
 public:
  /**
   * Cell wise verlet lists: For every cell, a vector of pairs. Each pair maps a particle to a vector of its neighbors.
   */
  using NeighborListsType = std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>;

  /**
   * Pairwise verlet lists: For every cell (cell1) and for each of its neighboring cells (cell2) a vector of pairs (a
   * pair for each particle). The pairs consist of a particle from cell1 and a vector of its (potential) partners from
   * cell2.
   */
  using PairwiseNeighborListsType =
      std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>;

  /**
   * This functor can generate verlet lists using the typical pairwise traversal.
   */
  class VerletListGeneratorFunctor : public Functor<Particle, VerletListGeneratorFunctor> {
   public:
    /**
     * Constructor
     * @param neighborLists a verletlist for each cell
     * @param particleToCellMap used to get the verletlist of a particle
     * @param cutoffskin cutoff + skin
     */
    VerletListGeneratorFunctor(NeighborListsType &neighborLists,
                               std::unordered_map<Particle *, std::pair<size_t, size_t>> &particleToCellMap,
                               double cutoffskin)
        : Functor<Particle, VerletListGeneratorFunctor>(0.),
          _neighborLists(neighborLists),
          _particleToCellMap(particleToCellMap),
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

    /**
     * @copydoc Functor::AoSFunctor()
     */
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

        auto &[cellIndex, particleIndex] = _particleToCellMap[&i];
        _neighborLists[cellIndex][particleIndex].second.push_back(&j);
      }
    }

   private:
    /**
     * For every cell, a vector of pairs. Each pair maps a particle to a vector of its neighbors.
     */
    NeighborListsType &_neighborLists;
    std::unordered_map<Particle *, std::pair<size_t, size_t>> &_particleToCellMap;
    double _cutoffskinsquared;
  };

  /**
   * This functor generates pairwise verlet lists (a verlet list attached to every pair of neighboring cells).
   */
  class PairwiseVerletListGeneratorFunctor : public Functor<Particle, PairwiseVerletListGeneratorFunctor> {
   public:
    /**
     * Constructor
     * @param neighborLists a verletlist for each cell
     * @param particleToCellMap used to get the verletlist of a particle
     * @param globalToLocalIndex mapping global index of cell2 to "local" index according to cell1's interactions
     * @param cutoffskin cutoff + skin
     */
    PairwiseVerletListGeneratorFunctor(PairwiseNeighborListsType &neighborLists,
                                       std::unordered_map<Particle *, std::pair<size_t, size_t>> &particleToCellMap,
                                       std::vector<std::unordered_map<size_t, size_t>> globalToLocalIndex,
                                       double cutoffskin)
        : Functor<Particle, PairwiseVerletListGeneratorFunctor>(0.),
          _neighborLists(neighborLists),
          _particleToCellMap(particleToCellMap),
          _globalToLocalIndex(globalToLocalIndex),
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

    /**
     * @copydoc Functor::AoSFunctor()
     */
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

        auto &[cellIndex, particleIndex] = _particleToCellMap[&i];
        auto &[cellIndexNeighbor, particleIndexNeighbor] = _particleToCellMap[&j];

        // if cell1 hasn't interacted with cell2 yet and there is no mapping from global to relative index for cell2,
        // add one
        auto iter = _globalToLocalIndex[cellIndex].find(cellIndexNeighbor);
        if (iter == _globalToLocalIndex[cellIndex].end()) {
          _globalToLocalIndex[cellIndex].insert({cellIndexNeighbor, _globalToLocalIndex[cellIndex].size()});
        }
        auto secondCellIndexInFirst = _globalToLocalIndex[cellIndex].at(cellIndexNeighbor);
        _neighborLists[cellIndex][secondCellIndexInFirst][particleIndex].second.push_back(&j);
      }
    }

   private:
    /**
     * Pairwise verlet lists: For every cell (cell1) and for each of its neighboring cells (cell2) a vector of pairs (a
     * pair for each particle). The pairs consist of a particle from cell1 and a vector of its (potential) partners from
     * cell2.
     */
    PairwiseNeighborListsType &_neighborLists;
    /**
     * Mapping of each particle to its corresponding cell and id within this cell.
     */
    std::unordered_map<Particle *, std::pair<size_t, size_t>> &_particleToCellMap;
    /**
     * For each cell1: a mapping of the "absolute" index of cell2 (in the base linked cells structure) to its "relative"
     * index (0 to 26) in cell1's neighbors.
     */
    std::vector<std::unordered_map<size_t, size_t>> _globalToLocalIndex;
    double _cutoffskinsquared;
  };

};  // class VerletListsCellsHelpers
}  // namespace autopas
