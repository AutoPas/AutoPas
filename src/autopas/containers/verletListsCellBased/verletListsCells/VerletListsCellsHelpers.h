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
   * Pairwise verlet lists: For every cell and for each of its neighboring cells a pair of particle and a vector of its potential partners is stored.
   */
  using PairwiseNeighborListsType = std::vector<std::vector<std::vector<std::pair<Particle *, std::vector<Particle *>>>>>;
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
        auto iter = _globalToLocalIndex[cellIndex].find(cellIndexNeighbor);
        if(iter == _globalToLocalIndex[cellIndex].end())
        {
          _globalToLocalIndex[cellIndex].insert({cellIndexNeighbor, _globalToLocalIndex[cellIndex].size()});
        }
        auto secondCellIndexInFirst = _globalToLocalIndex[cellIndex].at(cellIndexNeighbor);
        _neighborLists[cellIndex][secondCellIndexInFirst][particleIndex].second.push_back(&j);
      }
    }

   private:
    /**
     * For every cell and for every neghiboring cell of that cell, a vector of pairs.
     * Each pair maps a particle to a vector of its neighbors.
     */
    PairwiseNeighborListsType &_neighborLists;
    std::unordered_map<Particle *, std::pair<size_t, size_t>> &_particleToCellMap;
    //vector of cell1s, for each of theme a map from global cell2 index to
    std::vector<std::unordered_map<size_t, size_t>> _globalToLocalIndex;
    double _cutoffskinsquared;
  };

};  // class VerletListsCellsHelpers
}  // namespace autopas
