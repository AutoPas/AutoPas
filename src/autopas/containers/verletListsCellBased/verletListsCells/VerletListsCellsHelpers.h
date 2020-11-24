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

    using SoAArraysType = typename Particle::SoAArraysType;
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


    /**
 * SoAFunctor for verlet list generation. (single cell version)
 * @param soa the soa
 * @param newton3 whether to use newton 3
 */
    void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
      if (soa.getNumParticles() == 0) return;

      auto **const __restrict__ ptrptr = soa.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

      //index of cell1 is particleToCellMap of ptr1ptr, same for 2
      auto cell1 = _particleToCellMap.at(ptrptr[0]).first;

      auto iter = _globalToLocalIndex[cell1].find(cell1);
      if(iter == _globalToLocalIndex[cell1].end())
      {
        _globalToLocalIndex[cell1].insert({cell1, _globalToLocalIndex[cell1].size()});
      }
      auto localCell2Index = _globalToLocalIndex[cell1].at(cell1);

      auto &currentList = _neighborLists[cell1][localCell2Index];


      size_t numPart = soa.getNumParticles();
      for (unsigned int i = 0; i < numPart; ++i) {

        for (unsigned int j = i + 1; j < numPart; ++j) {
          const double drx = xptr[i] - xptr[j];
          const double dry = yptr[i] - yptr[j];
          const double drz = zptr[i] - zptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList[i].second.push_back(ptrptr[j]);
            if (not newton3) {
              // we need this here, as SoAFunctorSingle will only be called once for both newton3=true and false.
              currentList[j].second.push_back(ptrptr[i]);
            }
          }
        }
      }
    }

    /**
     * SoAFunctor for the verlet list generation. (two cell version)
     * @param soa1 soa of first cell
     * @param soa2 soa of second cell
     * @note newton3 is ignored here, as for newton3=false SoAFunctorPair(soa2, soa1) will also be called.
     */
    void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool /*newton3*/) override {
      if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

      auto **const __restrict__ ptr1ptr = soa1.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();

      auto **const __restrict__ ptr2ptr = soa2.template begin<Particle::AttributeNames::ptr>();
      double *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

      //index of cell1 is particleToCellMap of ptr1ptr, same for 2
      size_t cell1 = _particleToCellMap.at(ptr1ptr[0]).first;
      size_t cell2 = _particleToCellMap.at(ptr2ptr[0]).first;

      auto iter = _globalToLocalIndex[cell1].find(cell2);
      if(iter == _globalToLocalIndex[cell1].end())
      {
        _globalToLocalIndex[cell1].insert({cell2, _globalToLocalIndex[cell1].size()});
      }
      auto localCell2Index = _globalToLocalIndex[cell1].at(cell2);

      auto &currentList = _neighborLists[cell1][localCell2Index];

      size_t numPart1 = soa1.getNumParticles();
      //iterate through particles in cell1
      for (unsigned int i = 0; i < numPart1; ++i) {
        size_t numPart2 = soa2.getNumParticles();
        //iterate through particles in cell2
        for (unsigned int j = 0; j < numPart2; ++j) {
          const double drx = x1ptr[i] - x2ptr[j];
          const double dry = y1ptr[i] - y2ptr[j];
          const double drz = z1ptr[i] - z2ptr[j];

          const double drx2 = drx * drx;
          const double dry2 = dry * dry;
          const double drz2 = drz * drz;

          const double dr2 = drx2 + dry2 + drz2;

          if (dr2 < _cutoffskinsquared) {
            currentList[i].second.push_back(ptr2ptr[j]);
            //is i really the index of the particle?
          }
        }
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
