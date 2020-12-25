/**
 * @file VLCCellPairGeneratorFunctor.h
 * @author tirgendetwas
 * @date 05.12.2020
 */

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/pairwiseFunctors/Functor.h"

#pragma once

namespace autopas {

template <class Particle>
/**
 * This functor generates pairwise verlet lists (a verlet list attached to every pair of neighboring cells).
 */
class VLCCellPairGeneratorFunctor : public Functor<Particle, VLCCellPairGeneratorFunctor<Particle>> {
  using PairwiseNeighborListsType = typename VerletListsCellsHelpers<Particle>::PairwiseNeighborListsType;
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Constructor
   * @param neighborLists a verletlist for each cell
   * @param particleToCellMap used to get the verletlist of a particle
   * @param globalToLocalIndex mapping global index of cell2 to "local" index according to cell1's interactions
   * @param cutoffskin cutoff + skin
   */
  VLCCellPairGeneratorFunctor(PairwiseNeighborListsType &neighborLists,
                              std::unordered_map<Particle *, std::pair<size_t, size_t>> &particleToCellMap,
                              std::vector<std::unordered_map<size_t, size_t>> &globalToLocalIndex, double cutoffskin)
      : Functor<Particle, VLCCellPairGeneratorFunctor<Particle>>(0.),
        _neighborLists(neighborLists),
        _particleToCellMap(particleToCellMap),
        _globalToLocalIndex(globalToLocalIndex),
        _cutoffskinsquared(cutoffskin * cutoffskin) {}

  bool isRelevantForTuning() override { return false; }

  bool allowsNewton3() override {
    utils::ExceptionHandler::exception(
        "VLCAllCellsGeneratorFunctor::allowsNewton3() is not implemented, because it should not be called.");
    return true;
  }

  bool allowsNonNewton3() override {
    utils::ExceptionHandler::exception(
        "VLCAllCellsGeneratorFunctor::allowsNonNewton3() is not implemented, because it should not be called.");
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

  /**
   * @copydoc Functor::SoAFunctorSingle()
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    auto **const __restrict__ ptrptr = soa.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    // index of cell1 is particleToCellMap of ptr1ptr, same for 2
    auto cell = _particleToCellMap.at(ptrptr[0]).first;

    auto iter = _globalToLocalIndex[cell].find(cell);
    if (iter == _globalToLocalIndex[cell].end()) {
      _globalToLocalIndex[cell].insert({cell, _globalToLocalIndex[cell].size()});
    }
    auto localCell2Index = _globalToLocalIndex[cell].at(cell);

    auto &currentList = _neighborLists[cell][localCell2Index];

    // iterating over particle indices and accessing currentList at index i works
    // because the particles are iterated in the same order they are loaded in
    // which is the same order they were initialized when building the aosNeighborList
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
   * Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
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

    // index of cell1 is particleToCellMap of ptr1ptr, same for 2
    size_t cell1 = _particleToCellMap.at(ptr1ptr[0]).first;
    size_t cell2 = _particleToCellMap.at(ptr2ptr[0]).first;

    auto iter = _globalToLocalIndex[cell1].find(cell2);
    if (iter == _globalToLocalIndex[cell1].end()) {
      _globalToLocalIndex[cell1].insert({cell2, _globalToLocalIndex[cell1].size()});
    }
    auto localCell2Index = _globalToLocalIndex[cell1].at(cell2);

    auto &currentList = _neighborLists[cell1][localCell2Index];

    // iterating over particle indices and accessing currentList at index i works
    // because the particles are iterated in the same order they are loaded in
    // which is the same order they were initialized when building the aosNeighborList
    size_t numPart1 = soa1.getNumParticles();
    for (unsigned int i = 0; i < numPart1; ++i) {
      size_t numPart2 = soa2.getNumParticles();
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
          // is i really the index of the particle?
        }
      }
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 4>{
        Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 4>{
        Particle::AttributeNames::ptr, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() { return std::array<typename Particle::AttributeNames, 0>{}; }

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
   * index in cell1's neighbors.
   */
  std::vector<std::unordered_map<size_t, size_t>> &_globalToLocalIndex;
  double _cutoffskinsquared;
};

}  // namespace autopas
