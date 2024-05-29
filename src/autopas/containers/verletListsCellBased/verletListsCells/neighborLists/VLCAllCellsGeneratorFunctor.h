/**
 * @file VLCAllCellsGeneratorFunctor.h
 * @author tirgendetwas
 * @date 05.12.2020
 */
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

#pragma once

namespace autopas {

/**
 * This functor can generate verlet lists using the typical pairwise traversal.
 */
template <class Particle, enum TraversalOption::Value TraversalOptionEnum>
class VLCAllCellsGeneratorFunctor
    : public Functor<Particle, VLCAllCellsGeneratorFunctor<Particle, TraversalOptionEnum>> {
  using NeighborListsType = typename VerletListsCellsHelpers::AllCellsNeighborListsType<Particle>;
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Constructor
   * @param neighborLists a verletlist for each cell
   * @param particleToCellMap used to get the verletlist of a particle
   * @param cutoffSkin cutoff + skin
   * @param newListAllocationSize Default allocation size for new neighbor lists.
   * @param cellsPerDim Cells per dimension of the underlying Linked Cells container (incl. Halo).
   */
  VLCAllCellsGeneratorFunctor(NeighborListsType &neighborLists,
                              std::unordered_map<Particle *, std::pair<size_t, size_t>> &particleToCellMap,
                              double cutoffSkin, size_t newListAllocationSize,
                              const std::array<size_t, 3> &cellsPerDim = {})
      : Functor<Particle, VLCAllCellsGeneratorFunctor<Particle, TraversalOptionEnum>>(0.),
        _neighborLists(neighborLists),
        _particleToCellMap(particleToCellMap),
        _cutoffSkinSquared(cutoffSkin * cutoffSkin),
        _cellsPerDim(cellsPerDim),
        _newListAllocationSize(newListAllocationSize) {}

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

  /**
   * @copydoc Functor::AoSFunctor()
   */
  void AoSFunctor(Particle &i, Particle &j, bool) override {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const auto dist = i.getR() - j.getR();
    const double distSquare = utils::ArrayMath::dot(dist, dist);
    if (distSquare < _cutoffSkinSquared) {
      // this is thread safe, only if particle i is accessed by only one
      // thread at a time. which is ensured, as particle i resides in a
      // specific cell and each cell is only accessed by one thread at a time
      // (ensured by traversals)
      // also the list is not allowed to be resized!

      // Depending on the targeted traversal neighbor lists are built differently
      if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c08) {
        const auto &[cellIndexI, particleIndexI] = _particleToCellMap[&i];
        const auto &[cellIndexJ, particleIndexJ] = _particleToCellMap[&j];

        const auto cellIndex3DI = utils::ThreeDimensionalMapping::oneToThreeD(cellIndexI, _cellsPerDim);
        const auto cellIndex3DJ = utils::ThreeDimensionalMapping::oneToThreeD(cellIndexJ, _cellsPerDim);

        // WARNING: This is probably only valid for CSF==1
        const std::array<size_t, 3> cellIndex3DBaseCell = utils::ArrayMath::min(cellIndex3DI, cellIndex3DJ);
        const auto cellIndexBaseCell = utils::ThreeDimensionalMapping::threeToOneD(cellIndex3DBaseCell, _cellsPerDim);
        // If the first cell is also the base cell insert regularly
        if (cellIndexBaseCell == cellIndexI) {
          _neighborLists[cellIndexI][particleIndexI].second.push_back(&j);
        } else {
          // In the following two cases the list has to be in a different cell than cellIndex1:
          // 1. If the base cell is the cell of particle j
          //    To respect newton3 == disabled, we can't just do currentList[j].second.push_back(ptr1Ptr[i]).
          // 2. If base cell is neither cellIndex1 or cellIndex2
          auto &list = _neighborLists[cellIndexI];
          auto iterIandList =
              std::find_if(list.begin(), list.end(), [&](const std::pair<Particle *, std::vector<Particle *>> &pair) {
                const auto &[particle, neighbors] = pair;
                return particle == &i;
              });
          if (iterIandList != list.end()) {
            auto &[_, neighbors] = *iterIandList;
            neighbors.push_back(&j);
          } else {
            list.emplace_back(&i, std::vector<Particle *>{});
            list.back().second.reserve(_newListAllocationSize);
            list.back().second.push_back(&j);
          }
        }
      } else if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c01 or
                           TraversalOptionEnum == TraversalOption::Value::vlc_c18) {
        const auto &[cellIndex, particleIndex] = _particleToCellMap[&i];
        _neighborLists[cellIndex][particleIndex].second.push_back(&j);
      } else {
        utils::ExceptionHandler::exception(
            "VLCAllCellsGeneratorFunctor::AoSFunctor(): Encountered incompatible Traversal Option {}.",
            TraversalOptionEnum);
      }
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle()
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    if (soa.size() == 0) return;

    auto **const __restrict__ ptrPtr = soa.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict__ xPtr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yPtr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    // index of cellIndex is particleToCellMap of ptrPtr, same for 2
    const auto [cellIndex, _] = _particleToCellMap.at(ptrPtr[0]);

    // iterating over particle indices and accessing currentList at index i works
    // because the particles are iterated in the same order they are loaded in
    // which is the same order they were initialized when building the aosNeighborList
    const size_t numPart = soa.size();
    for (unsigned int i = 0; i < numPart; ++i) {
      for (unsigned int j = i + 1; j < numPart; ++j) {
        const double drx = xPtr[i] - xPtr[j];
        const double dry = yPtr[i] - yPtr[j];
        const double drz = zPtr[i] - zPtr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _cutoffSkinSquared) {
          auto &currentList = _neighborLists[cellIndex];
          currentList[i].second.push_back(ptrPtr[j]);
          if (not newton3) {
            // we need this here, as SoAFunctorSingle will only be called once for both newton3=true and false.
            currentList[j].second.push_back(ptrPtr[i]);
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
    if (soa1.size() == 0 or soa2.size() == 0) return;

    auto **const __restrict__ ptr1Ptr = soa1.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict__ x1Ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y1Ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z1Ptr = soa1.template begin<Particle::AttributeNames::posZ>();

    auto **const __restrict__ ptr2Ptr = soa2.template begin<Particle::AttributeNames::ptr>();
    double *const __restrict__ x2Ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y2Ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z2Ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    // index of cell1 is particleToCellMap of ptr1Ptr, same for 2
    const size_t cellIndex1 = _particleToCellMap.at(ptr1Ptr[0]).first;
    const size_t cellIndex2 = _particleToCellMap.at(ptr2Ptr[0]).first;

    const auto cellIndexBaseCell = [&]() {
      if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c01 or
                    TraversalOptionEnum == TraversalOption::Value::vlc_c18) {
        return cellIndex1;
      } else if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c08) {
        const auto cellIndex3D1 = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex1, _cellsPerDim);
        const auto cellIndex3D2 = utils::ThreeDimensionalMapping::oneToThreeD(cellIndex2, _cellsPerDim);
        const auto cellIndex3DBaseCell = utils::ArrayMath::min(cellIndex3D1, cellIndex3D2);
        return utils::ThreeDimensionalMapping::threeToOneD(cellIndex3DBaseCell, _cellsPerDim);
      } else {
        utils::ExceptionHandler::exception(
            "VLCAllCellsGeneratorFunctor::SoAFunctor(): Encountered incompatible Traversal Option {}.",
            TraversalOptionEnum);
        return 0ul;
      }
    }();

    auto &currentList = _neighborLists[cellIndexBaseCell];

    // iterating over particle indices and accessing currentList at index i works
    // because the particles are iterated in the same order they are loaded in
    // which is the same order they were initialized when building the aosNeighborList
    const size_t numPart1 = soa1.size();
    for (unsigned int i = 0; i < numPart1; ++i) {
      const size_t numPart2 = soa2.size();
      for (unsigned int j = 0; j < numPart2; ++j) {
        const double drx = x1Ptr[i] - x2Ptr[j];
        const double dry = y1Ptr[i] - y2Ptr[j];
        const double drz = z1Ptr[i] - z2Ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 < _cutoffSkinSquared) {
          if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c01 or
                        TraversalOptionEnum == TraversalOption::Value::vlc_c18) {
            currentList[i].second.push_back(ptr2Ptr[j]);
          } else if constexpr (TraversalOptionEnum == TraversalOption::Value::vlc_c08) {
            // depending on which cell is the base cell append the pointer to the respective list collection
            if (cellIndex1 == cellIndexBaseCell) {
              currentList[i].second.push_back(ptr2Ptr[j]);
            } else {
              // In the following two cases the list has to be in a different cell than cellIndex1:
              // 1. If the base cell is the cell of particle j
              //    To respect newton3 == disabled, we can't just do currentList[j].second.push_back(ptr1Ptr[i]).
              // 2. If base cell is neither cellIndex1 or cellIndex2
              auto iter = std::find_if(currentList.begin(), currentList.end(), [&](const auto &pair) {
                const auto &[particlePtr, list] = pair;
                return particlePtr == ptr1Ptr[i];
              });
              if (iter != currentList.end()) {
                iter->second.push_back(ptr2Ptr[j]);
              } else {
                currentList.emplace_back(ptr1Ptr[i], std::vector<Particle *>{});
                currentList.back().second.reserve(_newListAllocationSize);
                currentList.back().second.push_back(ptr2Ptr[j]);
              }
            }
          }
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
   * For every cell, a vector of pairs. Each pair maps a particle to a vector of its neighbors.
   */
  NeighborListsType &_neighborLists;
  std::unordered_map<Particle *, std::pair<size_t, size_t>> &_particleToCellMap;
  double _cutoffSkinSquared;
  /**
   * Cells per dimension of the underlying Linked Cells container.
   * Needed to calculate the base cell for vlc_c08.
   */
  std::array<size_t, 3> _cellsPerDim;

  /**
   * How many fields are reserved for newly inserted neighbor lists.
   */
  size_t _newListAllocationSize;
};

}  // namespace autopas
