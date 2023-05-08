/**
 * @file IteratorTestHelper.h
 * @author F. Gratl
 * @date 05.03.21
 */

#pragma once

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "testingHelpers/commonTypedefs.h"

namespace IteratorTestHelper {

/**
 * Inserts particles around all corners of the given AutoPas object at critical distances.
 * @tparam AutoPasT
 * @param autoPas
 * @param boxOfInterestMin
 * @param boxOfInterestMax
 * @return Tuple of four vectors containing IDs of added particles.
 *   1. All owned particles
 *   2. All halo particles.
 *   3. All owned particles in the box of interest.
 *   3. All halo particles in the box of interest.
 */
template <class AutoPasT>
auto fillContainerAroundBoundary(AutoPasT &autoPas, std::array<double, 3> boxOfInterestMin,
                                 std::array<double, 3> boxOfInterestMax) {
  constexpr size_t numParticles1dTotal = 10;

  auto cutoff = autoPas.getCutoff();
  auto skin = autoPas.getCurrentVerletSkin();

  // generator function for critical coordinates (along  one dimension)
  auto generateInteresting1DPositions = [&](double min, double max) -> auto {
    // ensure that all particles are at most skin away from halo!
    // interesting cases are:
    //   - outside of the halo by skin
    //   - edge of halo
    //   - in the halo
    //   - edge of actual domain
    //   - just inside the domain
    return std::array<double, numParticles1dTotal>{min - cutoff - skin + 1e-10,
                                                   min - cutoff,
                                                   min - skin / 4,
                                                   min,
                                                   min + skin / 4,
                                                   max - skin / 4,
                                                   max,
                                                   max + skin / 4,
                                                   max + cutoff,
                                                   max + cutoff + skin - 1e-10};
  };

  // fill container
  size_t id = 0;
  auto boxMin = autoPas.getBoxMin();
  auto boxMax = autoPas.getBoxMax();

  std::vector<size_t> particleIDsInterestHalo;
  std::vector<size_t> particleIDsInterestOwned;
  std::vector<size_t> particleIDsHalo;
  std::vector<size_t> particleIDsOwned;
  for (auto x : generateInteresting1DPositions(boxMin[0], boxMax[0])) {
    for (auto y : generateInteresting1DPositions(boxMin[1], boxMax[1])) {
      for (auto z : generateInteresting1DPositions(boxMin[2], boxMax[2])) {
        std::array<double, 3> pos{x, y, z};
        Molecule p(pos, {0., 0., 0.}, id++, 0);
        // add the particle as actual or halo particle
        if (autopas::utils::inBox(pos, boxMin, boxMax)) {
          autoPas.addParticle(p);
          particleIDsOwned.push_back(p.getID());
          if (autopas::utils::inBox(pos, boxOfInterestMin, boxOfInterestMax)) {
            particleIDsInterestOwned.push_back(p.getID());
          }
        } else {
          // AutoPas should set the ownership state of this particle to halo
          autoPas.addHaloParticle(p);
          particleIDsHalo.push_back(p.getID());
          if (autopas::utils::inBox(pos, boxOfInterestMin, boxOfInterestMax)) {
            particleIDsInterestHalo.push_back(p.getID());
          }
        }
      }
    }
  }

  // sanity check. Can not use ASSERT_EQ because this introduces a different return.
  EXPECT_EQ(particleIDsOwned.size() + particleIDsHalo.size(),
            numParticles1dTotal * numParticles1dTotal * numParticles1dTotal);
  // getNumberOfParticles works via counters in the logic handler
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), particleIDsOwned.size());
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::halo), particleIDsHalo.size());
  return std::make_tuple(particleIDsOwned, particleIDsHalo, particleIDsInterestOwned, particleIDsInterestHalo);
}

/**
 * Inserts particles around all corners of the given AutoPas object at critical distances.
 * @tparam AutoPasT
 * @param autoPas
 * @return Tuple of two vectors containing IDs of added particles. First for owned, second for halo particles.
 */
template <class AutoPasT>
auto fillContainerAroundBoundary(AutoPasT &autoPas) {
  std::array<double, 3> numericLimitMin{std::numeric_limits<double>::min(), std::numeric_limits<double>::min(),
                                        std::numeric_limits<double>::min()};
  std::array<double, 3> numericLimitMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                        std::numeric_limits<double>::max()};
  return fillContainerAroundBoundary(autoPas, numericLimitMin, numericLimitMax);
}

/**
 * Creats a grid of particles in the given AutoPas object.
 * Grid width is `sparsity * ( boxLength / ((cutoff + skin) * cellSizeFactor) )`.
 * E.g., for a sparsity of 1, 1 particle is inserted for every cell. For a sparsity of .5, 8 particles are inserted.
 * The lower corner of the grid is offset from boxMin by half the grid width in every dimension.
 * This way there should be one particle in every third Linked Cells cell.
 * @tparam AutoPasT
 * @param autoPas
 * @param sparsity
 * @return Vector of all particle IDs added.
 */
template <class AutoPasT>
auto fillContainerWithGrid(AutoPasT &autoPas, double sparsity) {
  auto cutoff = autoPas.getCutoff();
  auto skin = autoPas.getCurrentVerletSkin();
  auto cellSizeFactor = *(autoPas.getAllowedCellSizeFactors().getAll().begin());

  auto boxLength = autopas::utils::ArrayMath::sub(autoPas.getBoxMax(), autoPas.getBoxMin());

  auto gridWidth1D = (cutoff + skin) * cellSizeFactor;
  auto gridEdgesPerDim = autopas::utils::ArrayMath::mulScalar(boxLength, 1 / gridWidth1D);
  auto gridWidth3D = autopas::utils::ArrayMath::div(boxLength, gridEdgesPerDim);

  size_t id = 0;
  std::vector<size_t> particleIDs;
  for (double x = gridWidth3D[0] / 2; x < boxLength[0]; x += sparsity * gridWidth3D[0]) {
    for (double y = gridWidth3D[1] / 2; y < boxLength[1]; y += sparsity * gridWidth3D[1]) {
      for (double z = gridWidth3D[2] / 2; z < boxLength[2]; z += sparsity * gridWidth3D[2]) {
        std::array<double, 3> pos{x, y, z};
        Molecule p(pos, {0., 0., 0.}, id++, 0);
        autoPas.addParticle(p);
        particleIDs.push_back(p.getID());
      }
    }
  }

  return particleIDs;
}

template <class AutoPasT>
auto getHaloBoxMinMax(AutoPasT &autoPas) {
  const auto interactionLength = autoPas.getCutoff() + autoPas.getCurrentVerletSkin();
  // halo has width of interactionLength
  const auto haloBoxMin = autopas::utils::ArrayMath::subScalar(autoPas.getBoxMin(), interactionLength);
  const auto haloBoxMax = autopas::utils::ArrayMath::addScalar(autoPas.getBoxMax(), interactionLength);

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

/**
 * Creates a function to instantiate an iterator with the given properties and passes this function on to fun.
 * The iterator always covers the whole domain and, if necessary the halo.
 * This is necessary so that fun can decide for itself if it wants iterators to be created in an OpenMP region or not.
 * @tparam AutoPasT
 * @tparam F f(AutoPas, Iterator)
 * @param useRegionIterator
 * @param useConstIterator
 * @param behavior
 * @param autoPas
 * @param fun Function taking the AutoPas object and the generated iterator.
 */
template <bool useConstIterator, class AutoPasT, class F>
void provideIterator(AutoPasT &autoPas, autopas::IteratorBehavior behavior, bool useRegionIterator, F fun) {
  if (useRegionIterator) {
    std::array<double, 3> haloBoxMin, haloBoxMax;
    std::tie(haloBoxMin, haloBoxMax) = getHaloBoxMinMax(autoPas);
    if constexpr (useConstIterator) {
      const auto &autoPasRef = autoPas;
      auto getIter = [&]() -> typename AutoPasT::const_iterator_t {
        return autoPasRef.getRegionIterator(haloBoxMin, haloBoxMax, behavior);
      };
      fun(autoPasRef, getIter);
    } else {
      auto getIter = [&]() ->
          typename AutoPasT::iterator_t { return autoPas.getRegionIterator(haloBoxMin, haloBoxMax, behavior); };
      fun(autoPas, getIter);
    }
  } else {
    if constexpr (useConstIterator) {
      auto getIter = [&]() -> typename AutoPasT::const_iterator_t { return autoPas.cbegin(behavior); };
      fun(autoPas, getIter);
    } else {
      auto getIter = [&]() -> typename AutoPasT::iterator_t { return autoPas.begin(behavior); };
      fun(autoPas, getIter);
    }
  }
}

/**
 * Same as provideIterator(), but `useConstIterator` is a run-time variable.
 * @tparam useConstIterator
 * @tparam AutoPasT
 * @tparam F f(AutoPas, Iterator)
 * @param useRegionIterator
 * @param behavior
 * @param autoPas
 * @param fun Function taking the AutoPas object and the generated iterator.
 */
template <class AutoPasT, class F>
void provideIterator(bool useConstIterator, AutoPasT &autoPas, autopas::IteratorBehavior behavior,
                     bool useRegionIterator, F fun) {
  if (useConstIterator) {
    provideIterator<true>(autoPas, behavior, useRegionIterator, fun);
  } else {
    provideIterator<false>(autoPas, behavior, useRegionIterator, fun);
  }
}

/**
 * Creates a function to instantiate a region-iterator with the given properties and passes this function on to fun.
 * This is necessary so that fun can decide for itself if it wants iterators to be created in an OpenMP region or not.
 * @tparam useConstIterator
 * @tparam AutoPasT
 * @tparam F f(AutoPas, Iterator)
 * @param autoPas
 * @param behavior
 * @param boxMin
 * @param boxMax
 * @param fun Function taking the AutoPas object and the generated iterator.
 */
template <bool useConstIterator, class AutoPasT, class F>
void provideRegionIterator(AutoPasT &autoPas, autopas::IteratorBehavior behavior, const std::array<double, 3> &boxMin,
                           const std::array<double, 3> &boxMax, F fun) {
  if constexpr (useConstIterator) {
    const auto &autoPasRef = autoPas;
    auto getIter = [&]() ->
        typename AutoPasT::const_iterator_t { return autoPasRef.getRegionIterator(boxMin, boxMax, behavior); };
    fun(autoPasRef, getIter);
  } else {
    auto getIter = [&]() ->
        typename AutoPasT::iterator_t { return autoPas.getRegionIterator(boxMin, boxMax, behavior); };
    fun(autoPas, getIter);
  }
}

/**
 * Same as provideRegionIterator() but `useConstIterator` is a run-time variable.
 * @tparam AutoPasT
 * @tparam F f(AutoPas, Iterator)
 * @param useConstIterator
 * @param autoPas
 * @param behavior
 * @param boxMin
 * @param boxMax
 * @param fun Function taking the AutoPas object and the generated iterator.
 */
template <class AutoPasT, class F>
void provideRegionIterator(bool useConstIterator, AutoPasT &autoPas, autopas::IteratorBehavior behavior,
                           const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, F fun) {
  if (useConstIterator) {
    provideRegionIterator<true>(autoPas, behavior, boxMin, boxMax, fun);
  } else {
    provideRegionIterator<false>(autoPas, behavior, boxMin, boxMax, fun);
  }
}

/**
 * Apply an iterator, track what particle IDs are found and compare this to a vector of expected IDs
 * @tparam AutoPasT
 * @tparam IteratorT
 * @param autopas
 * @param iterator
 * @param particleIDsExpected
 */
template <class AutoPasT, class FgetIter>
void findParticles(AutoPasT &autopas, FgetIter getIter, const std::vector<size_t> &particleIDsExpected) {
  std::vector<size_t> particleIDsFound;

#ifdef AUTOPAS_OPENMP
  // aparently the version from WrapOpenMP.h can not be found
#pragma omp declare reduction( \
        vecMergeWorkaround : std::vector<size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel reduction(vecMergeWorkaround : particleIDsFound)
#endif
  {
    for (auto iterator = getIter(); iterator.isValid(); ++iterator) {
      auto id = iterator->getID();
      particleIDsFound.push_back(id);
    }
  }

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(particleIDsExpected));
}

/**
 * Generates a given amount of cells where only indicated cells contain a given amount of particles.
 * Cells can be considered to be on the main diagonal through 3D space. So the xyz coordinates of each cell's lower
 * corner are {cellID, cellID, cellID}. The particles are also placed along this line within their cells. Within each
 * cell the particles are placed equidistant around the center.
 * @param numCells
 * @param cellsToFill
 * @param particlesPerCell
 * @return Vector of generated and filled cells.
 */
static std::vector<FMCell> generateCellsWithPattern(const size_t numCells, const std::vector<size_t> &cellsToFill,
                                                    const size_t particlesPerCell) {
  constexpr double cellDiagonal = 1.;
  // distance between particles within one cell
  const double distBetweenParticles = cellDiagonal / (particlesPerCell + 1.);

  std::vector<FMCell> cells(numCells);
  size_t numParticlesAdded = 0;
  for (auto cellId : cellsToFill) {
    for (size_t i = 0; i < particlesPerCell; ++i) {
      auto position = cellId + distBetweenParticles * (i + 1.);
      Molecule m({position, position, position}, {0, 0, 0}, numParticlesAdded++, 0);
      cells[cellId].addParticle(m);
    }
  }
  return cells;
}
}  // namespace IteratorTestHelper