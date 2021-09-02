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

namespace ForEachTestHelper {

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
  auto skin = autoPas.getVerletSkin();

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
          autoPas.addOrUpdateHaloParticle(p);
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
 * Apply an iterator, track what particle IDs are found and compare this to a vector of expected IDs
 * @tparam AutoPasT
 * @tparam IteratorT
 * @param autopas
 * @param iterator
 * @param particleIDsExpected
 */
template <class AutoPasT, class Lambda>
void forEachParticleTest(AutoPasT &autopas, Lambda forEachInRegionLambda,
                         const std::vector<size_t> &particleIDsExpected) {
  std::vector<size_t> particleIDsFound;

  //#ifdef AUTOPAS_OPENMP
  //  // aparently the version from WrapOpenMP.h can not be found
  //#pragma omp declare reduction(vecMergeWorkaround : std::vector<size_t> : omp_out.insert(omp_out.end(),
  // omp_in.begin(), omp_in.end())) #pragma omp parallel reduction(vecMergeWorkaround : particleIDsFound) #endif
  {
    auto lambda = [&](auto &p) {
      auto id = p.getID();
      particleIDsFound.push_back(id);
    };
    forEachInRegionLambda(lambda);
  }

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(particleIDsExpected));
}

/**
 * Apply an iterator, track what particle IDs are found and compare this to a vector of expected IDs
 * @tparam AutoPasT
 * @tparam IteratorT
 * @param autopas
 * @param iterator
 * @param particleIDsExpected
 */
template <class AutoPasT, class Lambda>
void reduceParticlesTest(AutoPasT &autopas, Lambda reduceInRegionLambda,
                         const std::vector<size_t> &particleIDsExpected) {
  std::vector<size_t> particleIDsFound;
  size_t reductionValue = 0ul;

  //#ifdef AUTOPAS_OPENMP
  //  // aparently the version from WrapOpenMP.h can not be found
  //#pragma omp declare reduction(vecMergeWorkaround : std::vector<size_t> : omp_out.insert(omp_out.end(),
  // omp_in.begin(), omp_in.end())) #pragma omp parallel reduction(vecMergeWorkaround : particleIDsFound) #endif
  {
    auto lambda = [&](auto &p, size_t &rv) {
      auto id = p.getID();
      rv += id;
      particleIDsFound.push_back(id);
    };
    reduceInRegionLambda(lambda, reductionValue);
  }

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(particleIDsExpected));

  size_t expectedReductionValue = std::accumulate(particleIDsExpected.begin(), particleIDsExpected.end(), 0ul);
  EXPECT_EQ(reductionValue, expectedReductionValue);
}
}  // namespace ForEachTestHelper
