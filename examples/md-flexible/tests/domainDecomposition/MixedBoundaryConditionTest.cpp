/**
 * @file MixedBoundaryConditionTest.cpp
 * @author S. J. Newcome
 * @date 21/01/2022
 */
#include "MixedBoundaryConditionTest.h"

#include "autopas/AutoPasDecl.h"
#include "src/TypeDefinitions.h"
#include "src/domainDecomposition/DomainTools.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

extern template class autopas::AutoPas<ParticleType>;

auto MixedBoundaryConditionTest::setUpExpectations(
    const std::vector<std::array<double, 3>> &particlePositions, const std::array<double, 3> &boxMin,
    const std::array<double, 3> &boxMax, const double sigma, const double interactionLength,
    const std::array<options::BoundaryTypeOption, 3> &boundaryConditions) {
  const auto boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  auto expectedPositions = particlePositions;
  auto expectedHaloPositions = particlePositions;
  std::vector<std::array<double, 3>> expectedForces;
  expectedForces.resize(particlePositions.size());

  auto forceFromReflection = [&](const std::array<double, 3> position, const int dimensionOfBoundary,
                                 const bool isUpper) {
    const auto distanceToBoundary = isUpper ? boxMax[dimensionOfBoundary] - position[dimensionOfBoundary]
                                            : position[dimensionOfBoundary] - boxMin[dimensionOfBoundary];
    const auto distanceToMirrorParticle = distanceToBoundary * 2.;
    const auto distanceSquared = distanceToMirrorParticle * distanceToMirrorParticle;

    const auto inverseDistanceSquared = 1. / distanceSquared;
    const auto lj2 = sigma * sigma * inverseDistanceSquared;
    const auto lj6 = lj2 * lj2 * lj2;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto ljForceFactor = 24 * (lj12 + lj12m6) * inverseDistanceSquared;
    const auto force = ljForceFactor * distanceToMirrorParticle * (isUpper ? -1. : 1.);

    return force;
  };

  for (size_t id = 0; id < particlePositions.size(); ++id) {
    for (size_t dim = 0; dim < boundaryConditions.size(); ++dim) {
      switch (boundaryConditions[dim]) {
        case ::options::BoundaryTypeOption::periodic:
          // if a periodic boundary is crossed the position is shifted inside the domain
          if (particlePositions[id][dim] > boxMax[dim]) {
            expectedPositions[id][dim] -= boxLength[dim];
          } else if (particlePositions[id][dim] < boxMin[dim]) {
            expectedPositions[id][dim] += boxLength[dim];
          }
          // if particles are outside the domain, they will be shifted to the other side,
          // but a halo particle will spawn where they were
          if (particlePositions[id][dim] > (boxMax[dim]) or particlePositions[id][dim] < (boxMin[dim])) {
            expectedHaloPositions[id][dim] = particlePositions[id][dim];
          } else {
            // if near a periodic boundary create a halo particle
            if (particlePositions[id][dim] > (boxMax[dim] - interactionLength)) {
              expectedHaloPositions[id][dim] -= boxLength[dim];
            } else if (particlePositions[id][dim] < (boxMin[dim] + interactionLength)) {
              expectedHaloPositions[id][dim] += boxLength[dim];
            }
          }
          break;
        case ::options::BoundaryTypeOption::reflective:
          // if near a reflective boundary and flying towards it the velocity sign is flipped
          if (particlePositions[id][dim] < boxMin[dim] + sixthRootOfTwo * sigma / 2.) {
            expectedForces[id][dim] = forceFromReflection(particlePositions[id], (int)dim, false);
          } else if (particlePositions[id][dim] > boxMax[dim] - sixthRootOfTwo * sigma / 2.) {
            expectedForces[id][dim] = forceFromReflection(particlePositions[id], (int)dim, true);
          }
          break;
        case ::options::BoundaryTypeOption::none:
          break;
      }
      // if a expected position is outside the domain it is a halo position
      if (autopas::utils::notInBox(expectedPositions[id], boxMin, boxMax)) {
        expectedHaloPositions[id] = expectedPositions[id];
      }
    }
  }
  return std::make_tuple(expectedPositions, expectedHaloPositions, expectedForces);
}

void MixedBoundaryConditionTest::testFunction(const std::vector<std::array<double, 3>> &particlePositions,
                                              const std::array<options::BoundaryTypeOption, 3> &boundaryConditions) {
  const double cutoff = 0.3;
  const double sigma = 1.;

  // initialise AutoPas container & domainDecomposition
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();
  config.boxMin.value = {0., 0., 0.};
  config.boxMax.value = {5., 5., 5.};
  config.cutoff.value = cutoff;
  config.subdivideDimension.value = {true, true, true};
  config.boundaryOption.value = boundaryConditions;
  config.addParticleType(0, 1., sigma, 1.);

  const std::array<double, 3> boxLength = autopas::utils::ArrayMath::sub(config.boxMax.value, config.boxMin.value);
  RegularGridDecomposition domainDecomposition(config);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);
  auto particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff);

  autoPasContainer->setBoxMin(domainDecomposition.getLocalBoxMin());
  autoPasContainer->setBoxMax(domainDecomposition.getLocalBoxMax());
  autoPasContainer->setCutoff(config.cutoff.value);
  autoPasContainer->init();

  particlePropertiesLibrary->addType(0, 1., sigma, 1.);
  particlePropertiesLibrary->calculateMixingCoefficients();

  const auto &[expectedPositions, expectedHaloPositions, expectedForces] = setUpExpectations(
      particlePositions, config.boxMin.value, config.boxMax.value, sigma,
      config.cutoff.value + config.verletSkinRadiusPerTimestep.value * config.verletRebuildFrequencies.value->getMax(),
      config.boundaryOption.value);

  // particles need to be added at positions inside the domain
  // but also close to their designated positions, so they end up in the correct MPI rank
  std::vector<std::array<double, 3>> particlePositionsSafe;
  std::transform(particlePositions.cbegin(), particlePositions.cend(), std::back_inserter(particlePositionsSafe),
                 [&](const auto &pos) {
                   auto safePos = pos;
                   for (size_t i = 0; i < pos.size(); ++i) {
                     // if a position is outside the (global) domain shift it inside
                     if (pos[i] > domainDecomposition.getGlobalBoxMax()[i]) {
                       // next float inside the domain
                       safePos[i] = std::nextafter(domainDecomposition.getGlobalBoxMax()[i],
                                                   domainDecomposition.getGlobalBoxMax()[i] - 1);
                     } else if (pos[i] < domainDecomposition.getGlobalBoxMin()[i]) {
                       // next float inside the domain
                       safePos[i] = std::nextafter(domainDecomposition.getGlobalBoxMin()[i],
                                                   domainDecomposition.getGlobalBoxMin()[i] + 1);
                     }
                   }
                   return safePos;
                 });

  // add particles in safe positions, so they can be added to the correct MPI rank.
  for (size_t i = 0; i < particlePositions.size(); ++i) {
    if (not domainDecomposition.isInsideLocalDomain(particlePositionsSafe[i])) {
      continue;
    }
    ParticleType particle;
    particle.setID(i);
    particle.setR(particlePositionsSafe[i]);
    particle.setF({0., 0., 0.});
    autoPasContainer->addParticle(particle);
  }

#if not defined(AUTOPAS_INCLUDE_MPI)
  ASSERT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), particlePositions.size())
      << "Not all particles were added as owned.";
#endif

  // change particle positions so that some of them crossed boundaries
  for (auto particleIter = autoPasContainer->begin(autopas::IteratorBehavior::owned); particleIter.isValid();
       ++particleIter) {
    const auto id = particleIter->getID();
    particleIter->setR(particlePositions[id]);
  }

#if not defined(AUTOPAS_INCLUDE_MPI)
  EXPECT_EQ(expectedPositions.size(), autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned));
#endif

  auto emigrants = autoPasContainer->updateContainer();

  // Apply boundary conditions
  domainDecomposition.exchangeMigratingParticles(*autoPasContainer, emigrants);
  domainDecomposition.reflectParticlesAtBoundaries(*autoPasContainer, *particlePropertiesLibrary);
  domainDecomposition.exchangeHaloParticles(*autoPasContainer);

#if not defined(AUTOPAS_INCLUDE_MPI)
  // if there are no "none" boundaries the number of owned particles should stay constant
  if (std::none_of(boundaryConditions.cbegin(), boundaryConditions.cend(),
                   [](const auto &bc) { return bc == ::options::BoundaryTypeOption::none; })) {
    EXPECT_EQ(expectedPositions.size(), autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned));
  }
#endif

  // check owned and halo separately
  for (const auto &[iterBehavior, expectedPos, expectedForce] :
       std::vector<std::tuple<autopas::IteratorBehavior, decltype(expectedPositions), decltype(expectedForces)>>{
           {autopas::IteratorBehavior::owned, expectedPositions, expectedForces},
           {autopas::IteratorBehavior::halo, expectedHaloPositions, expectedForces}}) {
    for (auto particle = autoPasContainer->begin(iterBehavior); particle.isValid(); ++particle) {
      const auto id = particle->getID();
#if defined(AUTOPAS_INCLUDE_MPI)
      // In the MPI case only one rank can check each particle
      if (not domainDecomposition.isInsideLocalDomain(expectedPositions[id])) {
        continue;
      }
#endif
      const auto &position = particle->getR();
      const auto &force = particle->getF();
      for (size_t dim = 0; dim < position.size(); ++dim) {
        EXPECT_NEAR(position[dim], expectedPos[id][dim], 1e-13)
            << "Unexpected position in dim " << dim << " for " << iterBehavior.to_string() << particle->toString();
        EXPECT_NEAR(force[dim], expectedForce[id][dim], 1e-13)
            << "Unexpected force in dim " << dim << " for " << iterBehavior.to_string() << particle->toString();
      }
    }
  }
}

/**
 * Simple test designed to show that reflection at a corner between 2 reflective faces works correctly.
 * Places particles within reflective skin and tests for correct reflective behaviour including for:
 * - cases where the particle's velocity is away from the boundary
 * - cases where the particle is in a corner, including
 *   + reflective/reflective
 *   + reflective/periodic
 *   + reflective/reflective/periodic
 * - cases for both upper and lower boundaries
 */
TEST_F(MixedBoundaryConditionTest, testMixedReflection) {
  const std::array<options::BoundaryTypeOption, 3> boundaryConditions = {options::BoundaryTypeOption::periodic,
                                                                         options::BoundaryTypeOption::reflective,
                                                                         options::BoundaryTypeOption::reflective};
  const std::vector<std::array<double, 3>> particlePositions = {{1.0, 0.005, 1.0},     {1.0, 0.005, 0.005},
                                                                {4.995, 0.005, 0.005}, {4.0, 4.995, 4.0},
                                                                {4.0, 4.995, 4.995},   {0.005, 4.995, 4.995}};
  // particle 0 tests for reflection along a reflective boundary
  // particle 1 tests for no reflection along a reflective boundary when the direction is into the domain
  // particle 2 tests for reflection in the edge between two reflective boundaries
  // particle 3 tests for reflection in a single direction along a reflective/reflective edge when in one dimension the
  //    particle is travelling towards the boundary and the other away
  // particle 4 tests for reflection only in the directions with reflective boundaries in a periodic/refl/refl corner
  // particles 5-9 do the same as 0-4 respectively, except in the rightmost boundaries

  testFunction(particlePositions, boundaryConditions);
}

/**
 * Designed to test that exchangeMigratingParticles and exchangeHaloParticles in the mixed boundary case
 * Note: this is not designed to replace the more extensive tests in RegularGridDecompositionTest, but to test the
 * periodic BC in the mixed case.
 * Places particles in reflective skin and (primarily) outside of periodic boundary to test that particles
 * translated (as a result of the periodic boundary) and halo particles have the correct reflections.
 */
TEST_F(MixedBoundaryConditionTest, testPeriodic) {
  const std::array<options::BoundaryTypeOption, 3> boundaryConditions = {options::BoundaryTypeOption::periodic,
                                                                         options::BoundaryTypeOption::reflective,
                                                                         options::BoundaryTypeOption::reflective};

  const std::vector<std::array<double, 3>> particlePositions = {
      {-0.005, 0.005, 1.0}, {-0.005, 0.005, 2.0}, {-0.005, 0.005, 3.0}, {4.995, 4.995, 4.995}, {5.005, 4.995, 4.995}};
  const std::vector<std::array<double, 3>> particleVelocities = {
      {-1.0, -1.0, 0.0}, {-1.0, 1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

  // particle 0 tests that a particle that needs to be periodic translated in x, is also reflected correctly in y
  // particle 1 tests that a particle that needs to be periodic translated in x, whilst within the reflective skin of
  //    a y-boundary but moving away, is not reflected
  // particle 2 tests that a particle that needs to be periodic translated in x, but is moving towards the domain, is
  //    also correctly translated in x + reflected in y
  // particle 3 tests that a particle in a periodic/reflective/reflective corner produces a correctly reflected halo
  //    particle
  // particle 4 tests that a particle that needs is beyond the right x-boundary is also correctly reflected in y and z

  testFunction(particlePositions, boundaryConditions);
}

/**
 * Checks that the logic used to decide whether or not to implement a certain boundary condition is not only working in
 * the case where the boundaries are exclusively periodic or reflective. This should imply that there won't be any
 * problems with adding additional boundary condition types in the future.
 *
 * This test does this by confirming that no periodic or reflective behaviour occurs when no boundary conditions are
 * used. Two sets of particles are placed, one for the upper and lower boundaries, each of two particles - one just
 * inside the boundary, the other just outside - and nothing should change.
 */
TEST_F(MixedBoundaryConditionTest, testNoBoundary) {
  const std::array<options::BoundaryTypeOption, 3> boundaryConditions = {
      options::BoundaryTypeOption::none, options::BoundaryTypeOption::none, options::BoundaryTypeOption::none};

  const std::vector<std::array<double, 3>> particlePositions = {
      {-0.005, 2.5, 2.5}, {0.005, 2.5, 2.5}, {4.995, 2.5, 2.5}, {5.005, 2.5, 2.5}};
  const std::vector<std::array<double, 3>> particleVelocities = {
      {-1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};

  // particle 0 tests the lack of periodic translation in the left x-boundary
  // particle 1 tests the lack of reflection in the left x-boundary
  // particle 2 tests the lack of reflection in the right x-boundary
  // particle 3 tests the lack of periodic translation in the right x-boundary

  testFunction(particlePositions, boundaryConditions);
}