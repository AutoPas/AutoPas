
#include "tests/zonalMethods/MidpointTest.h"

#include "molecularDynamicsLibrary/LJFunctorAVX.h"
#include "src/configuration/MDFlexConfig.h"

void MidpointTest::initContainer(AutoPasType &autopas, std::vector<ParticleType> particles) {
  autopas.setBoxMin(_boxMin);
  autopas.setBoxMax(_boxMax);
  autopas.setCutoff(_cutoff);
  autopas.init();

  // insert particles
  for (auto &particle : particles) {
    autopas.addParticle(particle);
  }
}

std::shared_ptr<ParticlePropertiesLibraryType> MidpointTest::initializeParticlePropertiesLibrary() {
  std::shared_ptr<ParticlePropertiesLibraryType> _particlePropertiesLibrary =
      std::make_shared<ParticlePropertiesLibraryType>(_cutoff);

  std::map<size_t, double> massMap{{0, 1}};
  std::map<size_t, double> epsilonMap{{0, 1}}, sigmaMap{{0, 1}};
  std::map<size_t, double> nuMap{{0, 0.0073}};

  // initialize site Ids with mandatory mass parameter
  for (auto [siteTypeId, mass] : massMap) {
    _particlePropertiesLibrary->addSiteType(siteTypeId, massMap.at(siteTypeId));
  }
  // check size of LJ site parameter vectors match
  if (epsilonMap.size() != sigmaMap.size()) {
    throw std::runtime_error(
        "MDFlexConfig::initializeParticlePropertiesLibrary(): Number of LJ site-level properties differ! Potentially "
        "missing epsilon or sigma for some LJ sites.");
  }
  // initialize LJ parameters
  for (auto [siteTypeId, epsilon] : epsilonMap) {
    _particlePropertiesLibrary->addLJParametersToSite(siteTypeId, epsilon, sigmaMap.at(siteTypeId));
  }
  // initialize AT parameters
  for (auto [siteTypeId, nu] : nuMap) {
    _particlePropertiesLibrary->addATParametersToSite(siteTypeId, nu);
  }

  _particlePropertiesLibrary->calculateMixingCoefficients();

  return _particlePropertiesLibrary;
}

TEST_F(MidpointTest, testMidpointInitialization) {
  ASSERT_EQ(_exportRegions.size(), 26);
  ASSERT_EQ(_importRegions.size(), 26);
  ASSERT_EQ(_interactionZones.size(), 26);

  auto sum = 0;
  for (auto schedule : _interactionSchedule) {
    sum += schedule.second.size();
  }

  /**
   * This is a weak check, but it suffices for now.
   * (see Spahl 2016)
   * Number of total interactions =
   *  13 interactions with the neighbour on the other side +
   *  6 times 8 interaction (see Figure 3.5) +
   *  12 times 2 interaction (see Figure 3.6) =
   *  13 + 48 + 24 = 85
   */
  ASSERT_EQ(sum, 85);

  // some individual checks
  auto schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({-1, 0, 0})));
  ASSERT_EQ(schedule.size(), 9);
  std::vector<std::string> expected;
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      expected.push_back(std::to_string(convRelNeighboursToIndex({1, i, j})));
    }
  }
  std::sort(schedule.begin(), schedule.end());
  std::sort(expected.begin(), expected.end());
  ASSERT_EQ(schedule, expected);

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({1, 1, 0})));
  ASSERT_EQ(schedule.size(), 2);
  expected.clear();
  expected.push_back(std::to_string(convRelNeighboursToIndex({-1, -1, 1})));
  expected.push_back(std::to_string(convRelNeighboursToIndex({-1, -1, -1})));

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({-1, -1, -1})));
  ASSERT_EQ(schedule.size(), 1);

  schedule = _interactionSchedule.at(std::to_string(convRelNeighboursToIndex({1, 1, 1})));
  ASSERT_EQ(schedule.size(), 0);
}

TEST_F(MidpointTest, testMidpointInteractionSchedule) {
  std::vector<ParticleType> particles;
  std::map<size_t, std::array<double, 3>> expectedForce;
  size_t id = 0;
  // iterate over neighbours
  for (int x = -1; x < 2; ++x) {
    for (int y = -1; y < 2; ++y) {
      for (int z = -1; z < 2; ++z) {
        if (x == 0 && y == 0 && z == 0) {
          continue;
        }
        int d[3] = {x, y, z};
        // if neighbour is in the stencil
        // particle inside the domain
        double px = (x == -1) * (_boxMin[0] - 0.5) + (x == 0) * (_boxMin[0] + 5) + (x == 1) * (_boxMax[0] + 0.5);
        double py = (y == -1) * (_boxMin[0] - 0.5) + (y == 0) * (_boxMin[0] + 5) + (y == 1) * (_boxMax[0] + 0.5);
        double pz = (z == -1) * (_boxMin[0] - 0.5) + (z == 0) * (_boxMin[0] + 5) + (z == 1) * (_boxMax[0] + 0.5);

        ParticleType p({px, py, pz}, {0, 0, 0}, id, 0);
        p.setF({0, 0, 0});
        particles.push_back(p);

        ++id;
      }
    }
  }

  using LJFunctorTypeAVX = mdLib::LJFunctorAVX<ParticleType, true, true, autopas::FunctorN3Modes::Both,
                                               mdFlexibleTypeDefs::calcGlobals, mdFlexibleTypeDefs::countFLOPs>;
  auto particlePropertiesLibrary = initializeParticlePropertiesLibrary();
  LJFunctorTypeAVX ljFunctor(_cutoff, *particlePropertiesLibrary);
  // calculate expected forces
  for (size_t i = 0; i < particles.size(); ++i) {
    for (size_t j = i + 1; j < particles.size(); ++j) {
      using namespace autopas::utils::ArrayMath::literals;
      auto midpoint = (particles[i].getR() + particles[j].getR()) * 0.5;
      auto dist = midpoint - _homeBoxRegion._origin;
      bool calc = true;
      for (size_t i = 0; i < 3; ++i) {
        if (dist.at(i) < 0) {
          calc = false;
          break;
        }
        if (dist.at(i) >= _homeBoxRegion._size.at(i)) {
          calc = false;
          break;
        }
      }
      if (calc) {
        ljFunctor.AoSFunctor(particles[i], particles[j], true);
        if (expectedForce.find(particles[i].getID()) == expectedForce.end()) {
          expectedForce.insert_or_assign(particles[i].getID(), particles[i].getF());
        } else {
          expectedForce.at(particles[i].getID()) += particles[i].getF();
        }
        if (expectedForce.find(particles[j].getID()) == expectedForce.end()) {
          expectedForce.insert_or_assign(particles[j].getID(), particles[j].getF());
        } else {
          expectedForce.at(particles[j].getID()) -= particles[j].getF();
        }
      }
    }
  }

  // reset all forces to zero
  for (auto &p : particles) {
    p.setF({0, 0, 0});
  }

  // init container
  initContainer(_autopas, {});
  _autopas.addHaloParticles(particles);

  calculateExternalZonalInteractions(_autopas, particlePropertiesLibrary, _cutoff);

  for (auto iter = _autopas.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getF(), expectedForce.at(iter->getID()));
  }
}
