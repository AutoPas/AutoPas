/**
 * @file LinkedCellsTest.cpp
 * @author seckler
 * @date 27.04.18
 */

#include "LinkedCellsTest.h"

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayUtils.h"

TYPED_TEST_SUITE_P(LinkedCellsTest);

TYPED_TEST_P(LinkedCellsTest, testUpdateContainer) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> zero{0., 0., 0.};
  const std::array<double, 3> boxMin = zero;
  const std::array<double, 3> boxMax{4.5, 4.5, 4.5};
  // set values so we have 3x3x3 cells + halo = 5x5x5
  const double cutoff{1.0};
  const double skin{0.5};
  typename TestFixture::LinkedCellsType linkedCells(boxMin, boxMax, cutoff, skin);

  // create owned particles
  const std::vector<ParticleFP64> ownedParticles{
      // clang-format off
      {{0.5, 0.5, 0.5}, zero, 0},
      {{1.5, 1.5, 1.5}, zero, 1},
      {{1.6, 1.5, 1.5}, zero, 2},
      {{4.4, 1.5, 1.5}, zero, 3},
      {{4.0, 4.0, 4.0}, zero, 4},
      // clang-format on
  };

  // These are going to be halo particles
  const std::vector<ParticleFP64> haloParticles{
      {{-0.5, +1.5, +1.5}, zero, 5, autopas::OwnershipState::halo},
      {{+5.0, +1.5, +1.5}, zero, 6, autopas::OwnershipState::halo},
      {{+1.5, -0.5, +1.5}, zero, 7, autopas::OwnershipState::halo},
      {{+1.5, +1.5, -0.5}, zero, 8, autopas::OwnershipState::halo},
  };

  // calculate the cell IDs for each particle
  std::unordered_map<size_t, size_t> particleIdToCellIdMap;
  for (const auto &particleCollection : {ownedParticles, haloParticles}) {
    for (const auto &particle : particleCollection) {
      const auto cellID = linkedCells.getCellBlock().get1DIndexOfPosition(particle.getR());
      particleIdToCellIdMap[particle.getID()] = cellID;
    }
  }
  ASSERT_EQ(particleIdToCellIdMap.size(), ownedParticles.size() + haloParticles.size())
      << "There should be exactly one cellID per particleID.\n"
         "Either the test is set up wrong or get1DIndexOfPosition is broken.";

  // we insert owned and halo particles alternating. This way we can check if references are updated correctly when
  // using LinkedCellsReferences
  linkedCells.addParticle(ownedParticles[0]);
  linkedCells.addHaloParticle(haloParticles[0]);
  linkedCells.addParticle(ownedParticles[1]);
  linkedCells.addHaloParticle(haloParticles[1]);
  linkedCells.addParticle(ownedParticles[2]);
  linkedCells.addHaloParticle(haloParticles[2]);
  linkedCells.addParticle(ownedParticles[3]);
  linkedCells.addHaloParticle(haloParticles[3]);
  linkedCells.addParticle(ownedParticles[4]);

  this->checkParticleIDsInCells(
      linkedCells,
      {{particleIdToCellIdMap[8], {{8, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[0], {{0, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[7], {{7, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[5], {{5, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[3], {{3, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[1], {{1, autopas::OwnershipState::owned}, {2, autopas::OwnershipState::owned}}},
       {particleIdToCellIdMap[6], {{6, autopas::OwnershipState::halo}}},
       {particleIdToCellIdMap[4], {{4, autopas::OwnershipState::owned}}}},
      true, __LINE__);

  // // new locations for owned particles
  linkedCells.getCells()[particleIdToCellIdMap[0]].begin()->addR({+2.0, +0.0, +0.0});  // move to {1.5, 0.5, 0.5}
  linkedCells.getCells()[particleIdToCellIdMap[1]].begin()->addR({-1.0, -0.0, -0.0});  // move to {-0.5, 1.5, 1.5}
  linkedCells.getCells()[particleIdToCellIdMap[3]].begin()->addR({+0.2, +0.0, -1.0});  // move to {5.0, 1.5, 0.5}
  linkedCells.getCells()[particleIdToCellIdMap[4]].begin()->addR({-0.9, -2.0, -2.0});  // move to {1.6, 0.5, 0.5}

  std::vector<ParticleFP64> invalidParticles;
  EXPECT_NO_THROW(invalidParticles = linkedCells.updateContainer(this->_keepListsValid));
  EXPECT_EQ(invalidParticles.size(), 1);
  EXPECT_EQ(invalidParticles[0].getID(), 3);

  if (this->_keepListsValid) {
    // if the lists are kept valid, particles are NOT moved between cells!
    // halo particles should now be dummies
    // particle 3 should be a leaving particle and therefore a dummy
    this->checkParticleIDsInCells(
        linkedCells,
        {{particleIdToCellIdMap[8], {{8, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[0], {{0, autopas::OwnershipState::owned}}},
         {particleIdToCellIdMap[7], {{7, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[5], {{5, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[3], {{3, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[1], {{1, autopas::OwnershipState::owned}, {2, autopas::OwnershipState::owned}}},
         {particleIdToCellIdMap[6], {{6, autopas::OwnershipState::dummy}}},
         {particleIdToCellIdMap[4], {{4, autopas::OwnershipState::owned}}}},
        true, __LINE__);
  } else {
    // if the lists are not kept valid, particles should be moved between cells, so update the cells!
    // halo particles should be removed by updateContainer() at this point
    this->checkParticleIDsInCells(
        linkedCells,
        {
            {particleIdToCellIdMap[0] + 1, {{0, autopas::OwnershipState::owned}}},  // moved one cell to the right
            {particleIdToCellIdMap[1] - 1, {{1, autopas::OwnershipState::owned}}},  // moved one cell to the left
            {particleIdToCellIdMap[2], {{2, autopas::OwnershipState::owned}}},      // didn't change cell
            {particleIdToCellIdMap[4] - (0 + 1 * 5 + 1 * 5 * 5),
             {{4, autopas::OwnershipState::owned}}},  // moved one cell to the front and one down
        },
        false /*here, we do not know the order!*/, __LINE__);
  }
}

TYPED_TEST_P(LinkedCellsTest, testUpdateContainerCloseToBoundary) {
  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax{10., 10., 10.};
  const double cutoff{1.5};
  const double skin{1.};  // particles are moved by up to 0.5 and lists might be kept valid
  typename TestFixture::LinkedCellsType linkedCells(boxMin, boxMax, cutoff, skin);

  int id = 1;
  for (const double x : {0., 5., 9.999}) {
    for (const double y : {0., 5., 9.999}) {
      for (const double z : {0., 5., 9.999}) {
        const ParticleFP64 p({x, y, z}, {0., 0., 0.}, id++);
        EXPECT_NO_THROW(linkedCells.addParticle(p));  // inside, therefore ok!
      }
    }
  }
  std::set<unsigned long> movedIDs;
  // we move particles that are close to the boundary to outside the container and remember their IDs
  for (auto iter = linkedCells.begin(); iter.isValid(); ++iter) {
    for (unsigned short dim = 0; dim < 3; ++dim) {
      if (iter->getR()[dim] < 0.5) {
        auto r = iter->getR();
        // smallest double smaller than 0
        r[dim] = std::nexttoward(0., -1.);
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
      if (iter->getR()[dim] > 9.5) {
        auto r = iter->getR();
        r[dim] = 10.;
        iter->setR(r);
        movedIDs.insert(iter->getID());
      }
    }
  }

  // now update the container!
  const auto invalidParticles = linkedCells.updateContainer(this->_keepListsValid);
  // the particles should no longer be in the inner cells!
  for (auto iter = linkedCells.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    EXPECT_EQ(movedIDs.count(iter->getID()), 0)
        << "Particle " << iter->getID() << " at " << autopas::utils::ArrayUtils::to_string(iter->getR())
        << " is still in an inner cell although it was moved!";
  }

  // the particles should now be inside the invalidParticles vector!
  EXPECT_EQ(movedIDs.size(), invalidParticles.size());
  for (const auto &particle : invalidParticles) {
    EXPECT_EQ(movedIDs.count(particle.getID()), 1)
        << "Particle " << particle.getID() << " at " << autopas::utils::ArrayUtils::to_string(particle.getR())
        << " was not returned by updateContainer()!";
  }
}

REGISTER_TYPED_TEST_SUITE_P(LinkedCellsTest, testUpdateContainer, testUpdateContainerCloseToBoundary);

// Workaround for storing two types.
// Currently, clang produces bugs if one tries to store this using a tuple or a pair.
// This problem is described in P0641R2 and occurs if an explicitly defaulted constructor cannot be instantiated.
// This is fixed in c++20.
template <typename first, typename second>
struct two_values {
  using first_t = first;
  using second_t = second;
};

// defines the types of linkedCells and _keepListsValid
struct LC_KeepListsValid : two_values<autopas::LinkedCells<ParticleFP64>, std::true_type> {};
struct LC_DontKeepListsValid : two_values<autopas::LinkedCells<ParticleFP64>, std::false_type> {};
struct LCRef_KeepListsValid : two_values<autopas::LinkedCellsReferences<ParticleFP64>, std::true_type> {};
struct LCRef_DontKeepListsValid : two_values<autopas::LinkedCellsReferences<ParticleFP64>, std::false_type> {};

using MyTypes =
    ::testing::Types<LC_KeepListsValid, LC_DontKeepListsValid, LCRef_KeepListsValid, LCRef_DontKeepListsValid>;

/// @todo c++20: replace with:
// using MyTypes = ::testing::Types<std::tuple<autopas::LinkedCells<ParticleFP64>, std::true_type>,
//                                  std::tuple<autopas::LinkedCells<ParticleFP64>, std::false_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<ParticleFP64>, std::true_type>,
//                                  std::tuple<autopas::LinkedCellsReferences<ParticleFP64>, std::false_type> >;

INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, LinkedCellsTest, MyTypes);
