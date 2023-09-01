/**
 * @file CellFunctorTest.cpp
 * @author D. Martin
 * @date 29.08.23
 */

#include "CellFunctorTest.h"

typedef ::testing::Types<
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, false, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, false, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, false, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, false, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, true, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, true, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::aos, true, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::aos, true, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, false, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, false, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, false, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, false, true>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, true, false>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, true, false>>,
    std::pair<autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>, mdLib::LJFunctor<Molecule>,
                                             autopas::DataLayoutOption::soa, true, true>,
              CellFunctorWrapper<autopas::DataLayoutOption::soa, true, true>>>
    CellFTestingTypes;

template <class T>
std::tuple<double, double> ownedHaloInteractionHelper(T &cellFunctor, const autopas::OwnershipState os1,
                                                      const autopas::OwnershipState os2,
                                                      const autopas::DataLayoutOption dataLayout,
                                                      mdLib::LJFunctor<Molecule> &ljFunctor) {
  autopas::FullParticleCell<Molecule> cell1({1., 1., 1.});
  autopas::FullParticleCell<Molecule> cell2({1., 1., 1.});

  // insert one owned particle in cell1
  Molecule p1({0.6, 0.5, 0.5}, {0., 0., 0.}, 0);
  p1.setOwnershipState(os1);
  cell1.addParticle(p1);

  // insert one owned particle in cell2
  Molecule p2({1.4, 0.5, 0.5}, {0., 0., 0.}, 1);
  p2.setOwnershipState(os2);
  cell2.addParticle(p2);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ljFunctor.SoALoader(cell1, cell1._particleSoABuffer, 0);
    ljFunctor.SoALoader(cell2, cell2._particleSoABuffer, 0);
  }

  cellFunctor.processCellPair(cell1, cell2);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    ljFunctor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    ljFunctor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
  }

  return std::tuple<double, double>{std::abs(cell1[0].getF()[0]), std::abs(cell2[0].getF()[0])};
}

TYPED_TEST_CASE_P(CellFunctorTest);

TYPED_TEST_P(CellFunctorTest, testOwnedAndHaloCellInteraction) {
  const double cutoff = 1.;
  const double sigma = 1.;
  const double epsilon = 1.;

  using CellFT = typename TypeParam::first_type;
  using CFWT = typename TypeParam::second_type;

  for (auto ownerShipStateA : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
    for (auto ownerShipStateB : {autopas::OwnershipState::owned, autopas::OwnershipState::halo}) {
      mdLib::LJFunctor<Molecule> ljFunctor(cutoff);
      ljFunctor.setParticleProperties(sigma, epsilon);

      ljFunctor.initTraversal();

      CellFT cellFunctor(&ljFunctor, cutoff);
      CFWT cfw;

      std::tuple<double, double> result =
          ownedHaloInteractionHelper<CellFT>(cellFunctor, ownerShipStateA, ownerShipStateB, cfw.dataLayoutV, ljFunctor);

      ljFunctor.endTraversal(cfw.useNewton3V);

      if (ownerShipStateA == autopas::OwnershipState::halo and ownerShipStateB == autopas::OwnershipState::halo) {
        EXPECT_NEAR(std::get<0>(result), 0, 1e-13)
            << "Particle in cell 0 is experiencing force with both cells containing only halo particles";
        EXPECT_NEAR(std::get<1>(result), 0, 1e-13)
            << "Particle in cell 1 is experiencing force with with both cells containing only halo particles";
      } else if (ownerShipStateA == autopas::OwnershipState::owned) {
        // if cell1 contains owned particles forces need to be calculated
        EXPECT_TRUE(std::get<0>(result) > 0)
            << "Particle in cell 0 does not experience force with OwnershipState owned";
      } else if (ownerShipStateA == autopas::OwnershipState::halo and (not cfw.useNewton3V) and
                 (not cfw.bidirectionalV)) {
        // if cell1 is halo, and NoN3 and no bidirectional we can skip the interaction
        EXPECT_NEAR(std::get<0>(result), 0, 1e-13)
            << "Particle in cell 0 is experiencing force with OwnershipState halo, no newton3 and bidirectional off";
      } else {
        // in all other cases we expect forces on both particles
        EXPECT_TRUE(std::get<0>(result) > 0) << "Particle in cell 0 does not experience force";
        EXPECT_TRUE(std::get<1>(result) > 0) << "Particle in cell 1 does not experience force";
      }
    }
  }
}

REGISTER_TYPED_TEST_CASE_P(CellFunctorTest, testOwnedAndHaloCellInteraction);
INSTANTIATE_TYPED_TEST_CASE_P(TypedTest, CellFunctorTest, CellFTestingTypes);