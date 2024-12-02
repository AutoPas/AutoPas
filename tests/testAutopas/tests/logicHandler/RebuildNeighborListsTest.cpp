/**
* @file RebuildNeighborListsTest.cpp
* @author muehlhaeusser
* @date 29.11.2024
*/

#include "RebuildNeighborListsTest.h"

#include <unordered_map>
#include "autopas/LogicHandler.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/options/InteractionTypeOption.h"
#include "../../testingHelpers/commonTypedefs.h"

using ::testing::_;

std::set<autopas::Configuration> RebuildNeighborListsTest::getPairwiseConfigs() {
  const double _cellSizeFactor{1.};
  const autopas::Configuration _confVl_list_it_noN3{
      autopas::ContainerOption::verletLists,     _cellSizeFactor,
      autopas::TraversalOption::vl_list_iteration, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,     _cellSizeFactor,
      autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  return {_confLc_c01_noN3, _confVl_list_it_noN3};
}

std::set<autopas::Configuration> RebuildNeighborListsTest::getTriwiseConfigs() {
  const double _cellSizeFactor{1.};
  const autopas::Configuration _confVl_list_it_noN3_3B{
      autopas::ContainerOption::verletLists,     _cellSizeFactor,
      autopas::TraversalOption::vl_list_iteration, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};
  const autopas::Configuration _confLc_c01_noN3_3B{
      autopas::ContainerOption::linkedCells,     _cellSizeFactor,
      autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};
  const autopas::Configuration _confVl_pairlist_noN3_3B{
      autopas::ContainerOption::verletLists,     _cellSizeFactor,
      autopas::TraversalOption::vl_pair_list_iteration_3b, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};
  return {_confLc_c01_noN3_3B, _confVl_list_it_noN3_3B, _confVl_pairlist_noN3_3B};
}

TEST_P(RebuildNeighborListsTest, testRebuildDifferentContainerPairwiseTriwise) {
  ::testing::Sequence seq;

  // also check if rebuild is detected if next config is invalid

  const unsigned int verletRebuildFrequency = 3;
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 3.,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const auto [pairwiseConfig, triwiseConfig] = GetParam();

  const autopas::AutoTuner::SearchSpaceType searchSpace{
    std::get<0>(GetParam())
  };
  const autopas::AutoTuner::SearchSpaceType searchSpace3B{
    std::get<1>(GetParam())
  };

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  tunerMap.emplace(
      autopas::InteractionTypeOption::triwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace3B, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner2B = *tunerMap[autopas::InteractionTypeOption::pairwise];
  auto &autoTuner3B = *tunerMap[autopas::InteractionTypeOption::triwise];

  testing::NiceMock<MockPairwiseFunctor<Molecule>> pairwiseFunctor;
  testing::NiceMock<MockTriwiseFunctor<Molecule>> triwiseFunctor;

  auto p1 = Molecule{{5., 5., 5.},{0., 0., 0.}, 0};
  logicHandler.addParticle(p1);
  auto p2 = Molecule{{4., 5., 5.},{0., 0., 0.}, 0};
  logicHandler.addParticle(p2);
  auto p3 = Molecule{{5., 4., 5.},{0., 0., 0.}, 0};
  logicHandler.addParticle(p3);

  auto emptyVec = logicHandler.updateContainer();

  EXPECT_CALL(pairwiseFunctor, AoSFunctor(_, _, _)).Times(6);
  EXPECT_CALL(triwiseFunctor, AoSFunctor(_, _, _, _)).Times(3);

  logicHandler.computeInteractionsPipeline(&pairwiseFunctor, autopas::InteractionTypeOption::pairwise);
  logicHandler.computeInteractionsPipeline(&triwiseFunctor, autopas::InteractionTypeOption::triwise);
  ::testing::Mock::VerifyAndClearExpectations(&pairwiseFunctor);
  ::testing::Mock::VerifyAndClearExpectations(&triwiseFunctor);

  emptyVec = logicHandler.updateContainer();

  EXPECT_CALL(pairwiseFunctor, AoSFunctor(_, _, _)).Times(6);
  EXPECT_CALL(triwiseFunctor, AoSFunctor(_, _, _, _)).Times(3);

  logicHandler.computeInteractionsPipeline(&pairwiseFunctor, autopas::InteractionTypeOption::pairwise);
  logicHandler.computeInteractionsPipeline(&triwiseFunctor, autopas::InteractionTypeOption::triwise);
  ::testing::Mock::VerifyAndClearExpectations(&pairwiseFunctor);
  ::testing::Mock::VerifyAndClearExpectations(&triwiseFunctor);

  emptyVec = logicHandler.updateContainer();

  EXPECT_CALL(pairwiseFunctor, AoSFunctor(_, _, _)).Times(6);
  EXPECT_CALL(triwiseFunctor, AoSFunctor(_, _, _, _)).Times(3);

  logicHandler.computeInteractionsPipeline(&pairwiseFunctor, autopas::InteractionTypeOption::pairwise);
  logicHandler.computeInteractionsPipeline(&triwiseFunctor, autopas::InteractionTypeOption::triwise);
  ::testing::Mock::VerifyAndClearExpectations(&pairwiseFunctor);
  ::testing::Mock::VerifyAndClearExpectations(&triwiseFunctor);

  auto p4 = Molecule{{5., 5., 4.},{0., 0., 0.}, 0};
  logicHandler.addParticle(p4);

  emptyVec = logicHandler.updateContainer();
  EXPECT_CALL(pairwiseFunctor, AoSFunctor(_, _, _)).Times(12);
  EXPECT_CALL(triwiseFunctor, AoSFunctor(_, _, _, _)).Times(12);

  logicHandler.computeInteractionsPipeline(&pairwiseFunctor, autopas::InteractionTypeOption::pairwise);
  logicHandler.computeInteractionsPipeline(&triwiseFunctor, autopas::InteractionTypeOption::triwise);
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(Generated, RebuildNeighborListsTest,
                         Combine(ValuesIn(RebuildNeighborListsTest::getPairwiseConfigs()),
                                 ValuesIn(RebuildNeighborListsTest::getTriwiseConfigs())));