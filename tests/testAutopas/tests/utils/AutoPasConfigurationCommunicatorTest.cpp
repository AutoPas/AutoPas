/**
 * @file AutoPasConfigurationCommunicatorTest.cpp
 * @author muehlhaeusser
 * @date 16.05.24
 */
#include "AutoPasConfigurationCommunicatorTest.h"

#include <gmock/gmock-matchers.h>

#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"

/**
 * Takes a set of configurations, serializes them and deserializes them back.
 * Then checks if the deserialized configs are the same as before.
 * @param configs
 */
void AutoPasConfigurationCommunicatorTest::testConfigsCommunication(std::set<autopas::Configuration> &configs) {
  auto configsVector = std::vector<autopas::Configuration>(configs.begin(), configs.end());

  // Serialize
  auto byteVector = autopas::utils::AutoPasConfigurationCommunicator::serializeConfigurations(configsVector);
  // Deserialize
  auto configsDeserialized = autopas::utils::AutoPasConfigurationCommunicator::deserializeConfigurations(byteVector);

  EXPECT_THAT(configsVector, ::testing::ContainerEq(configsDeserialized));
}

/**
 * Creates a searchSpace (set of configurations) and tests if the configurations are the same after serializing and
 * deserializing them.
 */
TEST_F(AutoPasConfigurationCommunicatorTest, SerializationTest) {
  // Test pairwise configurations
  auto pairwiseSearchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, pairwiseTraversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options,
      &cellSizeFactors, autopas::InteractionTypeOption::pairwise);
  auto maxSearchSpaceSize = autopas::utils::AutoPasConfigurationCommunicator::getSearchSpaceSize(
      containerOptions, cellSizeFactors, pairwiseTraversalOptions, loadEstimatorOptions, dataLayoutOptions,
      newton3Options, autopas::InteractionTypeOption::pairwise);

  EXPECT_GE(maxSearchSpaceSize, pairwiseSearchSpace.size());
  testConfigsCommunication(pairwiseSearchSpace);

  // Test triwise configurations
  auto triwiseSearchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, triwiseTraversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options,
      &cellSizeFactors, autopas::InteractionTypeOption::triwise);
  maxSearchSpaceSize = autopas::utils::AutoPasConfigurationCommunicator::getSearchSpaceSize(
      containerOptions, cellSizeFactors, triwiseTraversalOptions, loadEstimatorOptions, dataLayoutOptions,
      newton3Options, autopas::InteractionTypeOption::triwise);

  EXPECT_GE(maxSearchSpaceSize, triwiseSearchSpace.size());
  testConfigsCommunication(triwiseSearchSpace);
}