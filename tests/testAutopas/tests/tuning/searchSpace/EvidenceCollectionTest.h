/**
 * @file EvidenceCollectionTest.h
 * @author F. Gratl
 * @date 05.07.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class EvidenceCollectionTest : public AutoPasTestBase {
 private:
  static constexpr int _threadCount = autopas::Configuration::ThreadCountNoTuning;

 public:
  static constexpr autopas::Configuration _configurationLC_C01 =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise, _threadCount);
  static constexpr autopas::Configuration _configurationLC_C08 =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise, _threadCount);

  static constexpr autopas::Configuration _configurationLC_Sliced =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise, _threadCount);
};
