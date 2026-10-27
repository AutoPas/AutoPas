/**
 * @file LogicHandlerTest.h
 * @author Manish
 * @date 13.05.24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/LogicHandler.h"
#include "autopas/tuning/TuningManager.h"
#include "testingHelpers/ArbitraryConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

class LogicHandlerTest : public AutoPasTestBase {
 public:
  std::unique_ptr<autopas::LogicHandler<Molecule>> _logicHandler;
  std::shared_ptr<autopas::TuningManager> _tuningManager;
  /**
   * Initializes _logicHandler and _tuningManager with a pairwise AutoTuner.
   * @param searchSpace The search space of the AutoTuner.
   */
  void initLogicHandler(const std::set<autopas::Configuration> &searchSpace = {
                            arbitraryConfigurations::_arbitrary_config_2B_0});
};
