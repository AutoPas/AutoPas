/**
 * @file LogicHandlerTest.h
 * @author Manish
 * @date 13.05.24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/LogicHandler.h"
#include "autopas/tuning/TunerManager.h"
#include "testingHelpers/commonTypedefs.h"

class LogicHandlerTest : public AutoPasTestBase {
 public:
  std::unique_ptr<autopas::LogicHandler<Molecule>> _logicHandler;
  std::shared_ptr<autopas::TunerManager> _tunerManager;
  void initLogicHandler();
};
