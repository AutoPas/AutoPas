/**
 * @file LogicHandlerTest.h
 * @author Manish
 * @date 13.05.24
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/LogicHandler.h"
#include "autopas/options/IteratorBehavior.h"
#include "testingHelpers/commonTypedefs.h"

class LogicHandlerTest : public AutoPasTestBase {
 public:
  std::shared_ptr<autopas::LogicHandler<Molecule>> _logicHandler;
  void initLogicHandler();
};
