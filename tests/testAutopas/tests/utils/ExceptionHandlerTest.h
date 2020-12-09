/**
 * @file ExceptionHandlerTest.h
 * @author seckler
 * @date 20.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"

namespace ExceptionHandlerTest {

class ExceptionHandlerTest : public AutoPasTestBase {
  void SetUp() override;

  void TearDown() override;
};

}  // end namespace ExceptionHandlerTest
