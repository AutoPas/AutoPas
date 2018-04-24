/**
 * @file ExceptionHandlerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>

#include "ExceptionHandlerTest.h"
#include "utils/ExceptionHandler.h"

using autopas::utils::ExceptionHandler;
using autopas::utils::ExceptionBehavior;

void ExceptionHandlerTest::SetUp() { autopas::logger::create(); }
void ExceptionHandlerTest::TearDown() { autopas::logger::unregister(); }

TEST_F(ExceptionHandlerTest, TestIgnore) {
  ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
  EXPECT_NO_THROW(ExceptionHandler::exception("testignore"));
}

TEST_F(ExceptionHandlerTest, TestThrow) {
  ExceptionHandler::setBehavior(ExceptionBehavior::throwException);
  EXPECT_ANY_THROW(ExceptionHandler::exception("testignore"));

  EXPECT_ANY_THROW(ExceptionHandler::exception(std::exception()));
}

TEST_F(ExceptionHandlerTest, TestAbort) {
  //  autopas::logger =
  //      std::unique_ptr<autopas::log::Logger>(new autopas::log::Logger());
  ExceptionHandler::setBehavior(ExceptionBehavior::printAbort);

  EXPECT_DEATH(ExceptionHandler::exception("testignore"), "");

  EXPECT_DEATH(ExceptionHandler::exception(std::exception()), "");
}

TEST_F(ExceptionHandlerTest, TestAbortCustom) {
  //  autopas::logger =
  //      std::unique_ptr<autopas::log::Logger>(new autopas::log::Logger());

  auto abortFunction = []() -> void {
    AutoPasLogger->error("TESTABORTCUSTOMCALL123");
    abort();
  };
  ExceptionHandler::setBehavior(ExceptionBehavior::printCustomAbortFunction);

  ExceptionHandler::setCustomAbortFunction(abortFunction);

  EXPECT_DEATH(ExceptionHandler::exception("testignore"), "");

  EXPECT_DEATH(ExceptionHandler::exception(std::exception()), "");
}