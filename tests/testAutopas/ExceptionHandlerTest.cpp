/**
 * @file ExceptionHandlerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>

#include "ExceptionHandlerTest.h"
#include "utils/ExceptionHandler.h"

using autopas::utils::ExceptionHandler;

void ExceptionHandlerTest::SetUp() { autopas::logger::create(); }
void ExceptionHandlerTest::TearDown() { autopas::logger::unregister(); }

TEST_F(ExceptionHandlerTest, TestIgnore) {
  ExceptionHandler exceptionHandler(
      ExceptionHandler::ExceptionBehavior::ignore);
  EXPECT_NO_THROW(exceptionHandler.exception("testignore"));
}

TEST_F(ExceptionHandlerTest, TestThrow) {
  ExceptionHandler exceptionHandler(
      ExceptionHandler::ExceptionBehavior::throwException);
  EXPECT_ANY_THROW(exceptionHandler.exception("testignore"));

  EXPECT_ANY_THROW(exceptionHandler.exception(std::exception()));
}

TEST_F(ExceptionHandlerTest, TestAbort) {
  //  autopas::logger =
  //      std::unique_ptr<autopas::log::Logger>(new autopas::log::Logger());
  ExceptionHandler exceptionHandler(
      ExceptionHandler::ExceptionBehavior::printAbort);

  EXPECT_DEATH(exceptionHandler.exception("testignore"), "");

  EXPECT_DEATH(exceptionHandler.exception(std::exception()), "");
}

TEST_F(ExceptionHandlerTest, TestAbortCustom) {
  //  autopas::logger =
  //      std::unique_ptr<autopas::log::Logger>(new autopas::log::Logger());

  auto abortFunction = []() -> void {
    AutoPasLogger->error("TESTABORTCUSTOMCALL123");
    abort();
  };
  ExceptionHandler exceptionHandler(
      ExceptionHandler::ExceptionBehavior::printCustomAbortFunction);

  exceptionHandler.setCustomAbortFunction(abortFunction);

  EXPECT_DEATH(exceptionHandler.exception("testignore"),
               "");

  EXPECT_DEATH(exceptionHandler.exception(std::exception()),
               "");
}