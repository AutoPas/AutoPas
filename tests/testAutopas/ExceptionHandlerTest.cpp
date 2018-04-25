/**
 * @file ExceptionHandlerTest.cpp
 * @author seckler
 * @date 17.04.18
 */

#include <gtest/gtest.h>

#include <omp.h>
#include "ExceptionHandlerTest.h"
#include "utils/ExceptionHandler.h"

using autopas::utils::ExceptionHandler;
using autopas::utils::ExceptionBehavior;


TEST_F(ExceptionHandlerTest, TestDefault) {
  EXPECT_ANY_THROW(ExceptionHandler::exception("testignore"));
  EXPECT_ANY_THROW(ExceptionHandler::exception(std::exception()));
  ExceptionHandler::setBehavior(ExceptionBehavior::printCustomAbortFunction);
  EXPECT_DEATH(ExceptionHandler::exception("testignore"), "");
  EXPECT_DEATH(ExceptionHandler::exception(std::exception()), "");
}


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

#ifdef _OPENMP
TEST_F(ExceptionHandlerTest, TestThreadSafe) {
  ASSERT_GT(omp_get_max_threads(), 1);
#pragma omp parallel
  {
    EXPECT_GT(omp_get_num_threads(), 1);
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::exception("testignore");
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::exception("testignore2");
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::exception("testignore3");
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::exception("testignore4");
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::exception("testignore5");
  };
#pragma omp parallel
  {
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
    ExceptionHandler::setBehavior(ExceptionBehavior::throwException);
    ExceptionHandler::setBehavior(ExceptionBehavior::printCustomAbortFunction);
    auto abortFunction = []() -> void {
      AutoPasLogger->error("TESTABORTCUSTOMCALL123");
      abort();
    };
    ExceptionHandler::setCustomAbortFunction(abortFunction);
    ExceptionHandler::setBehavior(ExceptionBehavior::throwException);
    ExceptionHandler::setBehavior(ExceptionBehavior::ignore);
  };
  ExceptionHandler::setBehavior(ExceptionBehavior::throwException);
}
#endif