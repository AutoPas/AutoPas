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

void ExceptionHandlerTest::SetUp() {
  autopas::logger::create();
}

void ExceptionHandlerTest::TearDown() {
  autopas::logger::unregister();
  // reset to default values
  ExceptionHandler::setBehavior(ExceptionBehavior::throwException);
  ExceptionHandler::setCustomAbortFunction(abort);
}


TEST_F(ExceptionHandlerTest, TestThrowCustom) {
  EXPECT_THROW(ExceptionHandler::exception(std::runtime_error("runtimeerror")), std::runtime_error);
}

TEST_F(ExceptionHandlerTest, TestDefault) {
  EXPECT_THROW(ExceptionHandler::exception("testthrow"), ExceptionHandler::AutoPasException);
  EXPECT_THROW(ExceptionHandler::exception(std::exception()), std::exception);
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

  EXPECT_THROW(ExceptionHandler::exception("testignore"), ExceptionHandler::AutoPasException);

  EXPECT_THROW(ExceptionHandler::exception(std::exception()), std::exception);
}

TEST_F(ExceptionHandlerTest, TestAbort) {
  ExceptionHandler::setBehavior(ExceptionBehavior::printAbort);

  EXPECT_DEATH(ExceptionHandler::exception("testignore"), "");

  EXPECT_DEATH(ExceptionHandler::exception(std::exception()), "");
}

TEST_F(ExceptionHandlerTest, TestAbortCustom) {

  auto abortFunction = []() -> void {
    AutoPasLogger->error("TESTABORTCUSTOMCALL123");
    abort();
  };
  ExceptionHandler::setBehavior(ExceptionBehavior::printCustomAbortFunction);

  ExceptionHandler::setCustomAbortFunction(abortFunction);

  EXPECT_DEATH(ExceptionHandler::exception("testignore"), "");

  EXPECT_DEATH(ExceptionHandler::exception(std::exception()), "");
}

TEST_F(ExceptionHandlerTest, TestTryRethrow) {
  try{
    throw std::runtime_error("me throwing things");
  } catch(std::exception& e){
    EXPECT_THROW(ExceptionHandler::rethrow(), std::runtime_error);
  }
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