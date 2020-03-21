/**
 * @file LJFunctorTestGlobals.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"

class LJFunctorTest : public AutoPasTestBase {
 public:
  LJFunctorTest() : AutoPasTestBase() {}

  enum InteractionType { own, pair, verlet };
  enum where_type { inside, boundary, outside };

  /**
   * Checks if the given function throws an exception containing "not implemented".
   * @tparam FunType Type of the given function.
   * @param f Code to be checked as a lambda.
   * @return Empty string if nothing was caught, the exception string if a matching exception was found.
   * If the exception does not match it is rethrown.
   */
  template <class FunType>
  static std::string shouldSkipIfNotImplemented(FunType &&f) {
    try {
      f();
    } catch (const autopas::utils::ExceptionHandler::AutoPasException &e) {
      // if the functor fails with an exception about "not implemented" do not treat this as a failure but skip
      if (std::string(e.what()).find("not implemented") != std::string::npos) {
        return std::string(e.what());
      } else {
        throw;  // rethrow original exception
      }
    }
    return "";
  }
};
