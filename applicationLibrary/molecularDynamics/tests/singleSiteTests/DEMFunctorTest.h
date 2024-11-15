/**
 * @file DEMFunctorTest.h
 * @author Joon Kim
 * @date 15.11.2024
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"
#include "discreteElementMethodLibrary/DEMFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class DEMFunctorTest : public AutoPasTestBase {
 public:
  DEMFunctorTest() : AutoPasTestBase() {}

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
// typedefs
template <bool mixing, bool globals>
using DEMFunctorGranular = demLib::DEMFunctor<GranularParticle, mixing, autopas::FunctorN3Modes::Both, globals>;

// struct aliasing for readable names
struct DEMFunNoMixNoGlob : public DEMFunctorGranular<false, false> {
  using DEMFunctorGranular<false, false>::DEMFunctor;
};
struct DEMFunMixNoGlob : public DEMFunctorGranular<true, false> {
  using DEMFunctorGranular<true, false>::DEMFunctor;
};
struct DEMFunNoMixGlob : public DEMFunctorGranular<false, true> {
  using DEMFunctorGranular<false, true>::DEMFunctor;
};
struct DEMFunMixGlob : public DEMFunctorGranular<true, true> {
  using DEMFunctorGranular<true, true>::DEMFunctor;
};