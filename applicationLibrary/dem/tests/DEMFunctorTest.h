/**
 * @file DEMFunctor.h
 * @author Joon Kim
 * @date 25.06.25
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"
#include "dem/demLibrary/DEMFunctor.h"
#include "dem/demLibrary/DEMParameters.h"
#include "dem/demLibrary/GranularDEMParticle.h"
#include "testingHelpers/commonTypedefs.h"

class DEMFunctorTest : public AutoPasTestBase {
 public:
  DEMFunctorTest() : AutoPasTestBase() {}

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

/**
 * Short for the AutoPas DEM particle type.
 */
using Granular = demLib::GranularDEMParticle;

// typedefs to hide clutter
template <bool mixing>
using DEMFunGran = demLib::DEMFunctor<Granular, mixing, autopas::FunctorN3Modes::Both, false, false>;

// struct aliasing for readable names
struct DEMFunMixGlob : public DEMFunGran<true> {
  using DEMFunGran<true>::DEMFunctor;
};

struct DEMFunNoMixGlob : public DEMFunGran<false> {
  using DEMFunGran<false>::DEMFunctor;
};