/**
 * @file KryptonExtendedATMFunctorTest.h
 * @author D. Martin
 * @date 15.12.25
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"
#include "molecularDynamicsLibrary/KryptonExtendedATMFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class KryptonExtendedATMFunctorTest : public AutoPasTestBase {
 public:
  KryptonExtendedATMFunctorTest() : AutoPasTestBase() {}

  enum InteractionType { own, pair12, pair21, triple, verlet };
  // Where to place 3 particles. Inside or outside the domain.
  enum where_type { allInside, ininout, inoutout, allOutside };

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

// typedefs to hide clutter
template <bool globals>
using KryptonATMFunMol = mdLib::KryptonExtendedATMFunctor<Molecule, autopas::FunctorN3Modes::Both, globals>;

// struct aliasing for readable names
struct KryptonATMFunNoMixNoGlob : public KryptonATMFunMol<false> {
  using KryptonATMFunMol<false>::KryptonExtendedATMFunctor;
};
struct KryptonATMFunNoMixGlob : public KryptonATMFunMol<true> {
  using KryptonATMFunMol<true>::KryptonExtendedATMFunctor;
};
