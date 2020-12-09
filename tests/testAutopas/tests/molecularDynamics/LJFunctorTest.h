/**
 * @file LJFunctorTest.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/utils/ExceptionHandler.h"
#include "testingHelpers/commonTypedefs.h"

namespace LJFunctorTest {

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

// typedefs to hide clutter
template <bool shift, bool mixing, bool globals>
using LJFunMol = autopas::LJFunctor<Molecule, shift, mixing, autopas::FunctorN3Modes::Both, globals>;
template <bool shift, bool mixing, bool globals>
using LJFunAVXMol = autopas::LJFunctorAVX<Molecule, shift, mixing, autopas::FunctorN3Modes::Both, globals>;

// struct aliasing for readable names
struct LJFunShiftMixNoGlob : public LJFunMol<true, true, false> {
  using LJFunMol<true, true, false>::LJFunctor;
};
struct LJFunShiftNoMixNoGlob : public LJFunMol<true, false, false> {
  using LJFunMol<true, false, false>::LJFunctor;
};
struct LJFunAVXShiftMixNoGlob : public LJFunAVXMol<true, true, false> {
  using LJFunAVXMol<true, true, false>::LJFunctorAVX;
};
struct LJFunAVXShiftNoMixNoGlob : public LJFunAVXMol<true, false, false> {
  using LJFunAVXMol<true, false, false>::LJFunctorAVX;
};
struct LJFunShiftMixGlob : public LJFunMol<true, true, true> {
  using LJFunMol<true, true, true>::LJFunctor;
};
struct LJFunShiftNoMixGlob : public LJFunMol<true, false, true> {
  using LJFunMol<true, false, true>::LJFunctor;
};
struct LJFunAVXShiftMixGlob : public LJFunAVXMol<true, true, true> {
  using LJFunAVXMol<true, true, true>::LJFunctorAVX;
};
struct LJFunAVXShiftNoMixGlob : public LJFunAVXMol<true, false, true> {
  using LJFunAVXMol<true, false, true>::LJFunctorAVX;
};
}  // end namespace LJFunctorTest
