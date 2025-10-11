/**
 * @file ATMFunctorTest.h
 * @author muehlhaeusser
 * @date 29.08.23
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ExceptionHandler.h"
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class ATMFunctorTest : public AutoPasTestBase {
 public:
  ATMFunctorTest() : AutoPasTestBase() {}

  enum SoAFunctorType { single, pair12, pair21, triple, verlet };
  static std::string to_string(SoAFunctorType type) {
    switch (type) {
      case SoAFunctorType::single:
        return "single";
      case SoAFunctorType::pair12:
        return "pair12";
      case SoAFunctorType::pair21:
        return "pair21";
      case SoAFunctorType::triple:
        return "triple";
      case SoAFunctorType::verlet:
        return "verlet";
      default:
        return "unknown";
    }
  }
  // Where to place 3 particles. Inside or outside the domain.
  enum where_type { allInside, ininout, inoutout, allOutside };
  static std::string to_string(where_type where) {
    switch (where) {
      case allInside:
        return "all inside";
      case ininout:
        return "in, in, out";
      case inoutout:
        return "in, out, out";
      case allOutside:
        return "all outside";
      default:
        return "unknown";
    }
  }
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
template <bool mixing, bool globals>
using ATMFunMol = mdLib::AxilrodTellerMutoFunctor<Molecule, mixing, autopas::FunctorN3Modes::Both, globals>;

// struct aliasing for readable names
struct ATMFunMixNoGlob : public ATMFunMol<true, false> {
  using ATMFunMol<true, false>::AxilrodTellerMutoFunctor;
};
struct ATMFunNoMixNoGlob : public ATMFunMol<false, false> {
  using ATMFunMol<false, false>::AxilrodTellerMutoFunctor;
};
struct ATMFunMixGlob : public ATMFunMol<true, true> {
  using ATMFunMol<true, true>::AxilrodTellerMutoFunctor;
};
struct ATMFunNoMixGlob : public ATMFunMol<false, true> {
  using ATMFunMol<false, true>::AxilrodTellerMutoFunctor;
};
