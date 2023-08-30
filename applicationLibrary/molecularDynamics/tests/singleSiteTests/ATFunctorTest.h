///**
// * @file ATFunctorTest.h
// * @author muehlhaeusser
// * @date 29.08.23
// */
//
//#pragma once
//
//#include <gtest/gtest.h>
//
//#include "AutoPasTestBase.h"
//#include "autopas/utils/ExceptionHandler.h"
//#include "molecularDynamicsLibrary/AxilrodTellerFunctor.h"
//#include "testingHelpers/commonTypedefs.h"
//
//class ATFunctorTest : public AutoPasTestBase {
// public:
//  ATFunctorTest() : AutoPasTestBase() {}
//
//  enum InteractionType { own, pair, verlet };
//  enum where_type { inside, boundary, outside };
//
//  /**
//   * Checks if the given function throws an exception containing "not implemented".
//   * @tparam FunType Type of the given function.
//   * @param f Code to be checked as a lambda.
//   * @return Empty string if nothing was caught, the exception string if a matching exception was found.
//   * If the exception does not match it is rethrown.
//   */
//  template <class FunType>
//  static std::string shouldSkipIfNotImplemented(FunType &&f) {
//    try {
//      f();
//    } catch (const autopas::utils::ExceptionHandler::AutoPasException &e) {
//      // if the functor fails with an exception about "not implemented" do not treat this as a failure but skip
//      if (std::string(e.what()).find("not implemented") != std::string::npos) {
//        return std::string(e.what());
//      } else {
//        throw;  // rethrow original exception
//      }
//    }
//    return "";
//  }
//};
//
//// typedefs to hide clutter
//template <bool shift, bool mixing, bool globals>
//using ATFunMol = mdLib::AxilrodTellerFunctor<Molecule, mixing, autopas::FunctorN3Modes::Both, globals>;
//
//// struct aliasing for readable names
//struct ATFunShiftMixNoGlob : public ATFunMol<true, true, false> {
//  using ATFunMol<true, true, false>::AxilrodTellerFunctor;
//};
//struct ATFunShiftNoMixNoGlob : public ATFunMol<true, false, false> {
//  using ATFunMol<true, false, false>::AxilrodTellerFunctor;
//};
//struct ATFunShiftMixGlob : public ATFunMol<true, true, true> {
//  using ATFunMol<true, true, true>::AxilrodTellerFunctor;
//};
//struct ATFunShiftNoMixGlob : public ATFunMol<true, false, true> {
//  using ATFunMol<true, false, true>::AxilrodTellerFunctor;
//};
