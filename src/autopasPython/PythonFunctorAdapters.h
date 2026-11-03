/**
 * @file PythonFunctorAdapters.h
 * @date 03.11.2026
 *
 * The CRTP functors which AutoPas sees for functors written in Python.
 */
#pragma once

#include <pybind11/pybind11.h>

#include <atomic>
#include <exception>
#include <mutex>
#include <string>
#include <type_traits>
#include <utility>

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/baseFunctors/TriwiseFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "src/TypeDefinitions.h"

namespace autopasPython {

/**
 * The connection from the OpenMP threads of a traversal back into the interpreter.
 *
 * All Python code runs with the GIL acquired, per kernel call, which also makes writing to particles from Python safe.
 * Traversals are started on one thread for performance, see PythonFunctor::compute(). This class does not rely on that
 * and the thread which starts a traversal must not hold the GIL, as other threads would deadlock otherwise.
 *
 * Exceptions must not leave an OpenMP region, as this terminates the process. So the first exception is stored,
 * all further kernel calls are skipped, and rethrow() raises it after the traversal has finished.
 */
class PythonDispatcher {
 public:
  /**
   * Fetches the Python methods from the object.
   * @note They are only held during a traversal, as a bound method refers to its object. Keeping it in here would
   * form a reference cycle which Python's garbage collector cannot see.
   * @param self The Python functor object.
   */
  void acquire(const pybind11::object &self) {
    const pybind11::gil_scoped_acquire gil{};
    _aos = self.attr("aos");
    _initTraversal = self.attr("init_traversal");
    _endTraversal = self.attr("end_traversal");
  }

  /**
   * Releases the methods fetched by acquire().
   */
  void release() {
    const pybind11::gil_scoped_acquire gil{};
    _aos = pybind11::none();
    _initTraversal = pybind11::none();
    _endTraversal = pybind11::none();
  }

  /**
   * Raises the exception caught during the traversal, if there was one.
   */
  void rethrow() {
    if (_failed.exchange(false)) {
      const std::scoped_lock lock{_mutex};
      std::rethrow_exception(std::exchange(_exception, nullptr));
    }
  }

  /**
   * Calls the kernel, passing particles by reference.
   * @tparam Args Types of the arguments.
   * @param args Particles and the newton3 flag.
   */
  template <class... Args>
  void callAoS(Args &&...args) {
    guarded([&]() { _aos(wrap(args)...); });
  }

  /**
   * Calls Python's init_traversal().
   */
  void callInitTraversal() {
    guarded([&]() { _initTraversal(); });
  }

  /**
   * Calls Python's end_traversal().
   * @param newton3
   */
  void callEndTraversal(bool newton3) {
    guarded([&]() { _endTraversal(newton3); });
  }

 private:
  /**
   * Particles are handed over without a copy, everything else by value.
   */
  template <class T>
  static pybind11::object wrap(T &value) {
    if constexpr (std::is_same_v<std::remove_cv_t<T>, ParticleType>) {
      return pybind11::cast(&value, pybind11::return_value_policy::reference);
    } else {
      return pybind11::cast(value);
    }
  }

  /**
   * Runs a Python call with the GIL held and without letting an exception escape.
   */
  template <class F>
  void guarded(F &&call) {
    if (_failed.load(std::memory_order_relaxed)) {
      return;
    }
    const pybind11::gil_scoped_acquire gil{};
    try {
      call();
    } catch (...) {
      const std::scoped_lock lock{_mutex};
      if (not _exception) {
        _exception = std::current_exception();
      }
      _failed.store(true);
    }
  }

  pybind11::object _aos{};
  pybind11::object _initTraversal{};
  pybind11::object _endTraversal{};
  std::atomic<bool> _failed{false};
  std::mutex _mutex{};
  std::exception_ptr _exception{};
};

/**
 * The properties of a Python functor, which AutoPas queries without the GIL.
 */
struct PythonFunctorTraits {
  /**
   * @copydoc FunctorHandle::name()
   */
  std::string name;
  /**
   * @copydoc FunctorHandle::allowsNewton3()
   */
  bool allowsNewton3;
  /**
   * @copydoc FunctorHandle::allowsNonNewton3()
   */
  bool allowsNonNewton3;
  /**
   * If the AutoTuner should measure this functor.
   */
  bool relevantForTuning;
};

/**
 * Pairwise functor which forwards to Python. Only provides the AoS kernel.
 */
class PythonPairwiseAdapter : public autopas::PairwiseFunctor<ParticleType, PythonPairwiseAdapter> {
 public:
  /**
   * Constructor.
   * @param cutoff Pairs further apart are filtered out before entering Python.
   * @param traits
   * @param dispatcher
   */
  PythonPairwiseAdapter(double cutoff, PythonFunctorTraits traits, PythonDispatcher &dispatcher)
      : autopas::PairwiseFunctor<ParticleType, PythonPairwiseAdapter>(cutoff),
        _cutoffSquared(cutoff * cutoff),
        _traits(std::move(traits)),
        _dispatcher(dispatcher) {}

  void AoSFunctor(ParticleType &i, ParticleType &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const auto dr = i.getR() - j.getR();
    if (autopas::utils::ArrayMath::dot(dr, dr) > _cutoffSquared) {
      return;
    }
    _dispatcher.callAoS(i, j, newton3);
  }

  void initTraversal() final { _dispatcher.callInitTraversal(); }

  void endTraversal(bool newton3) final { _dispatcher.callEndTraversal(newton3); }

  bool allowsNewton3() final { return _traits.allowsNewton3; }

  bool allowsNonNewton3() final { return _traits.allowsNonNewton3; }

  bool isVecPatternAllowed(const autopas::VectorizationPatternOption::Value) final { return true; }

  bool isRelevantForTuning() final { return _traits.relevantForTuning; }

  std::string getName() final { return _traits.name; }

 private:
  double _cutoffSquared;
  PythonFunctorTraits _traits;
  PythonDispatcher &_dispatcher;
};

/**
 * Triwise functor which forwards to Python. Only provides the AoS kernel.
 */
class PythonTriwiseAdapter : public autopas::TriwiseFunctor<ParticleType, PythonTriwiseAdapter> {
 public:
  /**
   * Constructor.
   * @param cutoff Triplets with any two particles further apart are filtered out before entering Python.
   * @param traits
   * @param dispatcher
   */
  PythonTriwiseAdapter(double cutoff, PythonFunctorTraits traits, PythonDispatcher &dispatcher)
      : autopas::TriwiseFunctor<ParticleType, PythonTriwiseAdapter>(cutoff),
        _cutoffSquared(cutoff * cutoff),
        _traits(std::move(traits)),
        _dispatcher(dispatcher) {}

  void AoSFunctor(ParticleType &i, ParticleType &j, ParticleType &k, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy() or k.isDummy()) {
      return;
    }
    const auto drIJ = i.getR() - j.getR();
    const auto drJK = j.getR() - k.getR();
    const auto drKI = k.getR() - i.getR();
    if (autopas::utils::ArrayMath::dot(drIJ, drIJ) > _cutoffSquared or
        autopas::utils::ArrayMath::dot(drJK, drJK) > _cutoffSquared or
        autopas::utils::ArrayMath::dot(drKI, drKI) > _cutoffSquared) {
      return;
    }
    _dispatcher.callAoS(i, j, k, newton3);
  }

  void initTraversal() final { _dispatcher.callInitTraversal(); }

  void endTraversal(bool newton3) final { _dispatcher.callEndTraversal(newton3); }

  bool allowsNewton3() final { return _traits.allowsNewton3; }

  bool allowsNonNewton3() final { return _traits.allowsNonNewton3; }

  bool isVecPatternAllowed(const autopas::VectorizationPatternOption::Value) final { return true; }

  bool isRelevantForTuning() final { return _traits.relevantForTuning; }

  std::string getName() final { return _traits.name; }

 private:
  double _cutoffSquared;
  PythonFunctorTraits _traits;
  PythonDispatcher &_dispatcher;
};

}  // namespace autopasPython
