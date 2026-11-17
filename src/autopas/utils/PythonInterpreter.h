/**
 * @file PythonInterpreter.h
 * @author S. Newcome
 * @date 21.10.2026
 */

#pragma once

#include <memory>

namespace autopas::utils {

/**
 * Holds and keeps alive the embedded Python interpreter, used for Python-based tuning strategies.
 *
 * This class is essentially a wrapper around the interpreter, so we don't need to expose pybind11 functionality to the
 * consumer of AutoPas, if they don't require it, and so the interface, regardless of if
 * AUTOPAS_ENABLE_PYTHON_BASED_TUNING is set ON or OFF, remains the same.
 */
class PythonInterpreter {
 public:
  /**
   * Destructor. Finalizes the interpreter only if the constructor started it. An interpreter that was already running
   * is left to its owner.
   */
  ~PythonInterpreter();

  /**
   * Copy constructor, deleted. There is only one interpreter per process.
   */
  PythonInterpreter(const PythonInterpreter &other) = delete;

  /**
   * Copy assignment operator, deleted. There is only one interpreter per process.
   * @return Reference to this.
   */
  PythonInterpreter &operator=(const PythonInterpreter &other) = delete;

  /**
   * Move constructor, deleted. The interpreter is shared through a std::shared_ptr, so it never needs to move.
   */
  PythonInterpreter(PythonInterpreter &&other) = delete;

  /**
   * Move assignment operator, deleted. Assigning over an object that started the interpreter would finalize it while
   * others still use it.
   * @return Reference to this.
   */
  PythonInterpreter &operator=(PythonInterpreter &&other) = delete;

  /**
   * Get the interpreter that all callers in this process share, and create it if none is alive. It is finalized when
   * the last caller releases it.
   *
   * To use Python in multiple AutoPas objects or for some separate use case, one interpreter must be kept alive.
   * AutoPas internally handles this if the multiple AutoPas objects are kept alive at the same time. But if creating
   * and destroying AutoPas objects serially (e.g. in tests), an interpreter must be otherwise kept alive independently
   * of AutoPas, by holding the shared pointer returned here.
   *
   * Causes an exception in the cases listed at the constructor.
   * @return Shared pointer to the interpreter.
   */
  static std::shared_ptr<PythonInterpreter> getSharedInterpreter();

 private:
  /**
   * Starts the embedded Python interpreter, or attaches to one that is already running in this process. Only an
   * interpreter started here is finalized on destruction.
   *
   * Causes an exception if AutoPas was compiled with AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF, or if a PythonInterpreter
   * already finalized the interpreter in this process.
   *
   * Private, so that the only way to get an interpreter is getSharedInterpreter().
   */
  PythonInterpreter();

  /**
   * Whether the constructor started the interpreter, so that the destructor has to finalize it.
   */
  [[maybe_unused]] bool _thisClassStartedInterpreter{false};
};

}  // namespace autopas::utils