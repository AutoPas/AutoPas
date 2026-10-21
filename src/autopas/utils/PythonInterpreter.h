/**
 * @file PythonInterpreter.h
 * @author S. Newcome
 * @date 21.10.2026
 */

#pragma once

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
   * Starts the embedded Python interpreter, or attaches to one that is already running in this process. Only an
   * interpreter started here is finalized on destruction.
   *
   * Causes an exception if AutoPas was compiled with AUTOPAS_ENABLE_PYTHON_BASED_TUNING=OFF.
   */
  PythonInterpreter();

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

 private:
  /**
   * Whether the constructor started the interpreter, so that the destructor has to finalize it.
   */
  [[maybe_unused]] bool _thisClassStartedInterpreter{false};
};

}  // namespace autopas::utils