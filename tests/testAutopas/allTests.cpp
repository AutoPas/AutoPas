/**
 * @file allTests.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>

#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
#include <memory>

#include "autopas/utils/PythonInterpreter.h"

/**
 * Keeps the shared Python interpreter alive for all tests. Python cannot be restarted easily after it is finalized, so
 * all tests have to use this one.
 */
class PythonInterpreterEnvironment : public testing::Environment {
 public:
  void SetUp() override { _pythonInterpreter = autopas::utils::PythonInterpreter::getSharedInterpreter(); }
  void TearDown() override { _pythonInterpreter.reset(); }

 private:
  /**
   * The interpreter, kept alive from before the first test until after the last one.
   */
  std::shared_ptr<autopas::utils::PythonInterpreter> _pythonInterpreter;
};
#endif

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";
#ifdef AUTOPAS_ENABLE_PYTHON_BASED_TUNING
  // gtest takes ownership of the environment
  testing::AddGlobalTestEnvironment(new PythonInterpreterEnvironment);
#endif

  // damit nur meine Tests durchgelaufen werden
  return RUN_ALL_TESTS();
}
