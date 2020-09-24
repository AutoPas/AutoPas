/**
 * @file ExecutionPolicies.h
 * @author seckler
 * @date 24.09.20
 */

#pragma once

namespace autopas {
/**
 * ExecutionPolicy defines how to execute the program.
 */
enum class ExecutionPolicy {
  /*
   * Sequential ExecutionPolicy
   */
  sequential,
  /*
   * Parallel ExecutionPolicy
   */
  parallel
};

}  // namespace autopas