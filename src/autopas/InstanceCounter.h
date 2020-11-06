/**
 * @file InstanceCounter.h
 * @author seckler
 * @date 06.11.20
 */

#pragma once

namespace autopas {

/**
 * Class to count autopas instances.
 */
struct InstanceCounter {
  /**
   * Instance counter to help track the number of autopas instances. Needed for correct management of the logger.
   */
  static inline unsigned int count{0};
};
}  // namespace autopas