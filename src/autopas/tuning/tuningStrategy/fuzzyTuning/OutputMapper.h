/**
 * @file OutputMapper.h
 * @author Manuel Lerchner
 * @date 29.05.24
 */

#pragma once

#include <map>
#include <memory>
#include <string>

namespace autopas::fuzzy_logic {

/**
 * Used to store the mapping information to transform the result of the fuzzy controller into a configuration suitable
 * for AutoPas.
 */
class OutputMapper {
 public:
  /**
   * Constructs a ConfigurationMapper with empty dimensions.
   */
  OutputMapper() = default;
};

}  // namespace autopas::fuzzy_logic
