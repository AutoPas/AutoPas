/**
 * @file TraversalBase.h
 * @author D. Martin
 * @date 8/1/24
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"

class TraversalInterface;

namespace autopas {

/**
 * This class serves as a common parent class for all traversals.
 */
class TraversalBase : public TraversalInterface {
 public:
  TraversalBase() = delete;

  TraversalBase(const DataLayoutOption::Value dataLayout, const bool useNewton3)
      : _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Return whether the traversal uses newton 3.
   * @return true iff newton 3 is used.
   */
  [[nodiscard]] bool getUseNewton3() const override { return _useNewton3; }

  /**
   * Return the data layout option
   * @return the data layout option
   */
  [[nodiscard]] DataLayoutOption getDataLayout() const override { return _dataLayout; }

 protected:
  const DataLayoutOption::Value _dataLayout;
  const bool _useNewton3;
};

}  // namespace autopas
