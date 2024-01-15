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

  /**
   * Constructor of the TraversalBase class.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
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
  /**
   * The datalayout used by this traversal.
   */
  const DataLayoutOption::Value _dataLayout;

  /**
   * If this traversal makes use of newton3.
   */
  const bool _useNewton3;
};

}  // namespace autopas
