/**
 * @file CellPairTraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

#include <selectors/TraversalSelector.h>
namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum TraversalOptions {
  c08 = 0,
  sliced = 1,
};

/**
 * This interface serves as a common parent class for all cell pair traversals.
 */
class CellPairTraversalInterface {
 public:

  /**
   * Return a enum representing the name of the traversal class.
   */
  virtual TraversalOptions getName() = 0;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  virtual bool isApplicable() = 0;
};

}  // namespace autopas
