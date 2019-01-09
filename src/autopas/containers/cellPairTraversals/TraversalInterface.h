/**
 * @file CellPairTraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum TraversalOptions {
  dummyTraversal = -1,
  c08 = 0,
  sliced = 1,
  c18 = 2,
  c01 = 3,
  directSumTraversal = 4,
};

/**
 * This interface serves as a common parent class for all cell pair traversals.
 */
class TraversalInterface {
 public:
  /**
   * Destructor of CellPairTraversal.
   */
  virtual ~TraversalInterface() = default;

  /**
   * Return a enum representing the name of the traversal class.
   * @return Enum representing traversal.
   */
  virtual TraversalOptions getTraversalType() = 0;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  virtual bool isApplicable() = 0;
};

}  // namespace autopas
