/**
 * @file TraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

#include "autopas/options/TraversalOptions.h"

namespace autopas {

/**
 * Provides a way to iterate over the possible choices of TraversalOption.
 */
static std::vector<TraversalOptions> allTraversalOptions = {TraversalOptions::c08,
                                                            TraversalOptions::sliced,
                                                            TraversalOptions::c18,
                                                            TraversalOptions::c01,
                                                            TraversalOptions::directSumTraversal,
                                                            TraversalOptions::slicedVerlet,
                                                            TraversalOptions::c18Verlet,
                                                            TraversalOptions::c01Verlet};

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
