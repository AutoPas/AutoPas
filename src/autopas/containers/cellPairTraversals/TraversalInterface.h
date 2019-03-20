/**
 * @file TraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {

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
  virtual TraversalOption getTraversalType() = 0;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  virtual bool isApplicable() = 0;

  /**
   * Returns the data Layout that has to be present in order to use the Traversal
   * @return Required DataLayout
   */
  virtual DataLayoutOption requiredDataLayout() = 0;
};

}  // namespace autopas
