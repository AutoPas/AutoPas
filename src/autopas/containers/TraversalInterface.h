/**
 * @file TraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all traversals.
 */
 template <InteractionTypeOption::Value interactionType>
class TraversalInterface {
 public:
  /**
   * Destructor of TraversalInterface
   */
  virtual ~TraversalInterface() = default;

  /**
   * Return a enum representing the name of the traversal class.
   * @return Enum representing traversal.
   */
  [[nodiscard]] virtual TraversalOption getTraversalType() const = 0;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  [[nodiscard]] virtual bool isApplicable() const = 0;

  /**
   * Return whether the traversal uses newton 3.
   * @return true iff newton 3 is used.
   */
  [[nodiscard]] virtual bool getUseNewton3() const = 0;

  /**
   * Return the data layout option
   * @return the data layout option
   */
  [[nodiscard]] virtual DataLayoutOption getDataLayout() const = 0;

  /**
   * Initializes the traversal. Should be called before traverse().
   */
  virtual void initTraversal() = 0;

  /**
   * Finalizes the traversal. Should be called after traverse().
   */
  virtual void endTraversal() = 0;

  /**
   * Traverses all particle pairs.
   */
  template <InteractionTypeOption::Value type = interactionType , std::enable_if<type == InteractionTypeOption::pairwise>> void traverseParticlePairs();

  template <std::size_t r = nrows, std::size_t c = ncols>
  typename std::enable_if<r == c>::type setIdentity ()
  /**
   * Traverses all particle triplets.
   */
  virtual void traverseParticleTriplets() = 0;
};

}  // namespace autopas
