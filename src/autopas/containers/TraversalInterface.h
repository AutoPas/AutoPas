/**
 * @file TraversalInterface.h
 * @author F. Gratl
 * @date 7/3/18
 */

#pragma once

#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {

/**
 * This interface serves as a common parent class for all traversals.
 */
class TraversalInterface {
 public:
  /**
   * Destructor of TraversalInterface
   */
  virtual ~TraversalInterface() = default;

  /**
   * Constructor of the TraversalInterface.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  TraversalInterface(DataLayoutOption dataLayout, bool useNewton3) : _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Return a enum representing the name of the traversal class.
   * @return Enum representing traversal.
   */
  [[nodiscard]] virtual TraversalOption getTraversalType() const = 0;

  /**
   * Checks if the traversal is applicable to the current state of the domain. This is designed to only return false if
   * the traversal is inapplicable for domain-related issues. Domain-independent factors in a traversal being
   * inapplicable should be handled by @link Configuration::hasCompatibleValues.
   * @return true iff the traversal is applicable to the domain.
   */
  [[nodiscard]] virtual bool isApplicableToDomain() const = 0;

  /**
   * Initializes the traversal. Should be called before traverse().
   */
  virtual void initTraversal() = 0;

  /**
   * Traverse the particles by pairs, triplets etc. as determined by the Functor type.
   */
  virtual void traverseParticles() = 0;

  /**
   * Finalizes the traversal. Should be called after traverse().
   */
  virtual void endTraversal() = 0;

  /**
   * Return whether the traversal uses newton 3.
   * @return true iff newton 3 is used.
   */
  [[nodiscard]] bool getUseNewton3() const { return _useNewton3; }

  /**
   * Return the data layout option
   * @return the data layout option
   */
  [[nodiscard]] DataLayoutOption getDataLayout() const { return _dataLayout; }

 protected:
  /**
   * The datalayout used by this traversal.
   */
  DataLayoutOption _dataLayout;

  /**
   * If this traversal makes use of newton3.
   */
  bool _useNewton3;
};
}  // namespace autopas
