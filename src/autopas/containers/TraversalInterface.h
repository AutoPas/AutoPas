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
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  [[nodiscard]] virtual bool isApplicable() const = 0;

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

  /**
   * This function os used in case halo region length is different from interactionLength, i.e. for Hierarchical Grid's
   * lower levels.
   * @param haloRegionLength length of the halo region, it is usually cutoff + skin
   */
  virtual void setHaloRegionLength(const double haloRegionLength) {}

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
