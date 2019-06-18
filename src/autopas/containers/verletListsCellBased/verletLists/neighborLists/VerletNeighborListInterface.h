/**
 * @file VerletNeighborListInterface.h
 *
 * @date 20.05.2019
 * @author humig
 */

#pragma once

#include <vector>
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {

/**
 * Interface for neighbor lists used by VarVerletLists.
 * @tparam Particle The particle type this neighbor list uses.
 */
template <class Particle>
class VerletNeighborListInterface {
 public:
  /**
   * virtual default destructor
   */
  virtual ~VerletNeighborListInterface() = default;

  /**
   * Returns the ContainerOption this neighbor list is for.
   * @return the ContainerOption this neighbor list is for.
   */
  virtual ContainerOption getContainerType() const = 0;

  /**
   * Builds the neighbor list from a LinkedCells object. This only builds the AoS.
   * @param linkedCells The linked cells to use for building the neighbor list.
   * @param useNewton3 If true, use newton 3 for the neighbor list.
   */
  virtual void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                             typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                                 bool useNewton3) = 0;

  /**
   * Generates the SoA from the AoS. @see buildNeighborList()
   */
  virtual void generateSoAFromAoS() = 0;

  /**
   * Returns whether the SoA is build and up to date with the AoS.
   * @return True, if the SoA is up to date with the AoS, false otherwise.
   */
  virtual bool isSoAListValid() const = 0;

  /**
   * Returns the number of neighbor pairs in the list.
   * @return the number of neighbor pairs in the list.
   */
  virtual long getNumberOfNeighborPairs() const = 0;
};

}  // namespace autopas