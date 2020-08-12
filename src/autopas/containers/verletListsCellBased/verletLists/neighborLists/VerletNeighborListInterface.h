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
template <class Particle, class ParticleCell>
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
  [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

  /**
   * Builds the neighbor list from a LinkedCells object. This only builds the AoS.
   * @param linkedCells The linked cells to use for building the neighbor list.
   * @param useNewton3 If true, use newton 3 for the neighbor list.
   */
  virtual void buildAoSNeighborList(
      LinkedCells<Particle, typename VerletListHelpers<Particle>::PositionSoAArraysType> &linkedCells, bool useNewton3) = 0;

  /**
   * Checks if the neighbor list contains all pairs that is should.
   *
   * This is very costly, comparable to rebuilding it.
   * @param useNewton3 If the neighbor list should use newton 3.
   * @param cutoff The cutoff. Two particles that are further away than this distance are not considered.
   * @return If the current neighbor list is valid.
   */
  virtual bool checkNeighborListValidity(bool useNewton3, double cutoff) = 0;

  /**
   * Generates the SoA from the AoS. @see buildAoSNeighborList()
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