/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <array>

#include "autopas/AutoPas.h"
#include "src/TypeDefinitions.h"

/**
 * An interface for domain decompositions which can be used in the simulation
 */
class DomainDecomposition {
 public:
  /**
   * Destructor.
   */
  virtual ~DomainDecomposition() = default;

  /**
   * Type for the AutoPas container
   */
  using SharedAutoPasContainer = std::shared_ptr<autopas::AutoPas<ParticleType>>;

  /**
   * Updates the domain decomposition to the current topology.
   */
  virtual void update(SharedAutoPasContainer &autoPasContainer, const double &work) = 0;

  /**
   * Returns the index of the local domain in the global domain context.
   * @return domain index.
   */
  virtual const int getDomainIndex() = 0;

  /**
   * Returns the minimum coordinates of the global domain.
   * @return bottom left front corner of the global domain.
   */
  virtual const std::array<double, 3> getGlobalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the global domain.
   * @return top right back corner of the global domain.
   */
  virtual const std::array<double, 3> getGlobalBoxMax() = 0;

  /**
   * Returns the minimum coordinates of the local domain.
   * @return bottom left front corner of the local domain.
   */
  virtual const std::array<double, 3> getLocalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the local domain.
   * @return top right back corner of the local domain.
   */
  virtual const std::array<double, 3> getLocalBoxMax() = 0;

  /**
   * Checks if the provided coordinates are located in the local domain.
   * @param coordinates: The coordinates in question.
   * @return true if the coordinates lie inside the local domain, false otherwise.
   */
  virtual bool isInsideLocalDomain(const std::array<double, 3> &coordinates) = 0;
};
