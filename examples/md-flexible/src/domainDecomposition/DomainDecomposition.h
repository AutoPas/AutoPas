/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <vector>

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
   * Updates the domain decomposition to the current topology.
   */
  virtual void update() = 0;

  /**
   * Returns the number of dimensions in this domain decomposition.
   */
  virtual const int getDimensionCount() = 0;

  /**
   * Returns the minimum coordinates of the global domain.
   */
  virtual std::vector<double> getGlobalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the global domain.
   */
  virtual std::vector<double> getGlobalBoxMax() = 0;

  /**
   * Returns the minimum coordinates of the local domain.
   */
  virtual std::vector<double> getLocalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the local domain.
   */
  virtual std::vector<double> getLocalBoxMax() = 0;

  /**
   * Checks if provided cooridnates are within the local domain.
   */
  virtual bool isInsideLocalDomain(const std::vector<double> &coordinates) = 0;
};
