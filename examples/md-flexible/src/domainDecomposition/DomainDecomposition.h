/**
 * @file DomainDecomposition.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

#include <array>
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
   * Returns the index of the local domain in the global domain context.
   */
  virtual const int getDomainIndex() = 0;

  /**
   * Returns the number of dimensions in this domain decomposition.
   */
  virtual const int getDimensionCount() = 0;

  /**
   * Returns the minimum coordinates of the global domain.
   */
  virtual const std::vector<double> getGlobalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the global domain.
   */
  virtual const std::vector<double> getGlobalBoxMax() = 0;

  /**
   * Returns the minimum coordinates of the local domain.
   */
  virtual const std::vector<double> getLocalBoxMin() = 0;

  /**
   * Returns the maximum coordinates of the local domain.
   */
  virtual const std::vector<double> getLocalBoxMax() = 0;

  /**
   * Checks if provided cooridnates are within the local domain.
   */
  virtual bool isInsideLocalDomain(const std::vector<double> &coordinates) = 0;

  /**
   * Checks if the provided coordinates are located in the local domain.
   * Instead of a vector, the coordinates are of type std::array<double, 3> to be compatible with AutoPas.
   */
  virtual bool isInsideLocalDomain(const std::array<double, 3> &coordinates) = 0;
};
