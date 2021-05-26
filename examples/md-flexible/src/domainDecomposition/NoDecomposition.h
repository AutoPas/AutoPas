/**
 * @file KdTree.h
 * @author J. Körner
 * @date 25.05.2021
 */
#pragma once

#include "DomainDecomposition.h"

class NoDecomposition final : public DomainDecomposition {
 public:
  /**
  * Constructor.
  * @param argc The argument count passed to the main function.
  * @param argv The argument vector passed to the main function.
  * @param dimensionCount The number of dimensions for this domain decomposition.
  * @param globalBoxMin The minimum coordinates of the global domain.
  * @param globalBoxMax The maximum coordinates of the global domain.
  */
  NoDecomposition(int argc, char **argv, const int &dimensionCount, const std::vector<double> &globalBoxMin,
               const std::vector<double> &globalBoxMax);

  /**
   * Destructor.
   */
  virtual ~NoDecomposition();

  /**
   * Updates the domain decomposition to the current topology.
   */
  void update() override;

  /**
   * Returns the number of dimesnions in the domain decomposition.
   */
  const int getDimensionCount() override { return _dimensionCount; }

  /**
   * Returns the minimum coordinates of global domain.
   */
  std::vector<double> getGlobalBoxMin() override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of global domain.
   */
  std::vector<double> getGlobalBoxMax() override { return _globalBoxMax; }

  /**
   * Returns the minimum coordinates of local domain.
   */
  std::vector<double> getLocalBoxMin() override { return _globalBoxMin; }

  /**
   * Returns the maximum coordinates of local domain.
   */
  std::vector<double> getLocalBoxMax() override { return _globalBoxMax; }

  /**
   * Checks if the provided coordinates are located in the local domain.
   */
  bool isInsideLocalDomain(const std::vector<double> &coordinates) override;

 private:
  /**
   * The number of dimensions in this decomposition.
   */
  int _dimensionCount;

  /**
   * The minimum coordinates of the global domain.
   */
  std::vector<double> _globalBoxMin;

  /**
   * The maximum coordinates of the global domain.
   */
  std::vector<double> _globalBoxMax;

  /**
   * Initializes the global domain coordinates.
   */
  void initializeGlobalBox(const std::vector<double> &globalBoxMin, const std::vector<double> &globalBoxMax);
};
