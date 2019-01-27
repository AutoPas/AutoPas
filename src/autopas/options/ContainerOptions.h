/**
 * @file ContainerOptions.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <vector>

namespace autopas {

/**
 * Possible choices for the particle container type.
 */
enum ContainerOptions {
  directSum = 0,
  linkedCells = 1,
  verletLists = 2,
  verletListsCells = 3,
  verletClusterLists = 4,
};

/**
 * Provides a way to iterate over the possible choices of ContainerOption.
 */
static const std::vector<ContainerOptions> allContainerOptions = {
    ContainerOptions::directSum,        ContainerOptions::linkedCells,        ContainerOptions::verletLists,
    ContainerOptions::verletListsCells, ContainerOptions::verletClusterLists,
};

}  // namespace autopas