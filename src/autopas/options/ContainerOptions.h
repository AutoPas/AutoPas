/**
 * @file ContainerOptions.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

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

}  // namespace autopas