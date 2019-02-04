/**
 * @file BlackBoxOptions.h
 * @author seckler
 * @date 04.02.19
 */

#pragma once

namespace autopas {

/**
 * Possible choices for the cell pair traversal.
 */
enum BlackBoxTraversalOption {
  normal,  ///< iterate over everything, this is the default option if blackbox is deactivated
  inner,   ///< iterate only over the inner parts of the domain
  outer    ///< iterate over the outer parts of the domain
};

}  // namespace autopas