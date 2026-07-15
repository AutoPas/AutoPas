/**
 * @file SortingThresholdInfoInterface.h
 * @date 14 Jul 2026
 */

#pragma once

namespace autopas {

/**
 * Polymorphic base for sorting-threshold storage.
 *
 * Lets CellFunctor (2-body) and CellFunctor3B (3-body) each declare setters that accept some threshold info
 * without being textually coupled to a concrete shape that belongs to the other. CellFunctor dynamic_casts the
 * received reference down to the concrete type it knows how to consume (SortingThresholdInfo2B).
 */
struct SortingThresholdInfoInterface {
  virtual ~SortingThresholdInfoInterface() = default;
};

}  // namespace autopas
