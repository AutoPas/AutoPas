/**
 * @file SortingThresholdInfoInterface.h
 * @date 14.07.2026
 * @author hmeyran
 */

#pragma once

namespace autopas {

/**
 * Polymorphic base for sorting-threshold storage.
 *
 * Lets CellFunctor (2-body) and CellFunctor3B (3-body) each declare threshold info and getters based on configurations
 * without being textually coupled to a concrete shape that belongs to the other. CellFunctor dynamic_casts the
 * received reference down to the concrete type it knows how to consume (SortingThresholdInfo2B).
 */
struct SortingThresholdInfoInterface {
  /**
   * default destructor to make this a pure virtual struct.
   */
  virtual ~SortingThresholdInfoInterface() = default;
};

}  // namespace autopas
