/**
 * @file IterationMeasurements.h
 * @authors muehlhaeusser
 * @date 05.06.2024
 */

#pragma once

namespace autopas {
/**
 * Struct to collect all sorts of measurements taken during a computeInteractions iteration.
 */
struct IterationMeasurements {
  /// Time
  long timeIteratePairwise{};
  long timeRemainderTraversal{};
  long timeRebuild{};
  long timeTotal{};
  /// Energy. See RaplMeter.h for the meaning of each field.
  bool energyMeasurementsPossible{false};
  double energyPsys{};
  double energyPkg{};
  double energyRam{};
  long energyTotal{};
};
}  // namespace autopas
