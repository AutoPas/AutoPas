/**
 * @file RebuildDecisionContext.h
 * @author Manish Mishra
 * @date 25.11.2026
 */

#pragma once

#include <algorithm>
#include <climits>
#include <iostream>
#include <limits>
#include <ostream>
#include <utility>
#include <vector>

#include "autopas/utils/DynamicRebuildDecisionWithFPB/RebuildEstimateManager.h"
#include "autopas/utils/DynamicRebuildDecisionWithFPB/RemainderTraversalEstimateManager.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"

namespace autopas::utils::DynamicRebuildDecisionWithFPB {

/**
 * Includes all methods and variables needed to make a rebuild decision by minimizing
 * the sum of rebuild and remainder traversal times over the simulation time.
 */
class RebuildDecisionContext {
 public:
  /**
   * Default constructor.
   */
  RebuildDecisionContext() : _rebuildTimeEstimator(), _remainderTraversalTimePredictor(){};

  /**
   * For rebuildNeighborTimeEstimate logging.
   * @return The estimated rebuild time based on the current samples.
   */
  [[nodiscard]] double getRebuildNeighborTimeEstimate() const { return _rebuildNeighborTimeEstimate; }

  /**
   * For remainderTraversalTimeEstimate logging.
   * @return The estimated remainder traversal time based on the current samples.
   */
  [[nodiscard]] double getRemainderTraversalTimeEstimate() const { return _remainderTraversalTimeEstimate; }

  /**
   * For particlesBufferEstimate logging.
   * @return The estimated number of particles in buffer in the next iteration based on the current samples.
   */
  [[nodiscard]] size_t getNumParticlesBufferEstimate() const { return _numParticlesBufferEstimate; }

  /**
   * Called after a neighbor list rebuild to reset remainder traversal estimator state and
   * register the latest rebuild time.
   *
   * @param rebuildTime Measured remainder traversal time
   */
  void afterRebuild(const long rebuildTime) {
    _rebuildTimeEstimator.addNewPoint(rebuildTime);
    _remainderTraversalTimePredictor.reset();
    _remainderTraversalTimeEstimate = 0;
  }

  /**
   * Called whenever autotuner changes the algorithmic configuration to ensure that rebuild estimate is accurate for the
   * currently used algorithmic configuration.
   */
  void afterConfigurationChange() { _rebuildTimeEstimator.reset(); }

  /**
   * The function is called after the remainder traversal has been called for all the functors, i.e. at the start of
   * next iteration in updateContainer(). It adds a new data point to the remainder traversal time predictor.
   *
   * Additionally, tracks the increase in the number of particles in the buffer due to other mechanisms (e.g. periodic/
   * MPI boundaries, domain resizing, etc.) to support future buffer size estimation.
   *
   * @param remainderTraversalTime Measured total remainder traversal time (over all computeInteraction calls, e.g.
   * 2B+3B)
   * @param numParticlesBuffer     Current number of particles in the particle buffer
   */
  void afterRemainderTraversal(const long remainderTraversalTime, const size_t numParticlesBuffer) {
    _lastRemainderTraversalTime = remainderTraversalTime;

    if (numParticlesBuffer > 0) {
      _lastIncreaseNumParticlesBuffer = numParticlesBuffer - _afterFastAddNumParticlesBuffer;
      if (_lastIncreaseNumParticlesBuffer < 0) {
        throw autopas::utils::ExceptionHandler::AutoPasException(
            "Number of particles in buffer should never decrease except due to rebuilding.");
      }
    }
    _remainderTraversalTimePredictor.addNewPoint(static_cast<long>(numParticlesBuffer), remainderTraversalTime);
  }

  /**
   * Updates the estimated number of particles in the buffer for the next iteration based on the current number of
   * particles and the increase in number of particles in the buffer from the last iteration.
   * This is called after adding the fast particles to the buffer.
   * @param currentNumParticlesBuffer Current number of particles in the particle buffer
   */
  void updateNumParticlesBufferEstimate(const long currentNumParticlesBuffer) {
    _afterFastAddNumParticlesBuffer = currentNumParticlesBuffer;
    _numParticlesBufferEstimate = _afterFastAddNumParticlesBuffer + _lastIncreaseNumParticlesBuffer;
  }

  /**
   * Determines whether the particle buffer has grown large enough such that
   * performing a rebuild is more advantageous than executing another remainder traversal.
   *
   * @param currIter iteration since the last rebuild.
   * @return True if a dynamic rebuild should be triggered due to particle buffer fullness.
   */
  bool decideToRebuildOnParticleBufferFullness(const unsigned int currIter) {
    bool doDynamicRebuild{false};
    const auto valueRebuild = _rebuildTimeEstimator.predict();

    const auto valueRemainder =
        _remainderTraversalTimePredictor.predict(static_cast<long>(_numParticlesBufferEstimate));

    // for logging and testing
    _rebuildNeighborTimeEstimate = valueRebuild;
    _remainderTraversalTimeEstimate = valueRemainder;

    if (std::isnan(valueRebuild)) {
      throw autopas::utils::ExceptionHandler::AutoPasException("No rebuild samples available.");
    }

    if (valueRemainder > 0) {
      // The rebuild criteria is derived in the paper "To be updated".
      const double remainderIncline =
          (valueRemainder * currIter -
           static_cast<double>(_remainderTraversalTimePredictor.getSumRemainderTraversalTime()));
      if (remainderIncline >= valueRebuild) {
        doDynamicRebuild = true;
      }
    }
    return doDynamicRebuild;
  }

 private:
  /**
   * Handles storing past rebuild times samples and predicting next rebuild time based on simple arithmetic mean of the
   * collected samples.
   */
  RebuildEstimateManager _rebuildTimeEstimator;

  /**
   * Stores (numParticlesBuffer, remainderTraversalTime) pairs for each iteration.
   *
   * Used to predict the remainder traversal time for non-zero buffer sizes,
   * assuming a linear relationship between the particle buffer size (predictor)
   * and the remainder traversal time (response).
   *
   * Has a collector for remainder traversal times for iterations where the particle buffer
   * size is zero. Zero buffer sizes distort linear regression and are therefore
   * handled separately using a mean estimator, because the sum of all times since the last rebuild are still needed.
   */
  RemainderTraversalEstimateManager _remainderTraversalTimePredictor;

  /**
   * Stores the remainder traversal time of the previous iteration.
   * Used by the first decision criterion to estimate the time increase
   * if no rebuild is performed.
   */
  long _lastRemainderTraversalTime{0};

  /// Number of particles in the particle buffer after the fast particles have been added.
  size_t _afterFastAddNumParticlesBuffer{0};

  /**
   * Stores the number of particles added to the buffer during a iteration,
   * due to migrating particles, for estimating the particle buffer filling for the next iteration.
   */
  size_t _lastIncreaseNumParticlesBuffer{0};

  /**
   * Estimated number of particles for the particle buffer based on the last increases
   * and the number of particles in the buffer after fast particles have been added.
   */
  size_t _numParticlesBufferEstimate{0};

  /**
   * Stores the estimated rebuild time for logging.
   */
  double _rebuildNeighborTimeEstimate{std::numeric_limits<double>::quiet_NaN()};

  /**
   * Stores the estimated remainder traversal time for logging.
   */
  double _remainderTraversalTimeEstimate{0.};
};
}  // namespace autopas::utils::DynamicRebuildDecisionWithFPB
