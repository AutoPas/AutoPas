#pragma once

#include <array>
#include <vector>

#include "autopas/utils/Timer.h"
#include "autopas/containers/P3M/P3M_shortRangeFunctor.h"

namespace autopas {

template <class Particle_Type>
class P3M_container;

template <class ParticleCell>
class P3MTraversalInterface {

  using ParticleType = typename ParticleCell::ParticleType;

 public:
  /**
   * Destructor
   */
  virtual ~P3MTraversalInterface() = default;

  
  virtual void set_p3m_traversal_parameters(unsigned int cao, std::array<unsigned int, 3> grid_dims, std::array<double, 3> grid_dist, const std::array<double, 3> &boxMin,
         std::vector<std::vector<double>> &selfForceCoeffs, P3M_container<ParticleType> *container, LCC08Traversal<ParticleCell, P3M_shortRangeFunctor<ParticleType>> *shortRangeTraversal) = 0;

  virtual void set_Timers(utils::Timer *fftTimer, utils::Timer *shortRangeTimer, utils::Timer *chargeAssignmentTimer, utils::Timer *forceInterpolationTimer) = 0;

  };
}