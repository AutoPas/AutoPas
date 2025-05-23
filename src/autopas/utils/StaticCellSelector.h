/**
 * @file StaticCellSelector.h
 * @author seckler
 * @date 02.05.19
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/linkedCells/LinkedCellsReferences.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"

namespace autopas::utils {


template<class Particle_T, class Functor>
std::unique_ptr<TraversalInterface> generateTraversalFromConfig(const Configuration& config, Functor &functor,
                                                                 const TraversalSelectorInfo &traversalInfo) {
  switch (config.container) {
    case ContainerOption::Value::linkedCellsReferences:
      return TraversalSelector<ReferenceParticleCell<Particle_T>>::template generateTraversal<Functor>(config.traversal, functor, traversalInfo, config.dataLayout, config.newton3);
    default:
      return TraversalSelector<FullParticleCell<Particle_T>>::template generateTraversal<Functor>(config.traversal, functor, traversalInfo, config.dataLayout, config.newton3);
  }
}
}  // namespace autopas::utils