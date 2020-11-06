/**
 * @file StaticSelectorMacros.h
 * @author seckler
 * @date 02.05.19
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/ParticleCell.h"
#include "autopas/cells/ReferenceParticleCell.h"

namespace autopas::utils {

/**
 *
 * @tparam F
 * @param theBool
 * @param func
 * @return the return value of func
 * @todo c++20 change to explicit template for templates.
 * @todo function calls also have to be changed to
 * ```
 * [&]<typename ParticleCellType>() { ... }
 * ```
 * instead of
 * ```
 * [&](auto particleCellDummy) {...}
 * ```
 */
template <typename ParticleType, typename F>
decltype(auto) withStaticCellType(autopas::CellType cellType, F func) {
  switch (cellType) {
    case autopas::CellType::ClusterTower:
      [[fallthrough]];
    case autopas::CellType::SortedCellView:
      [[fallthrough]];
    case autopas::CellType::IsNoCell:
      [[fallthrough]];
    case autopas::CellType::FullParticleCell:
      // todo c++20: return func.template operator()<autopas::FullParticleCell<ParticleType>>();
      return func(autopas::FullParticleCell<ParticleType>());
    case autopas::CellType::ReferenceParticleCell:
      // todo c++20: return func.template operator()<autopas::ReferenceParticleCell<ParticleType>>();
      return func(autopas::ReferenceParticleCell<ParticleType>());
    default:
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of of a Celltype not specified in TravelComparison::calculateForces. "
          "CelltypeEnum: {}",
          cellType);
  }
}
}  // namespace autopas::utils