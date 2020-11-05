/**
 * @file UnknowingCellBorderAndFlagManager.h
 * @author seckler
 * @date 25.03.20
 */

#pragma once
#include "autopas/containers/CellBorderAndFlagManager.h"

namespace autopas::internal {
/**
 * This is a FlagManager that does not know about the actual cells it contains.
 * It assumes that in every cell there can be any type of particle.
 */
class UnknowingCellBorderAndFlagManager : public CellBorderAndFlagManager {
  using index_t = std::size_t;

 public:
  bool cellCanContainHaloParticles(index_t index1d) const override { return true; }
  bool cellCanContainOwnedParticles(index_t index1d) const override { return true; }

  /**
   * Get the static instance of this class.
   * @return one instance.
   */
  static auto& get() {
    const static UnknowingCellBorderAndFlagManager unknowingCellBorderAndFlagManager;
    return unknowingCellBorderAndFlagManager;
  }
};


}  // namespace autopas::internal