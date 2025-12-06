/**
 * @file PseudoVerletLists.h
 *
 * @date 04 Dec 2025
 * @author Lars Doll
 */

#pragma once
#include <vector>

#include "autopas/cells/SortedCellView.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class Particle_T>
class PseudoVerletLists : public LinkedCells<Particle_T> {
  public:
  /**
   * Type of the Particle.
   */
  using ParticleType = Particle_T;

  /**
   * Type of the ParticleCell used by the underlying linked cells.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

  /**
 * Constructor of the PseudoVerletLists class
 * @param boxMin
 * @param boxMax
 * @param cutoff
 * @param skin
 * @param cellSizeFactor cell size factor relative to cutoff
 * @param sortingThreshold number of particles in two cells from which sorting should be performed
 * @param loadEstimator the load estimation algorithm for balanced traversals.
 * By default all applicable traversals are allowed.
 */
  PseudoVerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                    const double skin, const double cellSizeFactor = 1.0, const size_t sortingThreshold = 8,
                    LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : LinkedCells<Particle_T>(boxMin, boxMax, cutoff, skin, std::max(1.0, cellSizeFactor), sortingThreshold, loadEstimator) {
    if (cellSizeFactor < 1.0) {
      AutoPasLog(DEBUG, "VerletListsLinkedBase: CellSizeFactor smaller 1 detected. Set to 1.");
    }
    _orientationList();
    for (size_t i = 0; i < _directions.size(); ++i) {
      _directions[i] = utils::ArrayMath::normalize(_rawDirections[i]);
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    auto numCells = this->_cells.size();
    _orientationList.clear();
    _orientationList.resize(numCells);

    for (size_t i = 0; i < numCells; ++i) {
      for (auto & _direction : _directions) {
        _orientationList[i].emplace_back(SortedCellView{_linkedCells.getCells()[i], _direction});
      }
    }
  }

  protected:

  std::vector<std::vector<SortedCellView<ParticleCellType>>> _orientationList;

  LinkedCells<Particle_T> _linkedCells;

  LoadEstimatorOption _loadEstimator;

  std::array<std::array<double, 3>,13> _directions{};

  static constexpr std::array<std::array<double,3>, 13> _rawDirections = {{
    {{-1,  1, 0}},
    {{ 0,  1, 0}},
    {{ 1,  1, 0}},
    {{ 1,  0, 0}},
    {{-1,  1, 1}},
    {{ 0,  1, 1}},
    {{ 1,  1, 1}},
    {{-1,  0, 1}},
    {{ 0,  0, 1}},
    {{ 1,  0, 1}},
    {{-1, -1, 1}},
    {{ 0, -1, 1}},
    {{ 1, -1, 1}},
}};
};

} // namespace autopas