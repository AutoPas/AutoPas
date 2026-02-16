/**
 * @file PseudoVerletLists.h
 * @date 04.12.2025
 * @author Lars Doll
 */

#pragma once
#include <vector>

#include "autopas/cells/SortedCellView.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "traversals/PsVLTraversalInterface.h"

namespace autopas {

/**
 * PseudoVerletLists class. Pseudo Verlet lists store 13 SortedCellViews for each cell,
 * each of which projects the particles onto the cell connection axis to a neighboring cell.
 * By sorting the particles, the loop can be terminated earlier without needing to check
 * each particle to see if it is within the cutoff.
 * By introducing a skin, the sorted cell views do not have to be updated at each calculation step.
 * Based on: https://doi.org/10.1080/08927022.2012.762097
 * Currently only AoS and cell size factors >= 1 are supported
 * @tparam Particle_T type of the Particle
 */
template <class Particle_T>
class PseudoVerletLists : public VerletListsLinkedBase<Particle_T> {
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
   * Constructor of the PseudoVerletLists class.
   * \image html PseudoVerletLists.png "Projection and Sorting of the Particles in 2D"
   * The orientationList is built using sortedCellView.
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param cellSizeFactor cell size factor relative to cutoff.
   */
  PseudoVerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
                    const double skin, const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor) {
    for (size_t i = 0; i < _sortingDirections.size(); ++i) {
      _sortingDirections.emplace_back(utils::ArrayMath::normalize(getDirectionFromIndex(i)));
    }
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::pseudoVerletLists; }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    if (auto *pseudoVerletTraversalInterface = dynamic_cast<PsVLTraversalInterface<ParticleCellType> *>(traversal)) {
      pseudoVerletTraversalInterface->setOrientationList(this->_orientationList);
    } else {
      utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in PseudoVerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
   * Getter.
   * @return _orientationList
   */
  std::vector<std::vector<SortedCellView<ParticleCellType>>> &getOrientationList() { return _orientationList; }

  /**
   * Rebuilds the orientationList.
   * \image html DirectionsForPseudoVerletLists.png "Indices of the directions for the orientationList"
   * @param traversal from ParticleContainerInterface, but is not used here.
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    auto numCells = this->_linkedCells.getCells().size();
    _orientationList.clear();
    _orientationList.resize(numCells);

    for (size_t i = 0; i < numCells; ++i) {
      for (auto &_direction : _sortingDirections) {
        _orientationList[i].emplace_back(SortedCellView{this->_linkedCells.getCells()[i], _direction});
      }
    }
  }

  /**
   * Getter.
   * @return _cells from linkedCells
   */
  std::vector<ParticleCellType> getCells() { return this->_linkedCells.getCells(); }

  /**
   * Converts a 1D index to its corresponding direction.
   * @param index
   * @return unnormalized direction
   */
  static std::array<double, 3> getDirectionFromIndex(size_t index) {
    constexpr std::array<int, 3> dims{3, 3, 3};
    const auto vec = autopas::utils::ThreeDimensionalMapping::oneToThreeD(static_cast<int>(index) + 14, dims);

    return {static_cast<double>(vec[0] - 1), static_cast<double>(vec[1] - 1), static_cast<double>(vec[2] - 1)};
  }

  /**
   * Getter.
   * @return _sortingDirections
   */
  [[nodiscard]] std::vector<std::array<double, 3>> getDirections() const { return _sortingDirections; }

 protected:
  /**
   * Orientation List: or each cell, for each direction, a sortedCellView is stored, sorted along a
   * connecting vector between this cell center and the center of a neighboring cell.
   * Only the sortedCellView with directions to neighboring cells with a greater Index than the base cell are stored.
   */
  std::vector<std::vector<SortedCellView<ParticleCellType>>> _orientationList;

  /**
   * Stores the normalized directions to the neighboring cells. With an index greater than the base cell.
   */
  std::vector<std::array<double, 3>> _sortingDirections;
};

}  // namespace autopas