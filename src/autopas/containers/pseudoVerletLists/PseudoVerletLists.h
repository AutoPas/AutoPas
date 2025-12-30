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
#include "traversals/PsVLTraversalInterface.h"

namespace autopas {

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
     : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor){

    for (size_t i = 0; i < _directions.size(); ++i) {
      _directions[i] = utils::ArrayMath::normalize(_rawDirections[i]);
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
      utils::ExceptionHandler::exception("trying to use a traversal of wrong type in PseudoVerletLists::computeInteractions");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
  * Get the orientationList.
  * @return the orientationList.
  */
  std::vector<std::vector<SortedCellView<ParticleCellType>>> &getOrientationList() {return _orientationList;}

  /**
  * Rebuilds the orientationList.
  * \image html DirectionsForPseudoVerletLists.png "Indices of the directions for the orientationList"
  * @param traversal
  */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    auto numCells = this->_linkedCells.getCells().size();
    _orientationList.clear();
    _orientationList.resize(numCells);

    for (size_t i = 0; i < numCells; ++i) {
      for (auto & _direction : _directions) {
        _orientationList[i].emplace_back(SortedCellView{this->_linkedCells.getCells()[i], _direction});
      }
    }
  }

  std::vector<ParticleCellType> getCells() { return this->_linkedCells.getCells(); }

  protected:
  /**
   * Orientation List: For each cell, 13 sortedCellViews are stored, each of which sorts in the direction of the neighboring cell.
   */
  std::vector<std::vector<SortedCellView<ParticleCellType>>> _orientationList;

  /**
   * Stores the normalized directions to the neighboring cells.
   */
  std::array<std::array<double, 3>,13> _directions{};

  /**
  * Stores the directions to the neighboring cells.
  */
  static constexpr std::array<std::array<double,3>, 13> _rawDirections = {{
    {{1,  0, 0}},
    {{ -1, 1, 0}},
    {{ 0, 1, 0}},
    {{ 1, 1, 0}},
    {{-1, -1, 1}},
    {{ 0, -1, 1}},
    {{ 1, -1, 1}},
    {{-1, 0, 1}},
    {{0, 0, 1}},
    {{ 1, 0, 1}},
    {{-1, 1, 1}},
    {{ 0, 1, 1}},
    {{ 1, 1, 1}},
}};
};

} // namespace autopas