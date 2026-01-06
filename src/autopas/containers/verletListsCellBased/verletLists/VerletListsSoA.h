#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
//#include <likwid-marker.h>


namespace autopas {

/**
 * Verlet Lists SoA container.
 * The VerletLists class uses neighborhood lists to calculate pairwise
 * interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletListsSoA : public VerletListsLinkedBase<Particle_T> {

  /**
 * Type of the Particle.
 */
  using ParticleType = Particle_T;
  /**
   * Type of the ParticleCell used by the underlying linked cells.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

public:
  enum class BuildVerletListType {
    /**
     * Build it using AoS
     */
    VerletAoS,
    /**
     * Build it using AoS
     */
    VerletSoA,
  };

  /**
   * Constructor of the VerletListsSoA class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per timestep.
   * @param buildVerletListType Specifies how the verlet list should be build, see BuildVerletListType
   * @param cellSizeFactor cell size factor ralative to cutoff
   */
  VerletListsSoA(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

  /**
 * @copydoc ParticleContainerInterface::getContainerType()
 */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletListsSoA; }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<ParticleCellType> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _aosNeighborLists,
                                                         _soaNeighborLists);
    } else {
      utils::ExceptionHandler::exception("trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    // LIKWID_MARKER_START("force calculation");
    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
    // LIKWID_MARKER_STOP("force calculation");
  }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in computeInteractions()!
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();
    this->updateSoAVerletLists(traversal->getUseNewton3());
    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);
  }

 protected:
  /**
   * @param useNewton3
   */
  void updateSoAVerletLists(bool useNewton3) {
    generateSoANeighborLists();
    typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA f(this->_soaNeighborLists, this->_particlePtr2indexMap,
                                                                         this->getCutoff() + this->getVerletSkin());

    DataLayoutOption dataLayout;
    dataLayout = DataLayoutOption::soa;

    auto traversal =
        LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
            this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
    this->_linkedCells.computeInteractions(&traversal);

    this->_soaListIsValid = true;
  }

  /**
   * @return Number of particles in the container
   */
  size_t generateSoANeighborLists() {
    size_t numParticles = 0;
    size_t index = 0;
    this->_soaNeighborLists.clear();
    this->_particlePtr2indexMap.clear();
    // DON'T simply parallelize this loop!!! this needs modifications if you want to parallelize it!
    // We have to iterate also over dummy particles here to ensure a correct size of the arrays.
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles, ++index) {
      this->_particlePtr2indexMap[&(*iter)] = index;
    }
    this->_soaNeighborLists.resize(numParticles);

    for (auto &neighborlist : this->_soaNeighborLists) {
      neighborlist.clear();
      neighborlist.reserve(64);
    }

    return numParticles;
  }

private:
  /**
   * Neighbor Lists: Map of particle pointers to vector of particle pointers.
   */
  typename VerletListHelpers<Particle_T>::NeighborListAoSType _aosNeighborLists;

  /**
   * Mapping of every particle, represented by its pointer, to an index.
   * The index indexes all particles in the container.
   */
  std::unordered_map<const Particle_T *, size_t> _particlePtr2indexMap;

  /**
   * verlet list for SoA:
   * For every Particle, identified via the _particlePtr2indexMap, a vector of its neighbor indices is stored.
   */
  std::vector<std::vector<size_t, AlignedAllocator<size_t>>> _soaNeighborLists;

  /**
   * Shows if the SoA neighbor list is currently valid.
   */
  bool _soaListIsValid{false};

  /**
   * Specifies for what data layout the verlet lists are build.
   */
  BuildVerletListType _buildVerletListType;
};

};// namespace autopas
