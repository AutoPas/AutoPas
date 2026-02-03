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
      if (this->_preloadLJMixingPtr) {
        verletTraversalInterface->setPreloadMixingLJPtr(true);
      }
    } else {
      utils::ExceptionHandler::exception("trying to use a traversal of wrong type in VerletLists::computeInteractions");
    }

    auto *mortonTraversalInterface = dynamic_cast<MortonIndexTraversalInterface *> (traversal);
    if (mortonTraversalInterface && this->_useMortonIndex) {
      mortonTraversalInterface->setCellsByMortonIndex(this->_linkedCells.getCellsByMortonIndex());
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
    if (this->_useMortonIndex && this->_useLiveId && this->_reserveVLSizes) {
      generateSoANeighborListsWithMortonOrderAndReserving();
    } else if (this->_useMortonIndex && this->_useLiveId && !this->_reserveVLSizes) {
      generateSoANeighborListsWithMortonOrder();
    } else if (this->_useLiveId && this->_reserveVLSizes && !this->_useMortonIndex) {
      generateSoANeighborListsWithReserving();
    } else {
      generateSoANeighborLists();
    }
    typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA f(this->_soaNeighborLists, this->getCutoff() + this->getVerletSkin());

    DataLayoutOption dataLayout;
    dataLayout = DataLayoutOption::soa;

    auto traversal =
        LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
            this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);

    auto *mortonTraversalInterface = dynamic_cast<MortonIndexTraversalInterface *> (&traversal);
    if (mortonTraversalInterface && this->_useMortonIndex) {
      mortonTraversalInterface->setCellsByMortonIndex(this->_linkedCells.getCellsByMortonIndex());
    }

    this->_linkedCells.computeInteractions(&traversal);

    this->_soaListIsValid = true;
  }

  uint32_t generateSoANeighborLists() {
    uint32_t numParticles = 0;
    uint32_t index = 0;

    // DON'T simply parallelize this loop!!! this needs modifications if you want to parallelize it!
    // We have to iterate also over dummy particles here to ensure a correct size of the arrays.
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles, ++index) {

      iter->setLiveId(index);
    }

    this->_soaNeighborLists.clear();
    this->_soaNeighborLists.resize(numParticles);

    index = 0;

    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++index) {
      this->_soaNeighborLists[index].clear();
    }

    return numParticles;
  }

  /**
   * @return Number of particles in the container
   */
   uint32_t generateSoANeighborListsWithReserving() {
    uint32_t numParticles = 0;
    uint32_t index = 0;
    std::vector<uint16_t> oldVerletSizes;
    oldVerletSizes.reserve(_oldNumParticles);

     // DON'T simply parallelize this loop!!! this needs modifications if you want to parallelize it!
     // We have to iterate also over dummy particles here to ensure a correct size of the arrays.
     for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++numParticles, ++index) {

       oldVerletSizes.push_back(!_soaNeighborLists.empty() && !_soaNeighborLists[index].empty() ? _soaNeighborLists[index].size() : 64);
       iter->setLiveId(index);
     }

    _oldNumParticles = numParticles;

     this->_soaNeighborLists.clear();
     this->_soaNeighborLists.resize(numParticles);

     index = 0;

     for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter, ++index) {
       this->_soaNeighborLists[index].clear();
       this->_soaNeighborLists[index].reserve(oldVerletSizes[index]);
     }

     return numParticles;
   }

  uint32_t generateSoANeighborListsWithMortonOrder() {
    uint32_t numParticles = 0;
    uint32_t index = 0;

    auto &cells = this->_linkedCells.getCells();
    const auto cellsByMortonIndex = this->_linkedCells.getCellsByMortonIndex();

    for (size_t cellId : cellsByMortonIndex) {
      for (uint32_t i = 0; i < cells[cellId]._particles.size();  ++i, ++numParticles, ++index) {
        Particle_T &particleI = cells[cellId]._particles[i];

        particleI.setLiveId(index);
      }
    }


    this->_soaNeighborLists.clear();
    this->_soaNeighborLists.resize(numParticles);

    index = 0;

    for (uint32_t particleI = 0; particleI < numParticles; ++particleI) {
      _soaNeighborLists[particleI].clear();
    }

    return numParticles;
  }

  uint32_t generateSoANeighborListsWithMortonOrderAndReserving() {
    uint32_t numParticles = 0;
    uint32_t index = 0;
    std::vector<uint16_t> oldVerletSizes;
    oldVerletSizes.reserve(_oldNumParticles);

    auto &cells = this->_linkedCells.getCells();
    const auto cellsByMortonIndex = this->_linkedCells.getCellsByMortonIndex();

    for (size_t cellId : cellsByMortonIndex) {
      for (uint32_t i = 0; i < cells[cellId]._particles.size();  ++i, ++numParticles, ++index) {
        Particle_T &particleI = cells[cellId]._particles[i];

        oldVerletSizes.push_back( particleI.getLiveId() != std::numeric_limits<size_t>::max() && particleI.getLiveId() < _soaNeighborLists.size() ?
        oldVerletSizes.push_back( particleI.getLiveId() != std::numeric_limits<uint32_t>::max() && particleI.getLiveId() < _soaNeighborLists.size() ?
                       _soaNeighborLists[particleI.getLiveId()].size(): 64);
        particleI.setLiveId(index);
      }
    }

    _oldNumParticles = numParticles;

    // std::cout << "\n";

    this->_soaNeighborLists.clear();
    this->_soaNeighborLists.resize(numParticles);

    index = 0;

    for (uint32_t particleI = 0; particleI < numParticles; ++particleI) {
      _soaNeighborLists[particleI].clear();
      _soaNeighborLists[particleI].reserve(oldVerletSizes[particleI]);
    }

    return numParticles;
  }

private:
  /**
   * Neighbor Lists: Map of particle pointers to vector of particle pointers.
   */
  typename VerletListHelpers<Particle_T>::NeighborListAoSType _aosNeighborLists;

  /**
   * verlet list for SoA:
   * For every Particle, a vector of its neighbor indices is stored.
   */
  std::vector<std::vector<autopas::SoAIndexIntType, AlignedAllocator<autopas::SoAIndexIntType>>> _soaNeighborLists;

  /**
   * Shows if the SoA neighbor list is currently valid.
   */
  bool _soaListIsValid{false};

  /**
   * Specifies for what data layout the verlet lists are build.
   */
  BuildVerletListType _buildVerletListType;

  uint32_t _oldNumParticles = 0;
  size_t _oldNumParticles = 0;
};

};// namespace autopas
