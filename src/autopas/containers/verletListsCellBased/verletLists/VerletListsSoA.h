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
    if (mortonTraversalInterface && this->_orderCellsByMortonIndex) {
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
    selectVLGenerationFunctionAtRuntime();
    typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA f(this->_soaNeighborLists, this->getCutoff() + this->getVerletSkin());

    DataLayoutOption dataLayout;
    dataLayout = DataLayoutOption::soa;

    auto traversal =
        LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::VerletListGeneratorFunctorSoA>(
            this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), &f, this->getInteractionLength(),
            this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);

    auto *mortonTraversalInterface = dynamic_cast<MortonIndexTraversalInterface *> (&traversal);
    if (mortonTraversalInterface && this->_orderCellsByMortonIndex) {
      mortonTraversalInterface->setCellsByMortonIndex(this->_linkedCells.getCellsByMortonIndex());
    }

    this->_linkedCells.computeInteractions(&traversal);

    if (this->_sortVerletLists) {
      #pragma omp parallel for schedule(static)
      for (auto &nl : _soaNeighborLists) {
        std::ranges::sort(nl);
      }
    }

    this->_soaListIsValid = true;
  }

  size_t selectVLGenerationFunctionAtRuntime() {
    const uint8_t key =
        (static_cast<uint8_t>(this->_bucketSortParticles)      << 0) |
        (static_cast<uint8_t>(this->_orderCellsByMortonIndex)  << 1) |
        (static_cast<uint8_t>(this->_reserveVLSizes)           << 2);

    switch (key) {
      case 0: return generateVerletListsSoA_Options<false,false,false>();
      case 1: return generateVerletListsSoA_Options< true,false,false>();
      case 2: return generateVerletListsSoA_Options<false, true,false>();
      case 3: return generateVerletListsSoA_Options< true, true,false>();
      case 4: return generateVerletListsSoA_Options<false,false, true>();
      case 5: return generateVerletListsSoA_Options< true,false, true>();
      case 6: return generateVerletListsSoA_Options<false, true, true>();
      case 7: return generateVerletListsSoA_Options< true, true, true>();
      default: __builtin_unreachable();
    }
  }

  template<bool bucketSortParticles, bool orderCellsByMortonIndex, bool reserveVLSizes>
  size_t generateVerletListsSoA_Options() {
    size_t numParticles = 0;
    size_t index = 0;
    std::vector<size_t> oldVerletSizes;
    oldVerletSizes.reserve(_oldNumParticles);

    auto &cells = this->_linkedCells.getCells();
    #pragma omp parallel for schedule(static)
    for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
      auto &cellParticles = cells[cellId]._particles;
      const size_t particleNumCell = cellParticles.size();

      if constexpr (bucketSortParticles) {
        if (_sortingCounter == this->_sortingFrequency && cellParticles.size() > 1) {
          const auto [cellLowerCorner, cellUpperCorner] = this->_linkedCells.getCellBlock().getCellBoundingBox(cellId);
          const double midX = 0.5 * (cellLowerCorner[0] + cellUpperCorner[0]);
          const double midY = 0.5 * (cellLowerCorner[1] + cellUpperCorner[1]);
          const double midZ = 0.5 * (cellLowerCorner[2] + cellUpperCorner[2]);

          std::array<size_t, 8> particlesPerKey{};
          for (const auto &particle : cellParticles) {
            const auto &pos= particle.getR();
            const uint8_t key = static_cast<uint8_t>((pos[0] >= midX) | ((pos[1] >= midY) << 1) | ((pos[2] >= midZ) << 2));
            ++particlesPerKey[key];
          }

          std::array<std::vector<Particle_T>, 8> buckets{};
          for (uint8_t i = 0; i < 8; ++i) {
            buckets[i].clear();
            buckets[i].reserve(particlesPerKey[i]);
          }

          for (auto &particle : cellParticles) {
            const auto &pos = particle.getR();
            const uint8_t key = static_cast<uint8_t>((pos[0] >= midX) | ((pos[1] >= midY) << 1) | ((pos[2] >= midZ) << 2));
            buckets[key].push_back(std::move(particle));
          }

          cellParticles.clear();
          cellParticles.reserve(particleNumCell);
          for (uint8_t i = 0; i < 8; i++) {
            cellParticles.insert(cellParticles.end(),
              std::make_move_iterator(buckets[i].begin()),
              std::make_move_iterator(buckets[i].end()));
          }
        }
      }
    }

    if constexpr (orderCellsByMortonIndex) {
      const auto cellsByMortonIndex = this->_linkedCells.getCellsByMortonIndex();
      for (size_t cellId : cellsByMortonIndex) {
        cells[cellId]._particleSoABuffer.setParticlesIndexInSoAStart(index);
        for (size_t i = 0; i < cells[cellId]._particles.size();  ++i, ++numParticles, ++index) {
          Particle_T &particleI = cells[cellId]._particles[i];
          if constexpr(reserveVLSizes) {
            oldVerletSizes.push_back(particleI.getIndexInSoA() != std::numeric_limits<size_t>::max() && particleI.getIndexInSoA() < _soaNeighborLists.size() ?
                           _soaNeighborLists[particleI.getIndexInSoA()].size(): 64);
          }
          particleI.setIndexInSoA(index);
        }
      }
    } else {
      for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
        cells[cellId]._particleSoABuffer.setParticlesIndexInSoAStart(index);
        for (size_t i = 0; i < cells[cellId]._particles.size();  ++i, ++numParticles, ++index) {
          Particle_T &particleI = cells[cellId]._particles[i];
          if constexpr(reserveVLSizes) {
            oldVerletSizes.push_back(particleI.getIndexInSoA() != std::numeric_limits<size_t>::max() && particleI.getIndexInSoA() < _soaNeighborLists.size() ?
                           _soaNeighborLists[particleI.getIndexInSoA()].size(): 64);
          }
          particleI.setIndexInSoA(index);
        }
      }
    }

    _oldNumParticles = numParticles;

    this->_soaNeighborLists.clear();
    this->_soaNeighborLists.resize(numParticles);

    index = 0;

    for (size_t particleI = 0; particleI < numParticles; ++particleI) {
      _soaNeighborLists[particleI].clear();
      if constexpr(reserveVLSizes) {
        _soaNeighborLists[particleI].reserve(oldVerletSizes[particleI]);
      }
    }

    if constexpr (bucketSortParticles) {
      if (_sortingCounter == this->_sortingFrequency) {
        _sortingCounter = 1;
      } else {
        _sortingCounter++;
      }
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

  size_t _sortingCounter = this->_sortingFrequency;
};

};// namespace autopas
