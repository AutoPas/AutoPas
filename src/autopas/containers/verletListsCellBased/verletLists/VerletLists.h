/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "VerletListHelpers.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLListIterationTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/traversals/VLTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * Verlet Lists container.
 * The VerletLists class uses neighborhood lists to calculate pairwise or triwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of the interaction.
 * Cells are created using a cell size of at least cutoff + skin radius.
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletLists : public VerletListsLinkedBase<Particle_T> {
  /**
   * Type of the Particle.
   */
  using ParticleType = Particle_T;
  /**
   * Type of the ParticleCell used by the underlying linked cells.
   */
  using ParticleCellType = FullParticleCell<Particle_T>;

 public:
  /**
   * Enum that specifies how the verlet lists should be built
   */
  enum BuildVerletListType {
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
   * Constructor of the VerletLists class.
   * The neighbor lists are build using a search radius of cutoff + skin.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius per timestep.
   * @param buildVerletListType Specifies how the verlet list should be built, see BuildVerletListType
   * @param cellSizeFactor cell size factor relative to cutoff
   */
  VerletLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
              const double skin, const BuildVerletListType buildVerletListType = BuildVerletListType::VerletSoA,
              const double cellSizeFactor = 1.0)
      : VerletListsLinkedBase<Particle_T>(boxMin, boxMax, cutoff, skin, cellSizeFactor),
        _buildVerletListType(buildVerletListType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletLists; }

  void computeInteractions(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *verletTraversalInterface = dynamic_cast<VLTraversalInterface<ParticleCellType> *>(traversal);
    if (verletTraversalInterface) {
      verletTraversalInterface->setCellsAndNeighborLists(this->_linkedCells.getCells(), _crsNeighborList,
                                                         _indexToParticle, _soaNeighborLists, _aosNeighborPairsLists);
    } else {
      utils::ExceptionHandler::exception(
          "VerletLists::computeInteractions(): Trying to use a traversal of wrong type.");
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  /**
   * get the actual neighbor list
   * @return the neighbor list
   */
  typename VerletListHelpers<Particle_T>::CRSNeighborList &getVerletListsAoS() { return _crsNeighborList; }

  /**
   * Build the pair neighbor list if necessary without fully rebuilding the other neighbor lists.
   * @param traversal
   */
  void prepareForTraversal(TraversalInterface *traversal) override {
    if (traversal->getTraversalType() == TraversalOption::vl_pair_list_iteration) {
      if (not _pairListIsValid) {
        this->updatePairVerletListsAoS3B(this->_verletBuiltNewton3);
      }
    }
  }

  /**
   * Rebuilds the verlet lists, marks them valid and resets the internal counter.
   * @note This function will be called in computeInteractions()
   * @param traversal
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    _soaListIsValid = false;
    _pairListIsValid = false;
    const bool buildWithN3 = traversal->getUseNewton3();
    this->_verletBuiltNewton3 = buildWithN3;

    // Check for triwise traversals
    switch (traversal->getTraversalType()) {
      // Standard pairwise traversal
      case TraversalOption::vl_list_iteration: {
        this->updateVerletListsCRS<InteractionTypeOption::pairwise>(buildWithN3);
        break;
      }
      case TraversalOption::vl_list_intersection: {
        this->updateVerletListsCRS<InteractionTypeOption::triwise>(buildWithN3);
        sortCRSNeighborListsByParticlePtr();

        break;
      }
      case TraversalOption::vl_pair_list_iteration: {
        // build 3Body Verlet lists through VLIteration traversal
        this->updatePairVerletListsAoS3B(buildWithN3);
        break;
      }
      // Default builds normal neighbor lists including halo particles.
      default: {
        this->updateVerletListsCRS<InteractionTypeOption::triwise>(buildWithN3);
      }
    }

    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);

    if (traversal->getDataLayout() == DataLayoutOption::soa) {
      // only do this if we need it, i.e., if we are using soa!
      generateSoAListFromCRSNeighborLists();
    }
  }

 protected:
  void sortCRSNeighborListsByParticlePtr() {
    auto &offsets = _crsNeighborList.offsets();
    auto &neighbors = _crsNeighborList.neighbors();

    const auto indexLessByParticlePtr = [this](size_t lhs, size_t rhs) {
      return _indexToParticle[lhs] < _indexToParticle[rhs];
    };

    AUTOPAS_OPENMP(parallel for schedule(dynamic))
    for (size_t particleIndex = 0; particleIndex < _crsNeighborList.size(); ++particleIndex) {
      std::sort(neighbors.begin() + offsets[particleIndex], neighbors.begin() + offsets[particleIndex + 1],
                indexLessByParticlePtr);
    }
  }

  template <InteractionTypeOption::Value interactionType>
  void updateVerletListsCRS(bool useNewton3) {
    const auto numParticles = generateParticleIndexMap();

    _crsNeighborList.resizeParticles(numParticles);

    const double interactionLength = this->getInteractionLength();
    constexpr auto dataLayout = DataLayoutOption::aos;
    constexpr bool traverseHaloCells = (interactionType == InteractionTypeOption::triwise);

    {
      typename VerletListHelpers<Particle_T>::CRSNeighborCounterFunctor counter(
          _particlePtr2indexMap, _crsNeighborList.offsets(), interactionLength);

      auto traversal =
          LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::CRSNeighborCounterFunctor,
                         traverseHaloCells>(this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), counter,
                                            interactionLength, this->_linkedCells.getCellBlock().getCellLength(),
                                            dataLayout, useNewton3);
      this->_linkedCells.computeInteractions(&traversal);
    }

    std::inclusive_scan(_crsNeighborList.offsets().begin(), _crsNeighborList.offsets().end(),
                        _crsNeighborList.offsets().begin());

    _crsNeighborList.neighbors().resize(_crsNeighborList.offsets().back());

    auto writeOffsets = _crsNeighborList.offsets();

    {
      typename VerletListHelpers<Particle_T>::CRSNeighborFillFunctor filler(
          _particlePtr2indexMap, writeOffsets, _crsNeighborList.neighbors(), interactionLength);

      auto traversal = LCC08Traversal<ParticleCellType, typename VerletListHelpers<Particle_T>::CRSNeighborFillFunctor,
                                      traverseHaloCells>(
          this->_linkedCells.getCellBlock().getCellsPerDimensionWithHalo(), filler, interactionLength,
          this->_linkedCells.getCellBlock().getCellLength(), dataLayout, useNewton3);
      this->_linkedCells.computeInteractions(&traversal);
    }

#ifndef NDEBUG
    for (std::size_t i = 0; i < _crsNeighborList.size(); ++i) {
      if (writeOffsets[i] != _crsNeighborList.offsets()[i + 1]) {
        utils::ExceptionHandler::exception("CRS Verlet list fill mismatch for particle {}", i);
      }
    }
#endif
  }

  /**
   * Update the pair verlet lists for AoS usage
   * @param useNewton3
   */
  void updatePairVerletListsAoS3B(bool useNewton3) {
    updateVerletListsCRS<InteractionTypeOption::triwise>(false);

    const size_t numParticles = _indexToParticle.size();
    std::vector<std::vector<std::pair<Particle_T *, Particle_T *>>> tempPairLists(numParticles);

    const double interactionLength = this->getInteractionLength();
    typename VerletListHelpers<Particle_T>::PairVerletListGeneratorFunctor f(tempPairLists, _particlePtr2indexMap,
                                                                             interactionLength);

    DataLayoutOption dataLayout = DataLayoutOption::aos;
    auto traversal = VLListIterationTraversal<ParticleCellType,
                                              typename VerletListHelpers<Particle_T>::PairVerletListGeneratorFunctor>(
        f, dataLayout, useNewton3);
    this->computeInteractions(&traversal);

    // Flatten tempPairLists into _aosNeighborPairsLists
    _aosNeighborPairsLists.resizeParticles(numParticles);
    auto &offsets = _aosNeighborPairsLists.offsets();
    for (size_t i = 0; i < numParticles; ++i) {
      offsets[i + 1] = tempPairLists[i].size();
    }
    std::inclusive_scan(offsets.begin(), offsets.end(), offsets.begin());

    auto &neighborPairs = _aosNeighborPairsLists.neighborPairs();
    neighborPairs.resize(offsets.back());

        AUTOPAS_OPENMP(parallel for schedule(static))
        for (size_t i = 0; i < numParticles; ++i) {
          std::copy(tempPairLists[i].begin(), tempPairLists[i].end(), neighborPairs.begin() + offsets[i]);
        }

        _pairListIsValid = true;
  }

  /**
   * Clears and then generates the AoS neighbor pairs lists.
   */
  void generateAoSNeighborPairsLists() {
    _aosNeighborPairsLists.clear();
    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      // create the pair verlet list entries for all particles
      _aosNeighborPairsLists[&(*iter)];
    }
  }

  void generateSoAListFromCRSNeighborLists() {
    const size_t numParticles = _crsNeighborList.size();

    _soaNeighborLists.clear();
    _soaNeighborLists.resize(numParticles);

    size_t accumulatedListSize = 0;

    for (size_t particleIndex = 0; particleIndex < numParticles; ++particleIndex) {
      const auto crsNeighbors = _crsNeighborList.neighborsOf(particleIndex);

      accumulatedListSize += crsNeighbors.size();

      auto &soaList = _soaNeighborLists[particleIndex];
      soaList.resize(crsNeighbors.size());

      std::copy(crsNeighbors.begin(), crsNeighbors.end(), soaList.begin());
    }

    AutoPasLog(DEBUG,
               "VerletLists::generateSoAListFromCRSNeighborLists: average verlet list "
               "size is {}",
               numParticles == 0 ? 0.0 : static_cast<double>(accumulatedListSize) / numParticles);

    _soaListIsValid = true;
  }

  size_t generateParticleIndexMap() {
    _indexToParticle.clear();

    const auto estimatedParticles = this->getNumberOfParticles(IteratorBehavior::ownedOrHaloOrDummy);
    _indexToParticle.reserve(estimatedParticles);

    for (auto iter = this->begin(IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
      auto *particlePtr = &(*iter);
      _indexToParticle.emplace_back(particlePtr);
    }

    _particlePtr2indexMap.build(_indexToParticle);

    return _indexToParticle.size();
  }

 private:
  /**
   * Neighbor Pairs Lists: Map of particle pointers to vector of pairs of particle pointers. (To find triplets.)
   */
  VerletListHelpers<Particle_T>::CRSPairNeighborList _aosNeighborPairsLists;

  /**
   * Mapping of every particle, represented by its pointer, to an index.
   * The index indexes all particles in the container.
   */
  PointerToIndexMap _particlePtr2indexMap;

  std::vector<Particle_T *> _indexToParticle;
  VerletListHelpers<Particle_T>::CRSNeighborList _crsNeighborList;

  /**
   * verlet list for SoA:
   * For every Particle, identified via the _particlePtr2indexMap, a vector of its neighbor indices is stored.
   */
  std::vector<std::vector<size_t, AlignedAllocator<size_t>>> _soaNeighborLists;

  std::vector<std::vector<size_t>> _tempNeighbors;

  /**
   * Shows if the SoA neighbor list is currently valid.
   */
  bool _soaListIsValid{false};

  /**
   * Shows if the pair list for triwise interactions is currently valid.
   */
  bool _pairListIsValid{false};

  /**
   * Specifies for what data layout the verlet lists are build.
   */
  BuildVerletListType _buildVerletListType;
};

}  // namespace autopas
