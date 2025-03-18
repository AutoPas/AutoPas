/**
 * @file VLCNeighborListInterface.h
 * @author tirgendetwas
 * @date 27.10.20
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {
/**
 * Interface of neighbor lists to be used with VerletListsCells container.
 * @tparam ParticleT Type of particle to be used for the neighbor list.
 */
template <class ParticleT>
class VLCNeighborListInterface {
 public:
  /**
   * Virtual default destructor.
   */
  virtual ~VLCNeighborListInterface() = default;

  /**
   * Builds AoS neighbor list from underlying linked cells object.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff Cutoff radius.
   * @param skin Skin of the verlet list.
   * @param interactionLength Interaction length of the underlying linked cells object.
   * @param vlcTraversalOpt Traversal for which this neighbor list is built.
   * @param buildType Type of build functor to be used for the generation of the neighbor list.
   */
  virtual void buildAoSNeighborList(LinkedCells<ParticleT> &linkedCells, bool useNewton3, double cutoff, double skin,
                                    double interactionLength, const TraversalOption vlcTraversalOpt,
                                    typename VerletListsCellsHelpers::VLCBuildType buildType) = 0;

  /**
   * Gets the number of neighbors over all neighbor lists that belong to this particle.
   * @param particle
   * @return the size of the neighbor list(s) of this particle
   */
  virtual size_t getNumberOfPartners(const ParticleT *particle) const = 0;

  /**
   * Returns the container type of this neighbor list and the container it belongs to.
   * @return ContainerOption for this neighbor list and the container it belongs to.
   */
  [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

  /**
   * Generates neighbor list in SoA layout from available neighbor list in AoS layout.
   * Copies the structure of the AoS neighbor list and replaces the particle pointers with the global indices of the
   * particles.
   * @param linkedCells Underlying linked cells structure.
   */
  virtual void generateSoAFromAoS(LinkedCells<ParticleT> &linkedCells) = 0;

  /**
   * Loads cells into structure of arrays.
   * @tparam TFunctor
   * @param f Functor that handles the loading.
   * @return loaded structure of arrays
   */
  template <class TFunctor>
  auto *loadSoA(TFunctor *f) {
    _soa.clear();

    // First resize the SoA to the required number of elements to store. This avoids resizing successively the SoA in
    // SoALoader.
    auto &cells = _internalLinkedCells->getCells();
    std::vector<size_t> offsets(cells.size() + 1);
    std::inclusive_scan(
        cells.begin(), cells.end(), offsets.begin() + 1,
        [](const size_t &partialSum, const auto &cell) { return partialSum + cell.size(); }, 0);
    _soa.resizeArrays(offsets.back());

    AUTOPAS_OPENMP(parallel for)
    for (size_t i = 0; i < cells.size(); ++i) {
      f->SoALoader(cells[i], _soa, offsets[i], /*skipSoAResize*/ true);
    }

    return &_soa;
  }

  /**
   * Extracts cells from structure of arrays.
   * @tparam TFunctor
   * @param f Functor that handles the extraction.
   */
  template <class TFunctor>
  void extractSoA(TFunctor *f) {
    size_t offset = 0;
    for (auto &cell : _internalLinkedCells->getCells()) {
      f->SoAExtractor(cell, _soa, offset);
      offset += cell.size();
    }
  }

  /**
   * Assigns the current traversal to the correct traversal interface. The choice of traversal interface depends on the
   * type(s) of neighbor list allowed to use the current traversal. Currently VLCCellPairTraversalInterface handles the
   * traversals allowed solely for VLCCellPairNeighborList, while VLCTraversalInterface handles the traversals allowed
   * for both VLCCellPairNeighborList and VLCAllCellsNeighborList.
   * @param traversal the current traversal
   */
  virtual void setUpTraversal(TraversalInterface *traversal) = 0;

  /**
   * Set the Linked Cells Pointer for this List.
   * @param linkedCells
   */
  void setLinkedCellsPointer(LinkedCells<ParticleT> *linkedCells) { this->_internalLinkedCells = linkedCells; }

 protected:
  /**
   * Internal linked cells structure. Necessary for loading and extracting SoA.
   */
  LinkedCells<ParticleT> *_internalLinkedCells;

  /**
   * Structure of arrays necessary for SoA data layout.
   */
  SoA<typename ParticleT::SoAArraysType> _soa;

  /**
   * Creates and applies a generator functor for the building of the neighbor list.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff Cutoff radius.
   * @param skin Skin of the verlet list.
   * @param interactionLength Interaction length of the underlying linked cells object.
   * @param vlcTraversalOpt Traversal for which this neighbor list is built.
   * @param buildType Type of build functor to be used for the generation of the neighbor list.
   */
  virtual void applyBuildFunctor(LinkedCells<ParticleT> &linkedCells, bool useNewton3, double cutoff, double skin,
                                 double interactionLength, const TraversalOption &vlcTraversalOpt,
                                 typename VerletListsCellsHelpers::VLCBuildType buildType) = 0;
};

}  // namespace autopas
