/**
 * @file VLCNeighborListInterface.h
 * @author tirgendetwas
 * @date 27.10.20
 */

#pragma once

#include <numeric>
#include <vector>

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/SoA.h"

namespace autopas {
/**
 * Interface of neighbor lists to be used with VerletListsCells container.
 * @tparam Particle_T Type of particle to be used for the neighbor list.
 */
template <class Particle_T>
class VLCNeighborListInterface {
 public:
  /**
   * Virtual default destructor.
   */
  virtual ~VLCNeighborListInterface() = default;

  /**
   * Builds AoS neighbor list from underlying linked cells object.
   * @param vlcTraversalOpt Traversal for which this neighbor list is built.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   */
  virtual void buildAoSNeighborList(TraversalOption vlcTraversalOpt, LinkedCells<Particle_T> &linkedCells,
                                    bool useNewton3) = 0;

  /**
   * Gets the number of neighbors over all neighbor lists that belong to this particle.
   * @param particle
   * @return the size of the neighbor list(s) of this particle
   */
  virtual size_t getNumberOfPartners(const Particle_T *particle) const = 0;

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
  virtual void generateSoAFromAoS(LinkedCells<Particle_T> &linkedCells) = 0;

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
  void setLinkedCellsPointer(LinkedCells<Particle_T> *linkedCells) { this->_internalLinkedCells = linkedCells; }

 protected:
  /**
   * Internal linked cells structure. Necessary for loading and extracting SoA.
   */
  LinkedCells<Particle_T> *_internalLinkedCells;

  /**
   * Structure of arrays necessary for SoA data layout.
   */
  SoA<typename Particle_T::SoAArraysType> _soa;
};

}  // namespace autopas
