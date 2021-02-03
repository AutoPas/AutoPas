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
 * @tparam Particle Type of particle to be used for the neighbor list.
 */
template <class Particle>
class VLCNeighborListInterface {
 public:
  /**
   * Default destructor.
   */
  ~VLCNeighborListInterface() = default;

  /**
   * Builds AoS neighbor list from underlying linked cells object.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff Cutoff radius.
   * @param skin Skin of the verlet list.
   * @param interactionLength Interaction length of the underlying linked cells object.
   * @param buildTraversalOption Traversal option necessary for generator functor.
   * @param buildType Type of build functor to be used for the generation of the neighbor list.
   */
  virtual void buildAoSNeighborList(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                                    double interactionLength, const TraversalOption buildTraversalOption,
                                    typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) = 0;

  /**
   * Gets the number of neighbors over all neighbor lists that belong to this particle.
   * @param particle
   * @return the size of the neighbor list(s) of this particle
   */
  virtual const size_t getNumberOfPartners(const Particle *particle) const = 0;

  /**
   * Returns the container type of this neighbor list and the container it belongs to.
   * @return ContainerOption for this neighbor list and the container it belongs to.
   */
  [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

  /**
   * Generates neighbor list in SoA layout from available neighbor list in AoS layout.
   * @param linkedCells Underlying linked cells structure.
   */
  virtual void generateSoAFromAoS(LinkedCells<Particle> &linkedCells) = 0;

  /**
   * Loads cells into structure of arrays.
   * @tparam TFunctor
   * @param f Functor that handles the loading.
   * @return loaded structure of arrays
   */
  template <class TFunctor>
  auto *loadSoA(TFunctor *f) {
    _soa.clear();
    size_t offset = 0;
    for (auto &cell : _internalLinkedCells->getCells()) {
      f->SoALoader(cell, _soa, offset);
      offset += cell.numParticles();
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
      offset += cell.numParticles();
    }
  }

 protected:
  /**
   * Internal linked cells structure. Necessary for loading and extracting SoA.
   */
  LinkedCells<Particle> *_internalLinkedCells;

  /**
   * Structure of arrays necessary for SoA data layout.
   */
  SoA<typename Particle::SoAArraysType> _soa;

  /**
   * Creates and applies generator functor for the building of the neighbor list.
   * @param linkedCells Linked Cells object used to build the neighbor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff Cutoff radius.
   * @param skin Skin of the verlet list.
   * @param interactionLength Interaction length of the underlying linked cells object.
   * @param buildTraversalOption Traversal option necessary for generator functor.
   * @param buildType Type of build functor to be used for the generation of the neighbor list.
   * */
  virtual void applyBuildFunctor(LinkedCells<Particle> &linkedCells, bool useNewton3, double cutoff, double skin,
                                 double interactionLength, const TraversalOption buildTraversalOption,
                                 typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType) = 0;
};

}  // namespace autopas
