/**
 * @file ContainerOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the container choices.
 */
class ContainerOption : public Option<ContainerOption> {
 public:
  /**
   * Possible choices for the particle container type.
   */
  enum Value {
    /**
     * DirectSum : O(N^2) distance check of all particles and summation of those in cutoff.
     * Minimal overhead but bad complexity.
     */
    directSum,
    /**
     * LinkedCells : Segmentation of the domain into a regular cell grid. Only interactions with particles from
     * neighbor cells are considered. Good data locality and vectorizability but low hit rate of particles in cutoff.
     */
    linkedCells,
    /**
     * LinkedCellsReferences : Same algorithm as LinkedCells but stores all particles in one big vector. Cells only
     * contain references to this vector.
     */
    linkedCellsReferences,
    /**
     * VarVerletLists interface with neighbor list type VerletNeighborListAsBuild : Same algorithm as VerletLists.
     * Remembers which thread created the neighbor list of each particle to exploit this information to avoid data
     * races during the parallel force calculation.
     */
    varVerletListsAsBuild,
    /**
     * VerletClusterLists : Particles are grouped in clusters of fixed size. Similar to VerletLists for every cluster
     * a list of neighbor clusters is generated. Clusters always interact with whole clusters so vectorization is
     * possible.
     */
    verletClusterLists,
    /**
     * VerletLists : Built on top of LinkedCells, a neighbor list is generated for every particle and updated in
     * fixed intervals. Memory access, also in SoA mode is scattered but high hit rate of particles in cutoff.
     */
    verletLists,
    /**
     * VerletListsCells : Similar to VerletLists but Lists are associated with the underlying cells to achieve location
     * information. Parallelization options similar to LinkedCells.
     */
    verletListsCells,
    /**
     * PairwiseVerletLists : Also similar to VerletLists but the lists are associated to each pair of neighboring cells.
     * Improves data locality and cache efficiency.
     */
    pairwiseVerletLists,
  };

  /**
   * Constructor.
   */
  ContainerOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ContainerOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<ContainerOption> getDiscouragedOptions() {
    return {Value::directSum, Value::verletClusterCells, Value::linkedCellsReferences};
  }

  /**
   * Provides a way to iterate over the possible choices of ContainerOption.
   * @return map option -> string representation
   */
  static std::map<ContainerOption, std::string> getOptionNames() {
    return {
        {ContainerOption::directSum, "DirectSum"},
        {ContainerOption::linkedCells, "LinkedCells"},
        {ContainerOption::linkedCellsReferences, "LinkedCellsReferences"},
        {ContainerOption::verletLists, "VerletLists"},
        {ContainerOption::verletListsCells, "VerletListsCells"},
        {ContainerOption::verletClusterLists, "VerletClusterLists"},
        {ContainerOption::varVerletListsAsBuild, "VarVerletListsAsBuild"},
        {ContainerOption::verletClusterCells, "VerletClusterCells"},
        {ContainerOption::pairwiseVerletLists, "PairwiseVerletLists"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
