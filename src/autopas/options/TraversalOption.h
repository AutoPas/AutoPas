/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the traversal choices.
 */
class TraversalOption : public Option<TraversalOption> {
 public:
  /**
   * Possible choices for the cell traversal. Traversals marked with '+' can be used for both pairwise and triwise
   * interactions. Try to maintain lexicographic ordering.
   */
  enum Value {
    // DirectSum Traversals:
    /**
     * + DSSequentialTraversal : Sequential nested loop over all particles.
     */
    ds_sequential,

    // LinkedCell Traversals:
    /**
     * + LCC01Traversal : Every cell interacts with all neighbors. Is not compatible with Newton3 thus embarrassingly
     * parallel. Good load balancing and no overhead.
     */
    lc_c01,
    /**
     * LCC01CombinedSoATraversal : Same as LCC01Traversal but SoAs are combined into a circular buffer and the domain
     * is traversed line-wise.
     */
    lc_c01_combined_SoA,
    /**
     * LCC04Traversal : Four-way domain coloring using plus-shaped clusters of cells that are processed with the c08
     * base-step. Less scheduling overhead than LCC08Traversal because of fewer barriers but more coarse-grained.
     */
    lc_c04,
    /**
     * LCC04HCPTraversal : Same as LCC04Traversal but with only one block shape.
     */
    lc_c04_HCP,
    /**
     * LCC04CombinedSoATraversal : Combination of LCC08Traversal and the combined SoA variant of LCC01Traversal.
     * Stripe wise processing of the domain and combination of base plate of the c08 base-steps.
     */
    lc_c04_combined_SoA,
    /**
     * LCC08Traversal : More compact form of LCC18Traversal. Eight-way domain coloring in minimally small interaction
     * blocks. High degree of parallelism and good load balancing due to fine granularity.
     */
    lc_c08,
    /**
     * LCC18Traversal : More compact form of LCC01Traversal supporting Newton3 by only accessing forward neighbors.
     */
    lc_c18,
    /**
     * LCSlicedTraversal : 1D equidistant slicing of the domain with one slice per thread. One lock per slice interface.
     * Uses c08 base-step per cell. Minimal scheduling overhead at the cost of no load balancing at all.
     */
    lc_sliced,
    /**
     * LCSlicedBalancedTraversal : Same as lc_sliced but tries to balance slice thickness according to a given
     * LoadEstimatorOption.
     */
    lc_sliced_balanced,
    /**
     * LCSlicedC02Traversal : 1D slicing with as many slices of minimal thickness as possible. No locks but two way
     * coloring of slices.
     */
    lc_sliced_c02,

    // Octree Traversals:
    /**
     * OTC01Traversal : Simple DFS traversal without newton 3 optimization
     */
    ot_c01,
    /**
     * OTC18Traversal : DFS traversal with newton 3 optimization that checks whether a neighbor has already been
     * processed via ID comparison
     */
    ot_c18,

    // VerletClusterLists Traversals:
    /**
     * VCLC01BalancedTraversal : Assign a fixed set of towers to each thread balanced by number of contained clusters.
     * Does not support Newton3.
     */
    vcl_c01_balanced,
    /**
     * VCLC06Traversal : Six-way coloring of the 2D ClusterTower grid in the c18 base step style.
     * Rather coarse but dynamically balanced.
     */
    vcl_c06,
    /**
     * VCLClusterIterationTraversal : Dynamically schedule ClusterTower to threads.
     * Does not support Newton3.
     */
    vcl_cluster_iteration,
    /**
     * VCLSlicedTraversal : Equivalent to lc_sliced with slicing applied to the tower grid.
     */
    vcl_sliced,
    /**
     * VCLSlicedBalancedTraversal : Same as vcl_sliced but tries to balance slice thickness according to a given
     * LoadEstimatorOption.
     */
    vcl_sliced_balanced,
    /**
     * VCLSlicedC02Traversal : 1D slicing with as many slices of minimal thickness as possible. No locks but two way
     * coloring of slices.
     */
    vcl_sliced_c02,

    // VerletList Traversals:
    /**
     * VLListIterationTraversal : Distribute processing of neighbor lists dynamically to threads.
     * Does not support Newton3.
     */
    vl_list_iteration,

    // VerletListCells Traversals:
    /**
     * VLCC01Traversal : Equivalent to LCC01Traversal. Schedules all neighbor lists of one cell at once.
     * Does not support Newton3.
     */
    vlc_c01,
    /**
     * VLCC18Traversal : Equivalent to LCC18Traversal. Neighbor lists contain only forward neighbors.
     */
    vlc_c18,
    /**
     * VLCC08Traversal : Equivalent to LCC08Traversal. Base cell contains all neighbor lists of the base step.
     */
    vlc_c08,
    /**
     * VLCSlicedTraversal : Equivalent to LCSlicedTraversal but with a c18 base-step.
     */
    vlc_sliced,
    /**
     * VLCSlicedBalancedTraversal : Equivalent to LCSlicedBalancedTraversal but with a c18 base-step.
     * Tries to balance slice thickness according to a given LoadEstimatorOption.
     */
    vlc_sliced_balanced,
    /**
     * VLCSlicedC02Traversal : Equivalent to LCSlicedC02Traversal.
     * 1D slicing with as many slices of minimal thickness as possible. No locks but two-way coloring of slices.
     */
    vlc_sliced_c02,

    // PairwiseVerletLists Traversals - same traversals as VLC but with a new name for the pairwise container
    /**
     * VLPC01Traversal : Equivalent to LCC01Traversal. Schedules all neighbor lists of one cell at once.
     * Does not support Newton3.
     */
    vlp_c01,
    /**
     * VLPC18Traversal : Equivalent to LCC18Traversal. Neighbor lists contain only forward neighbors.
     */
    vlp_c18,
    /**
     * VLPSlicedTraversal : Equivalent to LCSlicedTraversal.
     */
    vlp_sliced,
    /**
     * VLPSlicedBalancedTraversal : Equivalent to LCSlicedBalancedTraversal.
     * Tries to balance slice thickness according to a given LoadEstimatorOption.
     */
    vlp_sliced_balanced,
    /**
     * VLPSlicedC02Traversal : Equivalent to LCSlicedC02Traversal.
     * 1D slicing with as many slices of minimal thickness as possible. No locks but two-way coloring of slices.
     */
    vlp_sliced_c02,

    /**
     * VLPCellPairC08Traversal : based on LCC08Traversal.
     * The pairwise neighbor list allows access to the relevant pairs of interacting particles for each pair of cells,
     * including the diagonal non-base pair of cells in the standard c08 step.
     */
    vlp_c08,

    // VarVerlet Traversals:
    /**
     * VVLAsBuildTraversal : Track which thread built what neighbor list and schedule them the same way for the pairwise
     * iteration. Provides some kind of load balancing if the force calculation is cheap but is sensitive to hardware
     * fluctuations.
     */
    vvl_as_built,
    // HierarchicalGrid Traversals:
    hgrid_test,
    hgrid_test2,
    hgrid_test3,
    hgrid_test4,
    hgrid_test5,
    hgrid_test6,
    /**
     * For each level, LCC08Traversal is used. For the cross-level interactions, for each level x only smaller levels
     * are iterated (newton3 on only). The cells on level x are iterated with colors (dynamic color count based on ratio
     * of cell lengths between level x and y) so that the cells on the lower level y
     * that are considered for each cell on level x do not intersect.
     * To reduce number of colors and increase memory efficiency, instead of only 1 upper
     * level cell a block of cells is assigned to a thread at a time. The size of block is calculated dynamically
     * by considering upper and lower cell lengths and number of threads. The number of blocks per color is at least
     * num_threads * 4 or 8, depending on the option.
     */
    hgrid_block4,
    hgrid_block8,
    /**
     * Similar to hgrid_block but instead of fully waiting for a color to end to start the next color, openmp task with
     * dependencies is used. The basic idea is that if the cells with the previous color around the cell is computed,
     * the cell with the next color can start computing. The numbers hgrid_taskX denote that the total number of
     * OpenMP tasks should be as close to X * num_threads as possible.
     */
    hgrid_task32,
    hgrid_task64,
    hgrid_task128,
    /**
     * Same as hgrid_block but with SoA cell to cell functor.
     */
    hgrid_block_soa_cell,
    /**
     * Same as hgrid_task but with SoA cell to cell functor.
     */
    hgrid_task_soa_cell,
  };

  /**
   * Constructor.
   */
  TraversalOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TraversalOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<TraversalOption> getDiscouragedOptions() {
    return {Value::ds_sequential, Value::vcl_cluster_iteration};
  }

  /**
   * Set of options that apply for pairwise interactions.
   * @return
   */
  static std::set<TraversalOption> getAllPairwiseOptions() { return getAllOptions(); }

  /**
   * Set of options that apply for triwise interactions.
   * @return
   */
  static std::set<TraversalOption> getAllTriwiseOptions() { return {Value::ds_sequential, Value::lc_c01}; }

  /**
   * Set of all pairwise traversals without discouraged options.
   * @return
   */
  static std::set<TraversalOption> getMostPairwiseOptions() {
    std::set<TraversalOption> mostPairwiseOptions;
    auto allOptions = getAllOptions();
    auto discouragedOptions = getDiscouragedOptions();
    std::set_difference(allOptions.begin(), allOptions.end(), discouragedOptions.begin(), discouragedOptions.end(),
                        std::inserter(mostPairwiseOptions, mostPairwiseOptions.begin()));
    return mostPairwiseOptions;
  }

  /**
   * Set of all triwise traversals without discouraged options.
   * @return
   */
  static std::set<TraversalOption> getMostTriwiseOptions() {
    std::set<TraversalOption> mostTriwiseOptions;
    auto allOptions = getAllTriwiseOptions();
    auto discouragedOptions = getDiscouragedOptions();
    std::set_difference(allOptions.begin(), allOptions.end(), discouragedOptions.begin(), discouragedOptions.end(),
                        std::inserter(mostTriwiseOptions, mostTriwiseOptions.begin()));
    return mostTriwiseOptions;
  }

  /**
   * Set of all options specific to an interaction type.
   * @param interactionType
   * @return
   */
  static std::set<TraversalOption> getAllOptionsOf(const autopas::InteractionTypeOption &interactionType) {
    switch (interactionType) {
      case autopas::InteractionTypeOption::pairwise:
        return getAllPairwiseOptions();
      case autopas::InteractionTypeOption::triwise:
        return getAllTriwiseOptions();
      default:
        return {};
    }
  }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<TraversalOption, std::string> getOptionNames() {
    return {
        // DirectSum Traversals:
        {TraversalOption::ds_sequential, "ds_sequential"},

        // LinkedCell Traversals:
        {TraversalOption::lc_sliced, "lc_sliced"},
        {TraversalOption::lc_sliced_balanced, "lc_sliced_balanced"},
        {TraversalOption::lc_sliced_c02, "lc_sliced_c02"},
        {TraversalOption::lc_c01, "lc_c01"},
        {TraversalOption::lc_c01_combined_SoA, "lc_c01_combined_SoA"},
        {TraversalOption::lc_c04, "lc_c04"},
        {TraversalOption::lc_c04_HCP, "lc_c04_HCP"},
        {TraversalOption::lc_c04_combined_SoA, "lc_c04_combined_SoA"},
        {TraversalOption::lc_c08, "lc_c08"},
        {TraversalOption::lc_c18, "lc_c18"},

        // VerletClusterLists Traversals:
        {TraversalOption::vcl_cluster_iteration, "vcl_cluster_iteration"},
        {TraversalOption::vcl_c06, "vcl_c06"},
        {TraversalOption::vcl_c01_balanced, "vcl_c01_balanced"},
        {TraversalOption::vcl_sliced, "vcl_sliced"},
        {TraversalOption::vcl_sliced_c02, "vcl_sliced_c02"},
        {TraversalOption::vcl_sliced_balanced, "vcl_sliced_balanced"},

        // VerletList Traversals:
        {TraversalOption::vl_list_iteration, "vl_list_iteration"},

        // VerletListCells Traversals:
        {TraversalOption::vlc_sliced, "vlc_sliced"},
        {TraversalOption::vlc_sliced_c02, "vlc_sliced_c02"},
        {TraversalOption::vlc_c18, "vlc_c18"},
        {TraversalOption::vlc_c01, "vlc_c01"},
        {TraversalOption::vlc_c08, "vlc_c08"},
        {TraversalOption::vlc_sliced_balanced, "vlc_sliced_balanced"},

        // VarVerlet Traversals:
        {TraversalOption::vvl_as_built, "vvl_as_built"},

        // PairwiseVerlet Traversals:
        {TraversalOption::vlp_sliced, "vlp_sliced"},
        {TraversalOption::vlp_sliced_c02, "vlp_sliced_c02"},
        {TraversalOption::vlp_c18, "vlp_c18"},
        {TraversalOption::vlp_c01, "vlp_c01"},
        {TraversalOption::vlp_sliced_balanced, "vlp_sliced_balanced"},
        {TraversalOption::vlp_c08, "vlp_c08"},

        // Octree Traversals:
        {TraversalOption::ot_c18, "ot_c18"},
        {TraversalOption::ot_c01, "ot_c01"},

        // HierarchicalGrid Traversals:
        {TraversalOption::hgrid_block_soa_cell, "hgrid_block_soa_cell"},
        {TraversalOption::hgrid_task_soa_cell, "hgrid_task_soa_cell"},
        {TraversalOption::hgrid_task32, "hgrid_task32"},
        {TraversalOption::hgrid_task64, "hgrid_task64"},
        {TraversalOption::hgrid_task128, "hgrid_task128"},
        {TraversalOption::hgrid_block4, "hgrid_block4"},
        {TraversalOption::hgrid_block8, "hgrid_block8"},
        {TraversalOption::hgrid_test, "hgrid_test"},
        {TraversalOption::hgrid_test2, "hgrid_test2"},
        {TraversalOption::hgrid_test3, "hgrid_test3"},
        {TraversalOption::hgrid_test4, "hgrid_test4"},
        {TraversalOption::hgrid_test5, "hgrid_test5"},
        {TraversalOption::hgrid_test6, "hgrid_test6"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
