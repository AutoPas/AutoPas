/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the traversal choices.
 */
class TraversalOption : public Option<TraversalOption> {
 public:
  /**
   * Possible choices for the cell pair traversal.
   */
  enum Value {
    // DirectSum Traversals:
    /**
     * DSSequentialTraversal : Sequential double loop over all particles.
     */
    ds_sequential,

    // LinkedCell Traversals:
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
    /**
     * LCC01Traversal : Every cell interacts with all neighbors. Is not compatible with Newton3 thus embarrassingly
     * parallel. Good load balancing and no overhead.
     */
    lc_c01,
    /**
     * LCC01Traversal : Same as LCC01Traversal but SoAs are combined into a circular buffer and the domain is traversed
     * line-wise.
     */
    lc_c01_combined_SoA,
    /**
     * LCC01CudaTraversal : CUDA version of LCC01Traversal.
     */
    lc_c01_cuda,
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

    // VerletClusterCells Traversals:
    /**
     * VCCClusterIterationCUDATraversal : CUDA. Concurrent processing of all clusters avoiding data races through
     * atomics.
     */
    vcc_cluster_iteration_cuda,

    // VerletClusterLists Traversals:
    /**
     * VCLClusterIterationTraversal : Schedule ClusterTower to threads.
     */
    vcl_cluster_iteration,
    /**
     * VCLC06Traversal : Six-way coloring of the 2D ClusterTower grid in the c18 base step style.
     * Rather coarse but dynamically balanced.
     */
    vcl_c06,
    /**
     * VCLC01BalancedTraversal : Assign a fixed set of towers to each thread balanced by number of contained clusters.
     * Does not support Newton3.
     */
    vcl_c01_balanced,
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
     * VCCSlicedC02Traversal : 1D slicing with as many slices of minimal thickness as possible. No locks but two way
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

    // VarVerlet Traversals:
    /**
     * VVLAsBuildTraversal : Track which thread built what neighbor list and schedule them the same way for the pairwise
     * iteration. Provides some kind of load balancing if the force calculation is cheap but is sensitive to hardware
     * fluctuations.
     */
    vvl_as_built,

    // PairwiseVerletLists Traversals - same traversals as VLC but with a new name for the pairwise container
    /**
     * VLCC01Traversal : Equivalent to LCC01Traversal. Schedules all neighbor lists of one cell at once.
     * Does not support Newton3.
     */
    vlp_c01,
    /**
     * VLCC18Traversal : Equivalent to LCC18Traversal. Neighbor lists contain only forward neighbors.
     */
    vlp_c18,
    /**
     * VLCSlicedTraversal : Equivalent to LCSlicedTraversal.
     */
    vlp_sliced,
    /**
     * VLCSlicedBalancedTraversal : Equivalent to LCSlicedBalancedTraversal.
     * Tries to balance slice thickness according to a given LoadEstimatorOption.
     */
    vlp_sliced_balanced,
    /**
     * VLCSlicedC02Traversal : Equivalent to LCSlicedC02Traversal.
     * 1D slicing with as many slices of minimal thickness as possible. No locks but two-way coloring of slices.
     */
    vlp_sliced_c02,

    // Octree Traversals
    // TODO(johannes): Documentation
    ot_naive,
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
        {TraversalOption::lc_c01_cuda, "lc_c01_cuda"},
        {TraversalOption::lc_c01_combined_SoA, "lc_c01_combined_SoA"},
        {TraversalOption::lc_c04, "lc_c04"},
        {TraversalOption::lc_c04_HCP, "lc_c04_HCP"},
        {TraversalOption::lc_c04_combined_SoA, "lc_c04_combined_SoA"},
        {TraversalOption::lc_c08, "lc_c08"},
        {TraversalOption::lc_c18, "lc_c18"},

        // VerletClusterCells Traversals:
        {TraversalOption::vcc_cluster_iteration_cuda, "vcc_cluster_iteration_cuda"},

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
        {TraversalOption::vlc_sliced_balanced, "vlc_sliced_balanced"},

        // VarVerlet Traversals:
        {TraversalOption::vvl_as_built, "vvl_as_built"},

        // PairwiseVerlet Traversals:
        {TraversalOption::vlp_sliced, "vlp_sliced"},
        {TraversalOption::vlp_sliced_c02, "vlp_sliced_c02"},
        {TraversalOption::vlp_c18, "vlp_c18"},
        {TraversalOption::vlp_c01, "vlp_c01"},
        {TraversalOption::vlp_sliced_balanced, "vlp_sliced_balanced"},

        // Octree Traversals:
        {TraversalOption::ot_naive, "ot_naive"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
