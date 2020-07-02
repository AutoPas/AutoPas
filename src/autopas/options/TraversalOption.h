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
    ds_sequential,

    // LinkedCell Traversals:
    lc_sliced,
    lc_sliced_balanced,
    lc_c01,
    lc_c01_combined_SoA,
    lc_c01_cuda,
    lc_c04,
    lc_c04_HCP,
    lc_c04_combined_SoA,
    lc_c08,
    lc_c18,

    // VerletClusterCells Traversals:
    vcc_cluster_iteration,

    // VerletClusterLists Traversals:
    vcl_cluster_iteration,
    vcl_c06,
    vcl_c01_balanced,

    // VerletList Traversals:
    vl_list_iteration,

    // VerletListCells Traversals:
    vlc_c01,
    vlc_c18,
    vlc_sliced,
    vlc_sliced_balanced,

    // VarVerlet Traversals:
    vvl_as_built,
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
        {TraversalOption::lc_c01, "lc_c01"},
        {TraversalOption::lc_c01_cuda, "lc_c01_cuda"},
        {TraversalOption::lc_c01_combined_SoA, "lc_c01_combined_SoA"},
        {TraversalOption::lc_c04, "lc_c04"},
        {TraversalOption::lc_c04_HCP, "lc_c04_HCP"},
        {TraversalOption::lc_c04_combined_SoA, "lc_c04_combined_SoA"},
        {TraversalOption::lc_c08, "lc_c08"},
        {TraversalOption::lc_c18, "lc_c18"},

        // VerletClusterCells Traversals:
        {TraversalOption::vcc_cluster_iteration, "vcc_cluster_iteration"},

        // VerletClusterLists Traversals:
        {TraversalOption::vcl_cluster_iteration, "vcl_cluster_iteration"},
        {TraversalOption::vcl_c06, "vcl_c06"},
        {TraversalOption::vcl_c01_balanced, "vcl_c01_balanced"},

        // VerletList Traversals:
        {TraversalOption::vl_list_iteration, "vl_list_iteration"},

        // VerletListCells Traversals:
        {TraversalOption::vlc_sliced, "vlc_sliced"},
        {TraversalOption::vlc_c18, "vlc_c18"},
        {TraversalOption::vlc_c01, "vlc_c01"},
        {TraversalOption::vlc_sliced_balanced, "vlc_sliced_balanced"},

        // VarVerlet Traversals:
        {TraversalOption::vvl_as_built, "vvl_as_built"},

    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
