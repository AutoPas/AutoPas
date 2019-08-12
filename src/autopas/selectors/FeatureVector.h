/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * FeatureVector representation of a Configuration
 */
class FeatureVector : public Configuration {
 public:
  /**
   * Number of tuneable dimensions
   */
  static constexpr size_t featureSpaceDims = 4;

  /**
   * Dimensions of a one-hot-encoded vector
   */
  inline static size_t oneHotDims =
      1 + allTraversalOptions.size() + allDataLayoutOptions.size() + allNewton3Options.size();

  /**
   * Default constructor. Results in invalid vector.
   */
  FeatureVector() : Configuration() {}

  /**
   * Constructor
   * @param _container
   * @param _traversal
   * @param _dataLayout
   * @param _newton3
   * @param _cellSizeFactor
   */
  FeatureVector(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                DataLayoutOption _dataLayout, Newton3Option _newton3)
      : Configuration(_container, _cellSizeFactor, _traversal, _dataLayout, _newton3) {}

  /**
   * Construct from Configuration
   * @param conf
   */
  FeatureVector(Configuration conf) : Configuration(conf) {}

  /**
   * Distance between two FeatureVectors
   * @param other
   * @return
   */
  FeatureVector operator-(const FeatureVector &other) const {
    ContainerOption co = ContainerOption((traversal == other.traversal) ? 0 : 1);
    double cfs = (cellSizeFactor - other.cellSizeFactor);
    TraversalOption to = TraversalOption((traversal == other.traversal) ? 0 : 1);
    DataLayoutOption dlo = DataLayoutOption((dataLayout == other.dataLayout) ? 0 : 1);
    Newton3Option n3o = Newton3Option((newton3 == other.newton3) ? 0 : 1);
    return FeatureVector(co, cfs, to, dlo, n3o);
  }

  /**
   * Cast to Eigen::VectorXd ignoring ContainerOption
   * @return
   */
  operator Eigen::VectorXd() const {
    Eigen::VectorXd result(featureSpaceDims);
    result << cellSizeFactor, static_cast<double>(traversal), static_cast<double>(dataLayout),
        static_cast<double>(newton3);

    return result;
  }

  /**
   * Encode to Eigen::VectorXd ignoring ContainerOption using one-hot-encoding
   * @return one-hot-encoded vector
   */
  Eigen::VectorXd oneHotEncode() const {
    std::vector<double> data;
    data.reserve(oneHotDims);

    data.push_back(cellSizeFactor);
    for (auto to : allTraversalOptions) {
      data.push_back((to == traversal) ? 1. : 0.);
    }
    for (auto dlo : allDataLayoutOptions) {
      data.push_back((dlo == dataLayout) ? 1. : 0.);
    }
    for (auto n3o : allNewton3Options) {
      data.push_back((n3o == newton3) ? 1. : 0.);
    }

    return Eigen::Map<Eigen::VectorXd>(data.data(), oneHotDims);
  }

  /**
   * Decode one-hot-encoded VectorXd to FeatureVector
   * @param vec one-hot-encoded vector
   * @return decoded FeatureVector
   */
  static FeatureVector oneHotDecode(Eigen::VectorXd vec) {
    if (static_cast<size_t>(vec.size()) != oneHotDims) {
      utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: Expected size {}, got {}", oneHotDims,
                                         vec.size());
    }

    size_t pos = 0;
    double cellSizeFactor = vec[pos++];

    // get traversal
    std::optional<TraversalOption> traversal = {};
    for (auto to : allTraversalOptions) {
      if (vec[pos++] == 1.) {
        if (traversal) {
          utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: More than one value for traversal equals 1.");
        }
        traversal = to;
      }
    }
    if (not traversal) {
      utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: All values for traversal equal 0.");
    }

    // get data layout
    std::optional<DataLayoutOption> dataLayout = {};
    for (auto dlo : allDataLayoutOptions) {
      if (vec[pos++] == 1.) {
        if (dataLayout) {
          utils::ExceptionHandler::exception(
              "FeatureVector.oneHotDecode: More than one value for dataLayout equals 1.");
        }
        dataLayout = dlo;
      }
    }
    if (not dataLayout) {
      utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: All values for dataLayout equal 0.");
    }

    // get newton3
    std::optional<Newton3Option> newton3 = {};
    for (auto n3o : allNewton3Options) {
      if (vec[pos++] == 1.) {
        if (newton3) {
          utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: More than one value for newton3 equals 1.");
        }
        newton3 = n3o;
      }
    }
    if (not newton3) {
      utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: All values for newton3 equal 0.");
    }

    return FeatureVector(ContainerOption(-1), cellSizeFactor, *traversal, *dataLayout, *newton3);
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace.
   * Container Option of samples are set to -1, because tuning currently
   * ignores this option.
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param traversals
   * @param dataLayouts
   * @param newton3
   * @return vector of sample featureVectors
   */
  static std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                                                      const std::set<TraversalOption> &traversals,
                                                      const std::set<DataLayoutOption> &dataLayouts,
                                                      const std::set<Newton3Option> &newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto tr = rng.uniformSample(traversals, n);
    auto dl = rng.uniformSample(dataLayouts, n);
    auto n3 = rng.uniformSample(newton3, n);

    std::vector<FeatureVector> result;
    for (unsigned i = 0; i < n; ++i) {
      result.emplace_back(ContainerOption(-1), csf[i], tr[i], dl[i], n3[i]);
    }

    return result;
  }
};
}  // namespace autopas
