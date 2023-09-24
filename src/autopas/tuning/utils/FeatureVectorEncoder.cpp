/**
 * @file FeatureVectorEncoder.cpp
 * @author F. Gratl
 * @date 17.11.22
 */

#include "FeatureVectorEncoder.h"

autopas::FeatureVectorEncoder::FeatureVectorEncoder() = default;

autopas::FeatureVectorEncoder::FeatureVectorEncoder(
    const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
    const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
    const autopas::NumberSet<double> &cellSizeFactors) {
  setAllowedOptions(containerTraversalEstimatorOptions, dataLayoutOptions, newton3Options, cellSizeFactors);
}

autopas::FeatureVectorEncoder::~FeatureVectorEncoder() = default;

void autopas::FeatureVectorEncoder::setAllowedOptions(
    const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
    const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
    const autopas::NumberSet<double> &cellSizeFactors) {
  _containerTraversalEstimatorOptions = containerTraversalEstimatorOptions;
  _dataLayoutOptions = dataLayoutOptions;
  _newton3Options = newton3Options;

  _oneHotDims = _containerTraversalEstimatorOptions.size() + _dataLayoutOptions.size() + _newton3Options.size() +
                tunableContinuousDims;

  _discreteRestrictions[static_cast<size_t>(DiscreteIndices::containerTraversalEstimator)] =
      _containerTraversalEstimatorOptions.size();
  _discreteRestrictions[static_cast<size_t>(DiscreteIndices::dataLayout)] = _dataLayoutOptions.size();
  _discreteRestrictions[static_cast<size_t>(DiscreteIndices::newton3)] = _newton3Options.size();

  _continuousRestrictions[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = cellSizeFactors.clone();
}

size_t autopas::FeatureVectorEncoder::getOneHotDims() const { return _oneHotDims; }

const std::array<int, autopas::FeatureVectorEncoder::tunableDiscreteDims>
    &autopas::FeatureVectorEncoder::getDiscreteRestrictions() const {
  return _discreteRestrictions;
}

Eigen::VectorXd autopas::FeatureVectorEncoder::oneHotEncode(const autopas::FeatureVector &vec) const {
  std::vector<double> data;
  data.reserve(_oneHotDims);

  auto [discreteValues, continuousValues] = convertToTunable(vec);

  // discrete values are encoded using one-hot-encoding
  for (size_t i = 0; i < tunableDiscreteDims; ++i) {
    for (int listIndex = 0; listIndex < _discreteRestrictions[i]; ++listIndex) {
      data.push_back(listIndex == discreteValues[i] ? 1. : 0.);
    }
  }

  // continuous values are simply copied
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    data.push_back(continuousValues[i]);
  }

  return Eigen::Map<Eigen::VectorXd>(data.data(), _oneHotDims);
}

autopas::FeatureVector autopas::FeatureVectorEncoder::oneHotDecode(const Eigen::VectorXd &vec) {
  if (static_cast<size_t>(vec.size()) != _oneHotDims) {
    utils::ExceptionHandler::exception("FeatureVectorEncoder.oneHotDecode: Expected size {}, got {}", _oneHotDims,
                                       vec.size());
  }

  size_t pos = 0;

  DiscreteDimensionType discreteValues;

  // extract each one-hot-encoded discrete option
  for (size_t i = 0; i < tunableDiscreteDims; ++i) {
    // for each dimension get index whose value equals 1.
    std::optional<int> value = {};
    for (int listIndex = 0; listIndex < _discreteRestrictions[i]; ++listIndex) {
      if (vec[pos++] == 1.) {
        if (value) {
          utils::ExceptionHandler::exception(
              "FeatureVectorEncoder.oneHotDecode: Vector encodes more than one option for tunable dimension {}. "
              "(More than one value equals 1.)",
              i);
        }
        value = listIndex;
      }
    }
    if (not value) {
      utils::ExceptionHandler::exception(
          "FeatureVectorEncoder.oneHotDecode: Vector encodes no option for tunable dimension {}. (All values equal "
          "0.)",
          i);
    }
    discreteValues[i] = *value;
  }

  ContinuousDimensionType continuousValues;
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    continuousValues[i] = vec[pos++];
  }

  return convertFromTunable(discreteValues, continuousValues);
}

std::pair<Eigen::VectorXi, Eigen::VectorXd> autopas::FeatureVectorEncoder::convertToCluster(
    const autopas::FeatureVector &vec, double iteration) const {
  auto [discreteValues, continuousValues] = convertToTunable(vec);
  Eigen::Map<Eigen::VectorXi> vecDiscrete(discreteValues.data(), tunableDiscreteDims);

  std::vector<double> continuousData;
  continuousData.reserve(tunableContinuousDims + 1);
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    continuousData.push_back(continuousValues[i]);
  }
  continuousData.push_back(iteration);
  Eigen::Map<Eigen::VectorXd> vecContinuous(continuousData.data(), tunableContinuousDims + 1);
  return std::make_pair(vecDiscrete, vecContinuous);
}

autopas::FeatureVector autopas::FeatureVectorEncoder::convertFromCluster(
    const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
  const auto &[vecDiscrete, vecContinuous] = vec;

  DiscreteDimensionType discreteValues;
  for (size_t i = 0; i < tunableDiscreteDims; ++i) {
    discreteValues[i] = vecDiscrete[i];
  }

  ContinuousDimensionType continuousValues;
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    continuousValues[i] = vecContinuous[i];
  }

  return convertFromTunable(discreteValues, continuousValues);
}

std::vector<std::pair<Eigen::VectorXi, double>> autopas::FeatureVectorEncoder::clusterNeighboursManhattan1(
    const Eigen::VectorXi &target) {
  std::vector<std::pair<Eigen::VectorXi, double>> result;
  // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
  // dimension)
  result.reserve(
      std::accumulate(_discreteRestrictions.begin(), _discreteRestrictions.end(), -_discreteRestrictions.size()));

  // for each dimension
  for (int i = 0; i < target.size(); ++i) {
    // initial value
    auto init = target[i];

    // for each possible value of that dimension
    for (int x = 0; x < _discreteRestrictions[i]; ++x) {
      // skip initial value
      if (x != init) {
        auto neighbour = target;
        neighbour[i] = x;
        result.emplace_back(std::move(neighbour), 1.);
      }
    }
  }
  return result;
}

std::vector<std::pair<Eigen::VectorXi, double>> autopas::FeatureVectorEncoder::clusterNeighboursManhattan1Container(
    const Eigen::VectorXi &target) {
  std::vector<std::pair<Eigen::VectorXi, double>> result;
  // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
  // dimension)
  result.reserve(
      std::accumulate(_discreteRestrictions.begin(), _discreteRestrictions.end(), -_discreteRestrictions.size()));

  auto targetContainer = std::get<0>(
      _containerTraversalEstimatorOptions[target[static_cast<int>(DiscreteIndices::containerTraversalEstimator)]]);

  // for each dimension
  for (int i = 0; i < target.size(); ++i) {
    // initial value
    auto init = target[i];

    // for each possible value of that dimension
    for (int x = 0; x < _discreteRestrictions[i]; ++x) {
      // skip initial value
      if (x != init) {
        auto neighbour = target;
        neighbour[i] = x;
        double weight = 1.;
        // check if container changed
        if (i == static_cast<int>(DiscreteIndices::containerTraversalEstimator)) {
          auto xContainer = std::get<0>(_containerTraversalEstimatorOptions[x]);
          if (targetContainer != xContainer) {
            // assign lower weight
            weight = 0.5;
          }
        }

        result.emplace_back(std::move(neighbour), weight);
      }
    }
  }
  return result;
}

std::vector<autopas::FeatureVector> autopas::FeatureVectorEncoder::lhsSampleFeatures(size_t n,
                                                                                     autopas::Random &rng) const {
  // create n samples for each discrete dimension.
  std::array<std::vector<size_t>, tunableDiscreteDims> lhsDiscreteSamples;
  for (size_t i = 0; i < tunableDiscreteDims; ++i) {
    lhsDiscreteSamples[i] = rng.uniformSample(0, _discreteRestrictions[i] - 1, n);
  }

  // create n samples for each continuous dimension.
  std::array<std::vector<double>, tunableContinuousDims> lhsContinuousSamples;
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
  }

  // create FeatureVectors from raw samples
  std::vector<FeatureVector> result;
  for (size_t i = 0; i < n; ++i) {
    DiscreteDimensionType discreteValues;
    for (size_t d = 0; d < tunableDiscreteDims; ++d) {
      discreteValues[d] = lhsDiscreteSamples[d][i];
    }

    ContinuousDimensionType continuousValues;
    for (size_t c = 0; c < tunableContinuousDims; ++c) {
      continuousValues[c] = lhsContinuousSamples[c][i];
    }

    result.emplace_back(convertFromTunable(discreteValues, continuousValues));
  }

  return result;
}

std::vector<Eigen::VectorXd> autopas::FeatureVectorEncoder::lhsSampleFeatureCluster(size_t n, autopas::Random &rng,
                                                                                    double iteration) const {
  // create n samples for each continuous dimension.
  std::array<std::vector<double>, tunableContinuousDims> lhsContinuousSamples;
  for (size_t i = 0; i < tunableContinuousDims; ++i) {
    lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
  }

  // create FeatureVectors from raw samples
  std::vector<Eigen::VectorXd> result;
  for (size_t i = 0; i < n; ++i) {
    std::array<double, tunableContinuousDims + 1> data;
    for (size_t c = 0; c < tunableContinuousDims; ++c) {
      data[c] = lhsContinuousSamples[c][i];
    }
    data[tunableContinuousDims] = iteration;
    result.emplace_back(Eigen::Map<Eigen::VectorXd>(data.data(), data.size()));
  }

  return result;
}

std::pair<autopas::FeatureVectorEncoder::DiscreteDimensionType, autopas::FeatureVectorEncoder::ContinuousDimensionType>
autopas::FeatureVectorEncoder::convertToTunable(const autopas::FeatureVector &vec) const {
  DiscreteDimensionType discreteValues;
  discreteValues[static_cast<size_t>(DiscreteIndices::containerTraversalEstimator)] =
      getIndex(_containerTraversalEstimatorOptions, std::make_tuple(vec.container, vec.traversal, vec.loadEstimator));
  discreteValues[static_cast<size_t>(DiscreteIndices::dataLayout)] = getIndex(_dataLayoutOptions, vec.dataLayout);
  discreteValues[static_cast<size_t>(DiscreteIndices::newton3)] = getIndex(_newton3Options, vec.newton3);

  ContinuousDimensionType continuousValues;
  continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = vec.cellSizeFactor;

  return std::make_pair(discreteValues, continuousValues);
}

autopas::FeatureVector autopas::FeatureVectorEncoder::convertFromTunable(
    const autopas::FeatureVectorEncoder::DiscreteDimensionType &discreteValues,
    const autopas::FeatureVectorEncoder::ContinuousDimensionType &continuousValues) const {
  const auto &[container, traversal, estimator] =
      _containerTraversalEstimatorOptions[discreteValues[static_cast<size_t>(
          DiscreteIndices::containerTraversalEstimator)]];
  auto dataLayout = _dataLayoutOptions[discreteValues[static_cast<size_t>(DiscreteIndices::dataLayout)]];
  auto newton3 = _newton3Options[discreteValues[static_cast<size_t>(DiscreteIndices::newton3)]];

  auto cellSizeFactor = continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)];

  return FeatureVector(container, cellSizeFactor, traversal, estimator, dataLayout, newton3);
}
