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
    const NumberSet<double> &cellSizeFactors, const std::vector<OpenMPKindOption> &ompKindOptions,
    const NumberSet<size_t> &ompChunkSizes, const InteractionTypeOption &interactionType)
    : _interactionType(interactionType) {
  setAllowedOptions(containerTraversalEstimatorOptions, dataLayoutOptions, newton3Options, cellSizeFactors,
    ompKindOptions, ompChunkSizes);
}

autopas::FeatureVectorEncoder::~FeatureVectorEncoder() = default;

void autopas::FeatureVectorEncoder::setAllowedOptions(
    const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
    const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
    const NumberSet<double> &cellSizeFactors, const std::vector<OpenMPKindOption> &ompKindOptions,
    const NumberSet<size_t> &ompChunkSizes) {
  _containerTraversalEstimatorOptions = containerTraversalEstimatorOptions;
  _dataLayoutOptions = dataLayoutOptions;
  _newton3Options = newton3Options;
  _openMPKindOptions = ompKindOptions;

  _oneHotDims = _containerTraversalEstimatorOptions.size() + _dataLayoutOptions.size() + _newton3Options.size() +
                _openMPKindOptions.size() + numTunableContinuousDims;

  _discreteRestrictions[static_cast<size_t>(OneHotEncodedIndices::containerTraversalEstimator)] =
      _containerTraversalEstimatorOptions.size();
  _discreteRestrictions[static_cast<size_t>(OneHotEncodedIndices::dataLayout)] = _dataLayoutOptions.size();
  _discreteRestrictions[static_cast<size_t>(OneHotEncodedIndices::newton3)] = _newton3Options.size();
  _discreteRestrictions[static_cast<size_t>(OneHotEncodedIndices::ompKind)] = _openMPKindOptions.size();

  _continuousRestrictions[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = cellSizeFactors.clone();

  // Convert ompChunkSizes from NumberSet<size_t> to NumberSet<double>
  const auto ompChunkSizesSize_t = ompChunkSizes.getAll();
  std::set<double> ompChunkSizesDouble;
  for (const auto &size : ompChunkSizesSize_t) {
    ompChunkSizesDouble.insert(static_cast<double>(size));
  }
  _continuousRestrictions[static_cast<size_t>(ContinuousIndices::ompChunkSize)] =
      std::make_unique<NumberSetFinite<double>>(ompChunkSizesDouble);
}

size_t autopas::FeatureVectorEncoder::getOneHotDims() const { return _oneHotDims; }

const std::array<int, autopas::FeatureVectorEncoder::numTunableOneHotIndices>
    &autopas::FeatureVectorEncoder::getDiscreteRestrictions() const {
  return _discreteRestrictions;
}

Eigen::VectorXd autopas::FeatureVectorEncoder::oneHotEncode(const autopas::FeatureVector &vec) const {
  std::vector<double> data;
  data.reserve(_oneHotDims);

  auto [discreteValues, continuousValues] = convertToTunable(vec);

  // discrete values are encoded using one-hot-encoding
  for (size_t i = 0; i < numTunableOneHotIndices; ++i) {
    for (int listIndex = 0; listIndex < _discreteRestrictions[i]; ++listIndex) {
      data.push_back(listIndex == discreteValues[i] ? 1. : 0.);
    }
  }

  // continuous values are simply copied
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
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

  OneHotIndicesType discreteValues;

  // extract each one-hot-encoded discrete option
  for (size_t i = 0; i < numTunableOneHotIndices; ++i) {
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

  ContinuousIndicesType continuousValues;
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
    continuousValues[i] = vec[pos++];
  }

  return convertFromTunable(discreteValues, continuousValues);
}

std::pair<Eigen::VectorXi, Eigen::VectorXd> autopas::FeatureVectorEncoder::convertToCluster(
    const autopas::FeatureVector &vec, double iteration) const {
  auto [discreteValues, continuousValues] = convertToTunable(vec);
  Eigen::Map<Eigen::VectorXi> vecDiscrete(discreteValues.data(), numTunableOneHotIndices);

  std::vector<double> continuousData;
  continuousData.reserve(numTunableContinuousDims + 1);
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
    continuousData.push_back(continuousValues[i]);
  }
  continuousData.push_back(iteration);
  Eigen::Map<Eigen::VectorXd> vecContinuous(continuousData.data(), numTunableContinuousDims + 1);
  return std::make_pair(vecDiscrete, vecContinuous);
}

autopas::FeatureVector autopas::FeatureVectorEncoder::convertFromCluster(
    const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
  const auto &[vecDiscrete, vecContinuous] = vec;

  OneHotIndicesType discreteValues;
  for (size_t i = 0; i < numTunableOneHotIndices; ++i) {
    discreteValues[i] = vecDiscrete[i];
  }

  ContinuousIndicesType continuousValues;
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
    continuousValues[i] = vecContinuous[i];
  }

  return convertFromTunable(discreteValues, continuousValues);
}

std::vector<std::pair<Eigen::VectorXi, double>> autopas::FeatureVectorEncoder::clusterNeighborsManhattan1(
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

std::vector<std::pair<Eigen::VectorXi, double>> autopas::FeatureVectorEncoder::clusterNeighborsManhattan1Container(
    const Eigen::VectorXi &target) {
  std::vector<std::pair<Eigen::VectorXi, double>> result;
  // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
  // dimension)
  result.reserve(
      std::accumulate(_discreteRestrictions.begin(), _discreteRestrictions.end(), -_discreteRestrictions.size()));

  auto targetContainer = std::get<0>(
      _containerTraversalEstimatorOptions[target[static_cast<int>(OneHotEncodedIndices::containerTraversalEstimator)]]);

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
        if (i == static_cast<int>(OneHotEncodedIndices::containerTraversalEstimator)) {
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
  std::array<std::vector<size_t>, numTunableOneHotIndices> lhsDiscreteSamples;
  for (size_t i = 0; i < numTunableOneHotIndices; ++i) {
    lhsDiscreteSamples[i] = rng.uniformSample(0, _discreteRestrictions[i] - 1, n);
  }

  // create n samples for each continuous dimension.
  std::array<std::vector<double>, numTunableContinuousDims> lhsContinuousSamples;
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
    lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
  }

  // create FeatureVectors from raw samples
  std::vector<FeatureVector> result;
  for (size_t i = 0; i < n; ++i) {
    OneHotIndicesType discreteValues;
    for (size_t d = 0; d < numTunableOneHotIndices; ++d) {
      discreteValues[d] = lhsDiscreteSamples[d][i];
    }

    ContinuousIndicesType continuousValues;
    for (size_t c = 0; c < numTunableContinuousDims; ++c) {
      continuousValues[c] = lhsContinuousSamples[c][i];
    }

    result.emplace_back(convertFromTunable(discreteValues, continuousValues));
  }

  return result;
}

std::vector<Eigen::VectorXd> autopas::FeatureVectorEncoder::lhsSampleFeatureCluster(size_t n, autopas::Random &rng,
                                                                                    double iteration) const {
  // create n samples for each continuous dimension.
  std::array<std::vector<double>, numTunableContinuousDims> lhsContinuousSamples;
  for (size_t i = 0; i < numTunableContinuousDims; ++i) {
    lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
  }

  // create FeatureVectors from raw samples
  std::vector<Eigen::VectorXd> result;
  for (size_t i = 0; i < n; ++i) {
    std::array<double, numTunableContinuousDims + 1> data;
    for (size_t c = 0; c < numTunableContinuousDims; ++c) {
      data[c] = lhsContinuousSamples[c][i];
    }
    data[numTunableContinuousDims] = iteration;
    result.emplace_back(Eigen::Map<Eigen::VectorXd>(data.data(), data.size()));
  }

  return result;
}

std::pair<autopas::FeatureVectorEncoder::OneHotIndicesType, autopas::FeatureVectorEncoder::ContinuousIndicesType>
autopas::FeatureVectorEncoder::convertToTunable(const autopas::FeatureVector &vec) const {
  OneHotIndicesType discreteValues;
  discreteValues[static_cast<size_t>(OneHotEncodedIndices::containerTraversalEstimator)] =
      getIndex(_containerTraversalEstimatorOptions, std::make_tuple(vec.container, vec.traversal, vec.loadEstimator));
  discreteValues[static_cast<size_t>(OneHotEncodedIndices::dataLayout)] = getIndex(_dataLayoutOptions, vec.dataLayout);
  discreteValues[static_cast<size_t>(OneHotEncodedIndices::newton3)] = getIndex(_newton3Options, vec.newton3);

  ContinuousIndicesType continuousValues;
  continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = vec.cellSizeFactor;

  return std::make_pair(discreteValues, continuousValues);
}

autopas::FeatureVector autopas::FeatureVectorEncoder::convertFromTunable(
    const autopas::FeatureVectorEncoder::OneHotIndicesType &discreteValues,
    const autopas::FeatureVectorEncoder::ContinuousIndicesType &continuousValues) const {
  const auto &[container, traversal, estimator] =
      _containerTraversalEstimatorOptions[discreteValues[static_cast<size_t>(
          OneHotEncodedIndices::containerTraversalEstimator)]];
  auto dataLayout = _dataLayoutOptions[discreteValues[static_cast<size_t>(OneHotEncodedIndices::dataLayout)]];
  auto newton3 = _newton3Options[discreteValues[static_cast<size_t>(OneHotEncodedIndices::newton3)]];
  auto ompKind = _openMPKindOptions[discreteValues[static_cast<size_t>(OneHotEncodedIndices::ompKind)]];

  auto cellSizeFactor = continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)];
  auto ompChunkSizeSigned = std::llround(continuousValues[static_cast<size_t>(ContinuousIndices::ompChunkSize)]);
  if (ompChunkSizeSigned <= 0) {
    AutoPasLog(WARN, "FeatureVectorEncoder is trying to convert a non-positive OMP chunk size to a size_t. The chunk size will be changed to 1.");
    ompChunkSizeSigned = 1;
  }


  return FeatureVector(container, cellSizeFactor, traversal, estimator, dataLayout, newton3, ompKind, ompChunkSizeSigned, _interactionType);
}
