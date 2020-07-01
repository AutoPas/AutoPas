/**
 * @file GaussianModelTypes.h
 * @author Jan Nguyen
 * @date 01.07.20
 */

#pragma once

#include <Eigen/Core>
#include <functional>
#include <utility>

/**
 * Aliases shared between GaussianModel based files.
 */
namespace autopas::GaussianModelTypes {

using VectorDiscrete = Eigen::VectorXi;
using VectorContinuous = Eigen::VectorXd;

/**
 * store pairs of vectors and corresponding acquisition
 */
using VectorAcquisition = std::pair<std::pair<VectorDiscrete, VectorContinuous>, double>;

/**
 * function that generate all neighbouring vectors of given vector
 */
using NeighbourFunction = std::function<std::vector<VectorDiscrete>(VectorDiscrete)>;

/**
 * for each vector store a vector of all neighbours and their corresponding weight
 */
using NeighboursWeights = std::vector<std::vector<std::pair<size_t, double>>>;

/**
 * Vector described by a discrete and a continuous part
 */
using VectorPairDiscreteContinuous = std::pair<VectorDiscrete, VectorContinuous>;

/**
 * function to convert a vector to a string.
 */
using VectorToStringFun = std::function<std::string(const VectorPairDiscreteContinuous &)>;

}  // namespace autopas::GaussianModelTypes
