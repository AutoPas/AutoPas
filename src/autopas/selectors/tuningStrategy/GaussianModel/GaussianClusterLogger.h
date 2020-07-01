/**
 * @file GaussianClusterLogger.h
 * @author Jan Nguyen
 * @date 30.06.20
 */

#pragma once

#include <functional>
#include <iostream>
#include <sstream>

#include "GaussianModelTypes.h"
#include "GaussianProcess.h"
#include "autopas/utils/Logger.h"

namespace autopas {

/**
 * Used to print out the clusters of GaussianClusters.
 * The resulting graph represents each cluster as a node and the weight between clusters as edges.
 * The graph is printed as two csv-files.
 */
class GaussianClusterLogger {
  const std::string node_start_marker = "GaussianCluster Graph: Nodes";
  const std::string edge_start_marker = "GaussianCluster Graph: Edges";
  const std::string end_marker = "GaussianCluster Graph: End";

 public:
  /**
   * Constructor
   * @param vecToStringFun function to convert vectors to readable string
   */
  GaussianClusterLogger(GaussianModelTypes::VectorToStringFun vecToStringFun) : _vecToStringFun(vecToStringFun) {
    if (autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace) {
      return;
    }

    _nodeStream << node_start_marker << std::endl;
    _nodeStream << "Label,mean,var,numEvidence" << std::endl;

    _edgeStream << edge_start_marker << std::endl;
    _edgeStream << "Source,Target,Weight" << std::endl;
  }

  /**
   * Add nodes and edges for given continous sample.
   * @param clusters all clusters
   * @param discreteVectorMap map to convert index to vector
   * @param currentContinous continuous sample
   * @param means predicted mean for each cluster
   * @param vars predicted variance for each cluster
   * @param neighbourWeights neighbours for each cluster
   */
  void add(const std::vector<GaussianProcess> &clusters,
           const std::vector<GaussianModelTypes::VectorDiscrete> &discreteVectorMap,
           const GaussianModelTypes::VectorContinuous &currentContinous, std::vector<double> means,
           std::vector<double> vars, const GaussianModelTypes::NeighboursWeights &neighbourWeights) {
    if (autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace) {
      return;
    }

    // log nodes
    for (size_t i = 0; i < clusters.size(); ++i) {
      std::string label = _vecToStringFun(std::make_pair(discreteVectorMap[i], currentContinous));
      size_t numEvidence = clusters[i].numEvidence();
      _nodeStream << "\"" << label << "\"," << means[i] << "," << vars[i] << "," << numEvidence << std::endl;
    }

    // log edges
    for (size_t i = 0; i < clusters.size(); ++i) {
      for (const auto &[n, weight] : neighbourWeights[i]) {
        if (weight > 0) {
          std::string source = _vecToStringFun(std::make_pair(discreteVectorMap[i], currentContinous));
          std::string target = _vecToStringFun(std::make_pair(discreteVectorMap[n], currentContinous));
          _edgeStream << "\"" << source << "\",\"" << target << "\"," << weight << std::endl;
        }
      }
    }
  }

  /**
   * Output graph stream in AutoPasLog trace.
   */
  void end() const {
    if (autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace) {
      return;
    }

    AutoPasLog(trace, _nodeStream.str());
    AutoPasLog(trace, _edgeStream.str());
    AutoPasLog(trace, end_marker);
  }

  /**
   * Stream for the nodes csv-file.
   */
  std::stringstream _nodeStream;
  /**
   * Stream for the edges csv-files
   */
  std::stringstream _edgeStream;

  /**
   * Function to convert vectors to strings.
   */
  GaussianModelTypes::VectorToStringFun _vecToStringFun;
};
}  // namespace autopas
