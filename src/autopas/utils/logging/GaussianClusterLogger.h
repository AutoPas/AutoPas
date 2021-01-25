/**
 * @file GaussianClusterLogger.h
 * @author Jan Nguyen
 * @date 30.06.20
 */

#pragma once

#include <unordered_set>

#include "selectors/tuningStrategy/GaussianModel/GaussianModelTypes.h"
#include "selectors/tuningStrategy/GaussianModel/GaussianProcess.h"

namespace autopas {

/**
 * Used to print out the clusters of GaussianClusters.
 * The resulting graph represents each cluster as a node and the weight between clusters as edges.
 * The graph is printed as two csv-files.
 *
 * By default logging the data is disabled. It can be enabled by setting the cmake variable AUTOPAS_Log_GaussianCluster
 * to ON.
 */
class GaussianClusterLogger {
  const std::string node_start_marker = "GaussianCluster Graph: Nodes";
  const std::string edge_start_marker = "GaussianCluster Graph: Edges";
  const std::string end_marker = "GaussianCluster Graph: End";

  /**
   * Mode of used string streams. Streams are open for input and append these to the end.
   */
  constexpr static auto streamMode = std::ios_base::out | std::ios_base::app;

 public:
  /**
   * Constructor
   * @param vecToStringFun function to convert vectors to readable string
   * @param outputType
   */
  GaussianClusterLogger(GaussianModelTypes::VectorToStringFun vecToStringFun);

  /**
   * Change the used function to convert from vector to string.
   * @param fun new converter
   */
  void setVectorToStringFun(const GaussianModelTypes::VectorToStringFun &fun);

  /**
   * Add nodes and edges for given continuous sample.
   * @param clusters all clusters
   * @param discreteVectorMap map to convert index to vector
   * @param currentContinous continuous sample
   * @param means predicted mean for each cluster
   * @param vars predicted variance for each cluster
   * @param neighbourWeights neighbours for each cluster
   */
  void add(const std::vector<GaussianProcess> &clusters,
           const std::vector<GaussianModelTypes::VectorDiscrete> &discreteVectorMap,
           const GaussianModelTypes::VectorContinuous &currentContinous, const std::vector<double> &means,
           const std::vector<double> &vars, const GaussianModelTypes::NeighboursWeights &neighbourWeights);

  /**
   * Dump all data accumulated by add() to the sink of this logger and clear all buffers.
   */
  void flush();

 private:
  /**
   * Reset the stream for both csv-files.
   */
  void reset();

  std::string _outputFileName;

  /**
   * Stream for the nodes csv-file.
   */
  std::stringstream _nodeStream;
  /**
   * Stream for the edges csv-files
   */
  std::stringstream _edgeStream;

  /**
   * Continuous samples already added to the graph.
   */
  std::vector<GaussianModelTypes::VectorContinuous> _currentContinuous{};

  /**
   * Function to convert vectors to strings.
   */
  GaussianModelTypes::VectorToStringFun _vecToStringFun;
};
}  // namespace autopas
