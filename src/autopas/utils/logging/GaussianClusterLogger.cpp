/**
 * @file GaussianClusterLogger.cpp
 * @author Jan Nguyen
 * @date 11.08.20
 */

#include "GaussianClusterLogger.h"

#include <fstream>

#include "utils/Timer.h"
#include "utils/logging/Logger.h"

namespace autopas {

GaussianClusterLogger::GaussianClusterLogger(GaussianModelTypes::VectorToStringFun vecToStringFun)
    : _nodeStream(streamMode), _edgeStream(streamMode), _vecToStringFun(std::move(vecToStringFun)) {
#ifdef AUTOPAS_Log_GaussianCluster
  _outputFileName = "gaussianCluster_graph_" + utils::Timer::getDateStamp() + ".out";

  reset();
#endif
}

void GaussianClusterLogger::setVectorToStringFun(const GaussianModelTypes::VectorToStringFun &fun) {
  _vecToStringFun = fun;
}

void GaussianClusterLogger::reset() {
#ifdef AUTOPAS_Log_GaussianCluster
  _nodeStream.str(node_start_marker);
  _nodeStream << std::endl << "Label,mean,var,numEvidence" << std::endl;
  _edgeStream.str(edge_start_marker);
  _edgeStream << std::endl << "Source,Target,Weight" << std::endl;

  _currentContinuous.clear();
#endif
}

void GaussianClusterLogger::add(const std::vector<GaussianProcess> &clusters,
                                const std::vector<GaussianModelTypes::VectorDiscrete> &discreteVectorMap,
                                const GaussianModelTypes::VectorContinuous &currentContinuous,
                                const std::vector<double> &means, const std::vector<double> &vars,
                                const GaussianModelTypes::NeighboursWeights &neighbourWeights) {
#ifdef AUTOPAS_Log_GaussianCluster
  // skip continuous values already added
  if (std::find(_currentContinuous.begin(), _currentContinuous.end(), currentContinuous) != _currentContinuous.end()) {
    return;
  }

  _currentContinuous.push_back(currentContinuous);

  // log nodes
  for (size_t i = 0; i < clusters.size(); ++i) {
    std::string label = _vecToStringFun(std::make_pair(discreteVectorMap[i], currentContinuous));
    size_t numEvidence = clusters[i].numEvidence();
    _nodeStream << "\"" << label << "\"," << means[i] << "," << vars[i] << "," << numEvidence << std::endl;
  }

  // log edges
  for (size_t i = 0; i < clusters.size(); ++i) {
    for (const auto &[n, _, weight] : neighbourWeights[i]) {
      if (weight > 0) {
        std::string source = _vecToStringFun(std::make_pair(discreteVectorMap[i], currentContinuous));
        std::string target = _vecToStringFun(std::make_pair(discreteVectorMap[n], currentContinuous));
        _edgeStream << "\"" << source << "\",\"" << target << "\"," << weight << std::endl;
      }
    }
  }
#endif
}

void GaussianClusterLogger::flush() {
#ifdef AUTOPAS_Log_GaussianCluster
  std::ofstream outputFile(_outputFileName, std::ofstream::out | std::ofstream::app);
  if (outputFile.is_open()) {
    outputFile << _nodeStream.str();
    outputFile << _edgeStream.str();
    outputFile << end_marker << std::endl;
  } else {
    AutoPasLog(error, "Output file {} is not open for writing!", _outputFileName);
  }
  outputFile.close();
  reset();
#endif
}
}  // namespace autopas
