/**
 * @file GaussianClusterLogger.cpp
 * @author Jan Nguyen
 * @date 11.08.20
 */

#include "GaussianClusterLogger.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "autopas/utils/Logger.h"

namespace autopas {

GaussianClusterLogger::GaussianClusterLogger(GaussianModelTypes::VectorToStringFun vecToStringFun,
                                             GaussianClusterLogger::OutputType outType)
    : _outType(outType),
      _nodeStream(std::ios_base::out | std::ios_base::app),
      _edgeStream(std::ios_base::out | std::ios_base::app),
      _vecToStringFun(std::move(vecToStringFun)) {
  if (generatesNoOutput()) {
    return;
  }

  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::ostringstream nowStrStr;
  tm unused;
  nowStrStr << std::put_time(localtime_r(&now, &unused), "%Y-%m-%d_%H-%M-%S");
  _outputFileName = "gaussian_graph_" + nowStrStr.str() + ".out";

  reset();
}

void GaussianClusterLogger::setVectorToStringFun(const GaussianModelTypes::VectorToStringFun &fun) {
  _vecToStringFun = fun;
}

bool GaussianClusterLogger::generatesNoOutput() const {
  return (_outType == OutputType::none or
          (_outType == OutputType::trace and autopas::Logger::get()->level() > autopas::Logger::LogLevel::trace));
}

void GaussianClusterLogger::reset() {
  _nodeStream.str(node_start_marker);
  _nodeStream << std::endl << "Label,mean,var,numEvidence" << std::endl;
  _edgeStream.str(edge_start_marker);
  _edgeStream << std::endl << "Source,Target,Weight" << std::endl;

  _currentContinuous.clear();
}

void GaussianClusterLogger::add(const std::vector<GaussianProcess> &clusters,
                                const std::vector<GaussianModelTypes::VectorDiscrete> &discreteVectorMap,
                                const GaussianModelTypes::VectorContinuous &currentContinuous,
                                std::vector<double> means, std::vector<double> vars,
                                const GaussianModelTypes::NeighboursWeights &neighbourWeights) {
  if (generatesNoOutput()) {
    return;
  }

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
}

void GaussianClusterLogger::end() {
  if (generatesNoOutput()) {
    return;
  }

  switch (_outType) {
    case OutputType::trace:
      AutoPasLog(trace, _nodeStream.str());
      AutoPasLog(trace, _edgeStream.str());
      AutoPasLog(trace, end_marker);
      break;
    case OutputType::file: {
      std::ofstream outputFile(_outputFileName, std::ofstream::out | std::ofstream::app);
      if (outputFile.is_open()) {
        outputFile << _nodeStream.str();
        outputFile << _edgeStream.str();
        outputFile << end_marker << std::endl;
      }
      outputFile.close();
    } break;
    case OutputType::none:
      utils::ExceptionHandler::exception("GaussianClusterLogger.end: Unexpected output type {}", _outType);
      break;
  }
  reset();
}
}  // namespace autopas
