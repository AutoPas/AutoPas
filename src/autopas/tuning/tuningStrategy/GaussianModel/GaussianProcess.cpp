/**
 * @file GaussianProcess.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "GaussianProcess.h"

autopas::GaussianProcess::GaussianProcess(size_t dims, double sigma, Random &rngRef)
    : _inputs(),
      _outputs(),
      _dims(dims),
      _evidenceMinValue(0),
      _evidenceMaxValue(0),
      _sigma(sigma),
      _hypers(),
      _rng(rngRef) {
  tuneHyperparameters();
}

autopas::GaussianProcess::~GaussianProcess() = default;

void autopas::GaussianProcess::setDimension(size_t dims) {
  _dims = dims;
  clear();
}

void autopas::GaussianProcess::clear() {
  _inputs.clear();
  _outputs = autopas::GaussianProcess::Vector::Zero(0);
  tuneHyperparameters();
}

size_t autopas::GaussianProcess::numEvidence() const { return _inputs.size(); }

std::pair<const std::vector<autopas::GaussianProcess::Vector> &, const autopas::GaussianProcess::Vector &>
autopas::GaussianProcess::getEvidence() const {
  return std::make_pair(std::cref(_inputs), std::cref(_outputs));
}

void autopas::GaussianProcess::addEvidence(const autopas::GaussianProcess::Vector &input, double output,
                                           bool tuneHypers) {
  if (static_cast<size_t>(input.size()) != _dims) {
    utils::ExceptionHandler::exception(
        "GaussianProcess.addEvidence: size of input {} does not match specified dimensions {}", input.size(), _dims);
  }

  if (_inputs.empty()) {
    // first evidence
    _evidenceMinValue = _evidenceMaxValue = output;
    _evidenceMaxVector = input;
  } else if (output < _evidenceMinValue) {
    _evidenceMinValue = output;
  } else if (output > _evidenceMaxValue) {
    _evidenceMaxValue = output;
    _evidenceMaxVector = input;
  }

  _inputs.push_back(input);
  long newSize = _inputs.size();

  // extend output vector
  _outputs.conservativeResize(newSize, Eigen::NoChange_t());
  _outputs(newSize - 1) = output;

  if (tuneHypers) {
    tuneHyperparameters();
  } else {
    // hyperparameters should be recalculated
    _hypers.clear();
  }
}

const autopas::GaussianProcess::Vector &autopas::GaussianProcess::getEvidenceMax() const {
  if (_inputs.empty()) {
    utils::ExceptionHandler::exception("GaussianProcess has no evidence");
  }

  return _evidenceMaxVector;
}

double autopas::GaussianProcess::predictMean(const autopas::GaussianProcess::Vector &input) const {
  if (static_cast<size_t>(input.size()) != _dims) {
    utils::ExceptionHandler::exception(
        "GaussianProcess.predictMean: size of input {} does not match specified dimensions {}", input.size(), _dims);
  }

  double result = 0.;
  if (_inputs.empty()) {
    // no evidence
    for (const auto &hyper : _hypers) {
      result += hyper.score * hyper.mean;
    }
  } else {
    for (const auto &hyper : _hypers) {
      result += hyper.score * (hyper.mean + kernelVector(input, hyper.theta, hyper.dimScales).dot(hyper.weights));
    }
  }

  return result;
}

double autopas::GaussianProcess::getDefaultVar() const {
  double result = 0.;
  for (const auto &hyper : _hypers) {
    result += hyper.score * hyper.theta;
  }
  return result;
}

double autopas::GaussianProcess::predictVar(const autopas::GaussianProcess::Vector &input) const {
  if (static_cast<size_t>(input.size()) != _dims) {
    utils::ExceptionHandler::exception(
        "GaussianProcess.predictVar: size of input {} does not match specified dimensions {}", input.size(), _dims);
  }

  double result = 0.;
  if (_inputs.empty()) {
    // no evidence
    return getDefaultVar();
  } else {
    for (const auto &hyper : _hypers) {
      auto kVec = kernelVector(input, hyper.theta, hyper.dimScales);
      result += hyper.score * (kernel(input, input, hyper.theta, hyper.dimScales) - kVec.dot(hyper.covMatInv * kVec));
    }
  }

  return result;
}

double autopas::GaussianProcess::predictOutputPDF(const autopas::GaussianProcess::Vector &input, double output) const {
  double stddev = std::sqrt(predictVar(input));
  double mean = predictMean(input);
  return utils::Math::normalPDF((mean - output) / stddev) / stddev;
}

double autopas::GaussianProcess::predictOutputScaledPDF(const autopas::GaussianProcess::Vector &input,
                                                        double output) const {
  double stddev = std::sqrt(predictVar(input));
  double mean = predictMean(input);
  return utils::Math::normalPDF((mean - output) / stddev) / utils::Math::normalScale;
}

double autopas::GaussianProcess::calcAcquisition(autopas::AcquisitionFunctionOption af,
                                                 const autopas::GaussianProcess::Vector &input) const {
  return AcquisitionFunction::calcAcquisition(af, predictMean(input), predictVar(input), _evidenceMaxValue);
}

autopas::GaussianProcess::Vector autopas::GaussianProcess::sampleAquisitionMax(
    autopas::AcquisitionFunctionOption af, const std::vector<autopas::GaussianProcess::Vector> &samples) const {
  size_t bestIdx = 0;
  double bestVal = calcAcquisition(af, samples[0]);

  // find optimum from samples
  for (size_t i = 1; i < samples.size(); ++i) {
    double val = calcAcquisition(af, samples[i]);

    if (val > bestVal) {
      bestIdx = i;
      bestVal = val;
    }
  }

  return samples[bestIdx];
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<autopas::GaussianProcess::Vector>>
autopas::GaussianProcess::generateHyperparameterSamples(size_t sampleSize, autopas::Random &rng, size_t dims,
                                                        double sigma, double evidenceMinValue,
                                                        double evidenceMaxValue) {
  // range of mean
  // inside bounds of evidence outputs
  NumberInterval<double> meanRange(evidenceMinValue, evidenceMaxValue);
  // range of theta
  // max sample stddev: (max - min)
  // max stddev from zero: abs(min) & abs(max)
  double thetaMax = std::pow(
      std::max({evidenceMaxValue - evidenceMinValue, std::abs(evidenceMinValue), std::abs(evidenceMaxValue)}), 2);
  // at least sigma
  thetaMax = std::max(thetaMax, sigma);
  NumberInterval<double> thetaRange(sigma, thetaMax);
  // range of dimScale
  // Assuming most distances are greater equal 1.
  // For a dimScale d > 5 + ln(thetaMax): theta * exp(-d r) < 1%.
  // So choosing a greater dimScale may lead to many kernels close to zero.
  // But if needed the upper bound can be increased.
  NumberInterval<double> dimScaleRange(0., 5. + std::max(0., std::log(thetaMax)));

  // generate mean
  auto sample_means = meanRange.uniformSample(sampleSize, rng);

  // generate theta
  auto sample_thetas = thetaRange.uniformSample(sampleSize, rng);

  // generate dimScale
  std::vector<std::vector<double>> sample_dimScaleData;
  sample_dimScaleData.reserve(dims);
  for (size_t d = 0; d < dims; ++d) {
    sample_dimScaleData.emplace_back(dimScaleRange.uniformSample(sampleSize, rng));
  }
  // convert dimScales to Vectors
  std::vector<autopas::GaussianProcess::Vector> sample_dimScales;
  sample_dimScales.reserve(sampleSize);
  for (size_t t = 0; t < sampleSize; ++t) {
    std::vector<double> dimScaleData;
    dimScaleData.reserve(dims);
    for (size_t d = 0; d < dims; ++d) {
      dimScaleData.push_back(sample_dimScaleData[d][t]);
    }
    sample_dimScales.emplace_back(utils::Math::makeVectorXd(dimScaleData));
  }

  return std::make_tuple(sample_means, sample_thetas, sample_dimScales);
}

std::vector<autopas::GaussianHyperparameters> &autopas::GaussianProcess::getHyperparameters() { return _hypers; }

void autopas::GaussianProcess::setHyperparameters(
    const std::vector<double> &sample_means, const std::vector<double> &sample_thetas,
    const std::vector<autopas::GaussianProcess::Vector> &sample_dimScales) {
  size_t hyperSize = sample_means.size();
  _hypers.clear();

  // initialize hyperparameter samples
  _hypers.reserve(hyperSize);
  for (size_t t = 0; t < hyperSize; ++t) {
    _hypers.emplace_back(sample_means[t], sample_thetas[t], sample_dimScales[t]);
  }

  // precalculate matrices for all hyperparameters
  // @TODO find sensible chunkSize
#ifdef AUTOPAS_OPENMP
  const size_t chunkSize = std::max(hyperSize / (autopas_get_num_threads() * 10), 1ul);
#pragma omp parallel for schedule(dynamic, chunkSize)
#endif
  for (size_t t = 0; t < hyperSize; ++t) {
    _hypers[t].precalculate(_sigma, _inputs, _outputs);
  }
}

void autopas::GaussianProcess::normalizeHyperparameters() {
  // sort by score
  std::sort(_hypers.begin(), _hypers.end(),
            [](const GaussianHyperparameters &h1, const GaussianHyperparameters &h2) { return h1.score > h2.score; });

  // only keep hp_size highest scores
  if (_hypers.size() > hp_size) {
    _hypers.erase(_hypers.begin() + hp_size, _hypers.end());
  }

  // normalize scores
  double scoreSum = 0.;
  for (auto &hyper : _hypers) {
    scoreSum += hyper.score;
  }
  if (scoreSum > 0) {
    for (auto &hyper : _hypers) {
      hyper.score /= scoreSum;
    }
  } else {
    // all scores are 0
    double uniformProbability = 1. / _hypers.size();
    for (auto &hyper : _hypers) {
      hyper.score = uniformProbability;
    }
  }
}

void autopas::GaussianProcess::tuneHyperparameters() {
  // number of evidence
  size_t newSize = _inputs.size();
  _hypers.clear();

  // if no evidence
  if (newSize == 0) {
    // use default values
    _hypers.emplace_back(0., 1., autopas::GaussianProcess::Vector::Ones(_dims));
    _hypers[0].precalculate(_sigma, _inputs, _outputs);
    _hypers[0].score = 1.;
    return;
  }

  auto [sample_means, sample_thetas, sample_dimScales] =
      generateHyperparameterSamples(hp_sample_size, _rng, _dims, _sigma, _evidenceMinValue, _evidenceMaxValue);
  setHyperparameters(sample_means, sample_thetas, sample_dimScales);
  normalizeHyperparameters();
}

double autopas::GaussianProcess::kernel(const autopas::GaussianProcess::Vector &input1,
                                        const autopas::GaussianProcess::Vector &input2, double theta,
                                        const autopas::GaussianProcess::Vector &dimScale) {
  double dot = 0;
  for (int i = 0; i < input1.size(); ++i) {
    double dist = input1[i] - input2[i];
    dist *= dist * dimScale[i];
    dot += dist;
  }
  return theta * std::exp(-dot);
}

autopas::GaussianProcess::Vector autopas::GaussianProcess::kernelVector(
    const autopas::GaussianProcess::Vector &input, double theta,
    const autopas::GaussianProcess::Vector &dimScale) const {
  std::vector<double> k(_inputs.size());
  for (size_t i = 0; i < k.size(); ++i) {
    k[i] = kernel(input, _inputs[i], theta, dimScale);
  }
  return utils::Math::makeVectorXd(k);
}
