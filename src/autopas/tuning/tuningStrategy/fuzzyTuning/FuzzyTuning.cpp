/**
 * @file FuzzyTuning.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzyTuning.h"

#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyControlSystem.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzyRule.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/FuzzySetFactory.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/fuzzyController/LinguisticVariable.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/xmlParser/XMLParser.h"

namespace autopas {

FuzzyTuning::FuzzyTuning(const std::string &ruleFileName) : _ruleFileName(ruleFileName) {}

bool FuzzyTuning::needsLiveInfo() const { return true; }

void FuzzyTuning::receiveLiveInfo(const LiveInfo &info) { _currentLiveInfo = info; }

void FuzzyTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  using namespace autopas::fuzzy_logic;

  parse(_ruleFileName);

  // Just a demo showing how to use the fuzzy logic controller

  auto service = LinguisticVariable("service", std::pair(0, 10));
  auto food = LinguisticVariable("food", std::pair(0, 10));
  auto tip = LinguisticVariable("tip", std::pair(0, 30));

  service.addLinguisticTerm(makeGaussian("poor", 0, 1.5));
  service.addLinguisticTerm(makeGaussian("good", 5, 1.5));
  service.addLinguisticTerm(makeGaussian("excellent", 10, 1.5));

  food.addLinguisticTerm(makeTrapezoid("rancid", 0, 0, 1, 3));
  food.addLinguisticTerm(makeTrapezoid("delicious", 7, 9, 10, 10));

  tip.addLinguisticTerm(makeTriangle("cheap", 0, 5, 10));
  tip.addLinguisticTerm(makeTriangle("average", 10, 15, 20));
  tip.addLinguisticTerm(makeTriangle("generous", 20, 25, 30));

  auto fcs = FuzzyControlSystem();
  fcs.addRule(FuzzyRule(service == "poor" || food == "rancid", tip == "cheap"));
  fcs.addRule(FuzzyRule(service == "good", tip == "average"));
  fcs.addRule(FuzzyRule(service == "excellent" || food == "delicious", tip == "generous"));

  // Make a prediction for a sample input
  FuzzySet::Data sample = {{"service", 1}, {"food", 3}};
  auto x2 = fcs.predict(sample);

  std::cout << "Predicted tip: " << x2 << std::endl;
}

void FuzzyTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                        const EvidenceCollection &evidenceCollection) {}

void FuzzyTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                      const EvidenceCollection &evidenceCollection) {}

TuningStrategyOption FuzzyTuning::getOptionType() { return TuningStrategyOption::fuzzyTuning; }

}  // namespace autopas