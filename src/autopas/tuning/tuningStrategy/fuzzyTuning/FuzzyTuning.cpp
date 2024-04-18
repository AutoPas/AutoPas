/**
 * @file FuzzyTuning.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzyTuning.h"

#include "CrispSet.h"
#include "FuzzyRule.h"
#include "FuzzySystem.h"
#include "LinguisticVariable.h"
#include "MembershipFunction.h"

namespace autopas {

FuzzyTuning::FuzzyTuning() {}

bool FuzzyTuning::needsLiveInfo() const { return true; }

void FuzzyTuning::receiveLiveInfo(const LiveInfo &info) { _currentLiveInfo = info; }

void FuzzyTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  using namespace fuzzy_logic;

  auto service_crisp = std::make_shared<CrispSet>(CrispSet({{"service", std::pair(0, 10)}}));
  auto food_crisp = std::make_shared<CrispSet>(CrispSet({{"food", std::pair(0, 10)}}));
  auto tip_crisp = std::make_shared<CrispSet>(CrispSet({{"tip", std::pair(0, 30)}}));

  auto service = LinguisticVariable(service_crisp);
  auto food = LinguisticVariable(food_crisp);
  auto tip = LinguisticVariable(tip_crisp);

  auto s1 = MembershipFunction::makeGaussian("poor", 0, 1.5);
  auto s2 = MembershipFunction::makeGaussian("good", 5, 1.5);
  auto s3 = MembershipFunction::makeGaussian("excellent", 10, 1.5);
  auto f1 = MembershipFunction::makeTrapezoid("rancid", 0, 0, 1, 3);
  auto f2 = MembershipFunction::makeTrapezoid("delicious", 7, 9, 10, 10);
  auto f3 = MembershipFunction::makeGaussian("delicious", 10, 1.5);
  auto t1 = MembershipFunction::makeTriangle("cheap", 0, 5, 10);
  auto t2 = MembershipFunction::makeTriangle("average", 10, 15, 20);
  auto t3 = MembershipFunction::makeTriangle("generous", 20, 25, 30);

  service.addLinguisticTerm(std::make_shared<FuzzySet>(s1));
  service.addLinguisticTerm(std::make_shared<FuzzySet>(s2));
  service.addLinguisticTerm(std::make_shared<FuzzySet>(s3));

  food.addLinguisticTerm(std::make_shared<FuzzySet>(f1));
  food.addLinguisticTerm(std::make_shared<FuzzySet>(f2));
  food.addLinguisticTerm(std::make_shared<FuzzySet>(f3));

  tip.addLinguisticTerm(std::make_shared<FuzzySet>(t1));
  tip.addLinguisticTerm(std::make_shared<FuzzySet>(t2));
  tip.addLinguisticTerm(std::make_shared<FuzzySet>(t3));

  auto fs1 = (service == "poor") || (food == "rancid");
  auto r1 = FuzzyRule(fs1, tip == "cheap");

  auto fs2 = (service == "good");
  auto r2 = FuzzyRule(fs2, tip == "average");

  auto fs3 = (service == "excellent") || (food == "delicious");
  auto r3 = FuzzyRule(fs3, tip == "generous");

  auto fs = FuzzySystem();
  fs.addRule(r1);
  fs.addRule(r2);
  fs.addRule(r3);

  const std::map<std::string, double> &configMap = {{"service", 1}, {"food", 3}};

  auto x2 = fs.predict(configMap, 100);

  std::cout << "Predicted tip: " << x2 << std::endl;
}

void FuzzyTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                        const EvidenceCollection &evidenceCollection) {}

void FuzzyTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                      const EvidenceCollection &evidenceCollection) {}

TuningStrategyOption FuzzyTuning::getOptionType() { return TuningStrategyOption::fuzzyTuning; }

}  // namespace autopas