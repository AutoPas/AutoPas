/**
 * @file RuleBasedTuning.cpp
 * @author F. Gratl
 * @date 28.07.23
 */

#include "RuleBasedTuning.h"

#include <sys/stat.h>

namespace autopas {

RuleBasedTuning::RuleBasedTuning(const std::set<Configuration> &searchSpace, bool verifyModeEnabled,
                                 std::string ruleFileName, RuleBasedTuning::PrintTuningErrorFunType tuningErrorPrinter)
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
    : _verifyModeEnabled(verifyModeEnabled),
      _tuningErrorPrinter(std::move(tuningErrorPrinter)),
      _searchSpace(searchSpace),
      _originalSearchSpace(searchSpace.begin(), searchSpace.end()),
      _ruleFileName(std::move(ruleFileName)) {
  // Check if the given rule file exists and throw if not
  struct stat buffer {};
  if (stat(_ruleFileName.c_str(), &buffer) != 0) {
    utils::ExceptionHandler::exception("Rule file {} does not exist!", _ruleFileName);
  }
  // By default, dump the rules for reproducibility reasons.
  int myRank{};
  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    AutoPasLog(INFO, "Rule File {}:\n{}", _ruleFileName, rulesToString(_ruleFileName));
  }
#else
{
  autopas::utils::ExceptionHandler::exception(
      "RuleBasedTuning constructed but AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=OFF! ");
#endif
}

bool RuleBasedTuning::needsLiveInfo() const { return true; }

void RuleBasedTuning::receiveLiveInfo(const LiveInfo &info) { _currentLiveInfo = info; }

void RuleBasedTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  _tuningTime += evidence.value;
  _tuningTimeLifetime += evidence.value;
  _traversalTimes[configuration] = evidence.value;
  if (_verifyModeEnabled) {
    verifyCurrentConfigTime(configuration);
    if (_removedConfigurations.find(configuration) != _removedConfigurations.end()) {
      _wouldHaveSkippedTuningTime += evidence.value;
      _wouldHaveSkippedTuningTimeLifetime += evidence.value;
    }
  }
#endif
}

bool RuleBasedTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                            const EvidenceCollection &evidenceCollection) {
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  if (_verifyModeEnabled and _tuningTime > 0) {
    // It is fine to leave this log statement on Info level because it is behind if (_verificationModeEnabled).
    AutoPasLog(INFO, "Rules would have saved {} ns removing {}/{} configurations. ({}% of total tuning time)",
               _wouldHaveSkippedTuningTime, _removedConfigurations.size(), _originalSearchSpace.size(),
               // TODO: This lambda could probably be replaced by some formatting parameters in the fmt string.
               [&]() {
                 const auto percent =
                     ((static_cast<double>(_wouldHaveSkippedTuningTime) / static_cast<double>(_tuningTime)) * 100);
                 const auto percentRounded = std::round(percent * 100) / 100;
                 return percentRounded;
               }());
  }

  _tuningTime = 0;
  _wouldHaveSkippedTuningTime = 0;
  _removedConfigurations.clear();
  _rulesTooHarsh = false;
  optimizeSuggestions(configQueue, evidenceCollection);
#endif

  // RuleBasedTuning does no intentional config wipes to stop the tuning phase
  return false;
}

long RuleBasedTuning::getLifetimeWouldHaveSkippedTuningTime() const { return _wouldHaveSkippedTuningTimeLifetime; }

long RuleBasedTuning::getLifetimeTuningTime() const { return _tuningTimeLifetime; }

bool RuleBasedTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                          const EvidenceCollection &evidenceCollection) {
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  _lastApplicableConfigurationOrders = applyRules(configQueue);

  // Don't apply rules if they would wipe the queue and nothing has been tested yet.
  if (_rulesTooHarsh or (_searchSpace.empty() and _tuningTime == 0)) {
    _rulesTooHarsh = true;
    AutoPasLog(WARN, "Rules would remove all available options! Not applying them until next reset.");
    return false;
  }

  if (not _verifyModeEnabled) {
    configQueue.clear();
    std::copy(_searchSpace.rbegin(), _searchSpace.rend(), std::back_inserter(configQueue));
  }
#endif
  // RuleBasedTuning does no intentional config wipes to stop the tuning phase
  return false;
}

TuningStrategyOption RuleBasedTuning::getOptionType() const { return TuningStrategyOption::ruleBasedTuning; }

#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING

std::string RuleBasedTuning::rulesToString(const std::string &filePath) const {
  std::ifstream ruleFile(filePath);
  std::string line;
  std::stringstream ruleFileStringStream;
  while (std::getline(ruleFile, line)) {
    // any string that has a '#' only preceded by space characters.
    const std::regex rgxContainsComment("^[ \t]*#.*");

    if (line.empty() or std::regex_match(line, rgxContainsComment)) {
      continue;
    }
    ruleFileStringStream << line << "\n";
  }
  return ruleFileStringStream.str();
}

void RuleBasedTuning::verifyCurrentConfigTime(const Configuration &configuration) const {
  for (const auto &order : _lastApplicableConfigurationOrders) {
    bool shouldBeBetter{};
    if (order.smaller.matches(configuration)) {
      shouldBeBetter = true;
    } else if (order.greater.matches(configuration)) {
      shouldBeBetter = false;
    } else {
      continue;
    }

    const auto currentConfigTime = _traversalTimes.at(configuration);
    const auto &comparePattern = shouldBeBetter ? order.greater : order.smaller;
    for (const auto &[otherConfig, time] : _traversalTimes) {
      bool error = false;
      if (comparePattern.matches(otherConfig) and order.haveEqualSameProperties(configuration, otherConfig)) {
        if ((shouldBeBetter and time > currentConfigTime) or (not shouldBeBetter and time < currentConfigTime)) {
          error = true;
        }
      }
      if (error) {
        if (shouldBeBetter) {
          _tuningErrorPrinter(order, configuration, currentConfigTime, otherConfig, time, _currentLiveInfo);
        } else {
          _tuningErrorPrinter(order, otherConfig, time, configuration, currentConfigTime, _currentLiveInfo);
        }
      }
    }
  }
}

std::vector<RuleSyntax::ConfigurationOrder> RuleBasedTuning::applyRules(const std::vector<Configuration> &searchSpace) {
  AutoPasLog(DEBUG, _currentLiveInfo.toString());

  std::vector<RuleVM::MemoryCell> initialStack;
  std::vector<std::pair<std::string, RuleSyntax::Define>> defines{};
  for (const auto &[name, value] : _currentLiveInfo.get()) {
    initialStack.emplace_back(value);
    defines.push_back({name, {name, std::make_shared<RuleSyntax::Literal>(value)}});
  }

  std::ifstream ifs{_ruleFileName};
  const std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));

  RuleSyntax::RuleBasedProgramParser parser{defines};
  auto [programTree, context] = parser.parse(content);

  const auto generatedProgram = programTree.generateCode(context);

  RuleVM vm{};
  const auto removePatterns = vm.execute(generatedProgram, initialStack);

  AutoPasLog(DEBUG, "Remove patterns (Count {}):", removePatterns.size());
  std::vector<ConfigurationPattern> toRemovePatterns{};
  std::vector<RuleSyntax::ConfigurationOrder> applicableConfigurationOrders{};
  for (const auto &patternIdx : removePatterns) {
    const auto pattern = context.smallerConfigurationPatternByIndex(patternIdx);
    toRemovePatterns.push_back(pattern);
    AutoPasLog(DEBUG, "Remove {}", pattern.toString());

    applicableConfigurationOrders.push_back(context.getConfigurationOrders().at(patternIdx));
  }

  std::list<Configuration> newSearchSpace{searchSpace.begin(), searchSpace.end()};
  newSearchSpace.remove_if([&](const Configuration &configuration) {
    const bool remove = std::any_of(toRemovePatterns.begin(), toRemovePatterns.end(),
                                    [&configuration](const auto &pattern) { return pattern.matches(configuration); });
    if (remove) {
      _removedConfigurations.insert(configuration);
    }
    return remove;
  });
  // log number of removed configurations in this call to applyRules.
  AutoPasLog(DEBUG, "Rules remove {} out of {} configurations", searchSpace.size() - newSearchSpace.size(),
             searchSpace.size());
  if (not _verifyModeEnabled) {
    _searchSpace.clear();
    _searchSpace.insert(newSearchSpace.begin(), newSearchSpace.end());
  }

  return applicableConfigurationOrders;
}

#endif

}  // namespace autopas
