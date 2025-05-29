/**
 * @file FuzzyTuning.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzyTuning.h"

#include <sys/stat.h>

#include <list>
#include <numeric>
#include <utility>

#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
#include "autopas/tuning/tuningStrategy/fuzzyTuning/OutputMapper.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/FuzzyRuleErrorListener.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/TranslationVisitor.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageLexer.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/parser_generated/autopas_generated_fuzzy_rule_syntax/FuzzyLanguageParser.h"
#include "autopas/utils/WrapMPI.h"
#endif

namespace autopas {

FuzzyTuning::FuzzyTuning(std::string fuzzyRuleFileName) : _fuzzyRuleFileName(std::move(fuzzyRuleFileName)) {
  // Check if the given rule file exists and throw if not
  struct stat buffer {};
  if (stat(_fuzzyRuleFileName.c_str(), &buffer) != 0) {
    utils::ExceptionHandler::exception(
        "Rule file '{}' does not exist! Please provide a valid rule file using the "
        "command line option '--fuzzy-rule-file'.",
        _fuzzyRuleFileName);
  }

#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING

  auto [fuzzyControlSettings, linguisticVariables, outputMappings, fuzzyControlSystems] = parse(_fuzzyRuleFileName);

  // By default, dump the rules for reproducibility reasons.
  std::string settingsStr = std::accumulate(
      fuzzyControlSettings->begin(), fuzzyControlSettings->end(), std::string("\nSettings:\n"),
      [](const std::string &acc, const auto &b) { return acc + "\t" + b.first + ": " + b.second + "\n"; });

  std::string linguisticVariablesStr =
      std::accumulate(linguisticVariables.begin(), linguisticVariables.end(), std::string("\n"),
                      [](const std::string &acc, const std::shared_ptr<FuzzyLogic::LinguisticVariable> &b) {
                        return acc + std::string(*b);
                      });

  std::string outputMappingsStr = std::accumulate(
      outputMappings.begin(), outputMappings.end(), std::string("\n"),
      [](const std::string &acc, const std::pair<const std::string, std::shared_ptr<FuzzyLogic::OutputMapper>> &b) {
        return acc + std::string(*b.second);
      });

  std::string fuzzyControlSystemsStr =
      std::accumulate(fuzzyControlSystems.begin(), fuzzyControlSystems.end(), std::string("\n"),
                      [](const std::string &acc,
                         const std::pair<const std::string, std::shared_ptr<FuzzyLogic::FuzzyControlSystem>> &b) {
                        return acc + std::string(*b.second);
                      });
  int myRank{};
  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  if (myRank == 0) {
    AutoPasLog(INFO, "FuzzyTuning initialized with the following Rules:");
    AutoPasLog(INFO, "{}", settingsStr);
    AutoPasLog(INFO, "{}", linguisticVariablesStr);
    AutoPasLog(INFO, "{}", outputMappingsStr);
    AutoPasLog(INFO, "{}", fuzzyControlSystemsStr);
  }
  _fuzzyControlSettings = fuzzyControlSettings;
  _fuzzyControlSystems = fuzzyControlSystems;
  _outputMappings = outputMappings;

#else
  autopas::utils::ExceptionHandler::exception(
      "FuzzyTuning constructed but AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=OFF! ");
#endif
}

bool FuzzyTuning::needsLiveInfo() const { return true; }

void FuzzyTuning::receiveLiveInfo(const LiveInfo &info) {
  // clear the current live info
  _currentLiveInfo.clear();

  // only filter out the numeric values
  for (const auto &infoEntry : info.get()) {
    const auto &name = infoEntry.first;
    const auto &entry = infoEntry.second;
    std::visit(
        [&](auto value) {
          if constexpr (std::is_arithmetic_v<decltype(value)>) {
            _currentLiveInfo[name] = value;
          }
        },
        entry);
  }
}

void FuzzyTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {}

bool FuzzyTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                        const EvidenceCollection &evidenceCollection) {
#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
  using namespace autopas::FuzzyLogic;

  // Runs the FuzzyTuning strategy for the first time. Since the decision of the Fuzzy Tuner doesn't depend on the
  // evidence, and it is expected that the LiveInfo doesn't change that much during a tuning phase. We can run the
  // FuzzyTuning strategy ONLY once at the beginning of a tuning phase.

  std::string outputInterpretation = _fuzzyControlSettings->at("interpretOutputAs");

  // Treat all fuzzy control systems as different systems and filter out configurations that are never predicted by any
  // of the systems
  if (outputInterpretation == "IndividualSystems") {
    updateQueueInterpretOutputAsIndividualSystems(configQueue);
  } else if (outputInterpretation == "Suitability") {
    updateQueueInterpretOutputAsSuitability(configQueue);
  } else {
    utils::ExceptionHandler::exception("FuzzyTuning: Unknown output interpretation: {}", outputInterpretation);
  }
#else
  autopas::utils::ExceptionHandler::exception(
      "FuzzyTuning constructed but AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=OFF! ");
#endif
  return false;
}

bool FuzzyTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                      const EvidenceCollection &evidenceCollection) {
  // This is left empty because the FuzzyTuning strategy is only run once at the beginning of a tuning phase.
  AutoPasLog(TRACE,
             "FuzzyTuning: optimizeSuggestions called. This method does nothing, as reset() already updated the "
             "config-queue");
  return false;
}

TuningStrategyOption FuzzyTuning::getOptionType() const { return TuningStrategyOption::fuzzyTuning; }

#ifdef AUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING
void FuzzyTuning::updateQueueInterpretOutputAsIndividualSystems(std::vector<Configuration> &configQueue) {
  using namespace autopas::FuzzyLogic;

  std::list<Configuration> newSearchSpace{configQueue.begin(), configQueue.end()};

  // For each fuzzy control system
  for (const auto &[name, fcs] : _fuzzyControlSystems) {
    // 1. Predict
    auto prediction = fcs->predict(_currentLiveInfo);

    // 2. Get closest configuration patterns
    const auto &configuration_patterns = _outputMappings.at(name)->getClosestConfigurationPatterns(prediction);

    AutoPasLog(DEBUG, "FuzzySystem '{}' predicted {}. Mapping to: {}", name, prediction, [&configuration_patterns]() {
      std::string result;
      for (const auto &pattern : configuration_patterns) {
        result += "[" + pattern.toString() + "] ";
      }
      return result;
    }());

    // 3. Filter out all configurations that do not match the patterns
    newSearchSpace.remove_if([&](const Configuration &configuration) {
      return std::none_of(configuration_patterns.begin(), configuration_patterns.end(),
                          [&configuration](const auto &pattern) { return pattern.matches(configuration); });
    });
  }

  AutoPasLog(DEBUG, "Fuzzy tuning selected {} out of {} possible configurations: {}", newSearchSpace.size(),
             configQueue.size(), [&newSearchSpace]() {
               std::string result;
               for (const auto &config : newSearchSpace) {
                 result += config.toString() + ", ";
               }
               return result;
             }());
  configQueue.clear();
  std::copy(newSearchSpace.rbegin(), newSearchSpace.rend(), std::back_inserter(configQueue));
}

void FuzzyTuning::updateQueueInterpretOutputAsSuitability(std::vector<Configuration> &configQueue) {
  using namespace autopas::FuzzyLogic;

  std::list<Configuration> newSearchSpace{configQueue.begin(), configQueue.end()};

  std::vector<std::pair<ConfigurationPattern, double>> configSuitabilities;

  // For each fuzzy control system
  for (const auto &[name, fcs] : _fuzzyControlSystems) {
    // 1. Predict
    auto prediction = fcs->predict(_currentLiveInfo);

    // 2. Get single configuration pattern for each system
    const auto &configurations = _outputMappings.at(name)->getClosestConfigurationPatterns(prediction);
    if (configurations.size() != 1) {
      utils::ExceptionHandler::exception(
          "FuzzyTuning: Expected exactly one configuration pattern for suitability "
          "value, but got {}. In system {} with prediction {}",
          configurations.size(), name, prediction);
    }

    // 3. Insert suitability values
    configSuitabilities.emplace_back(configurations[0], prediction);

    AutoPasLog(DEBUG, "FuzzySystem predicted suitability {} for pattern {}", prediction, configurations[0].toString());
  }

  // Sort the suitability values
  std::sort(configSuitabilities.begin(), configSuitabilities.end(),
            [](const auto &a, const auto &b) { return a.second > b.second; });

  // get the best configurations and everything that is within 10% of the best
  const double THRESHOLD = std::stod(_fuzzyControlSettings->at("suitabilityThreshold"));

  const auto bestSuitability = configSuitabilities.front().second;
  std::vector<ConfigurationPattern> bestConfigurations;
  for (const auto &[pattern, suitability] : configSuitabilities) {
    if (suitability < bestSuitability * (1 - THRESHOLD)) {
      break;
    }

    AutoPasLog(DEBUG,
               "FuzzyTuning: Adding configuration pattern {} with suitability {} since its within {}% of the "
               "best suitability {}",
               pattern.toString(), suitability, THRESHOLD * 100, bestSuitability);
    bestConfigurations.push_back(pattern);
  }

  // 4. Only keep configurations that are in the best configurations
  newSearchSpace.remove_if([&bestConfigurations](const Configuration &configuration) {
    return std::none_of(bestConfigurations.begin(), bestConfigurations.end(),
                        [&configuration](const auto &pattern) { return pattern.matches(configuration); });
  });

  AutoPasLog(DEBUG, "Fuzzy tuning selected {} out of {} possible configurations", newSearchSpace.size(),
             configQueue.size());

  configQueue.clear();
  std::copy(newSearchSpace.rbegin(), newSearchSpace.rend(), std::back_inserter(configQueue));
}

std::tuple<std::shared_ptr<FuzzyControlSettings>, std::vector<std::shared_ptr<LinguisticVariable>>,
           std::map<std::string, std::shared_ptr<OutputMapper>>,
           std::map<std::string, std::shared_ptr<FuzzyControlSystem>>>
FuzzyTuning::parse(const std::string &fuzzyRuleFilename) {
  using namespace antlr4;
  using namespace autopas::FuzzyLogic;

  try {
    // Read file
    std::ifstream stream(fuzzyRuleFilename);
    ANTLRInputStream input(stream);
    FuzzyLanguageLexer lexer(&input);

    // Lexer
    FuzzyRuleErrorListener ErrorListener;
    lexer.removeErrorListeners();
    lexer.addErrorListener(&ErrorListener);

    CommonTokenStream tokens(&lexer);

    // Parser
    FuzzyLanguageParser parser(&tokens);
    parser.removeErrorListeners();
    parser.addErrorListener(&ErrorListener);

    tree::ParseTree *tree = parser.rule_file();

    // Translation
    TranslationVisitor visitor;
    auto fuzzy_rule_program = std::any_cast<
        std::tuple<std::shared_ptr<FuzzyControlSettings>, std::vector<std::shared_ptr<LinguisticVariable>>,
                   std::map<std::string, std::shared_ptr<OutputMapper>>,
                   std::map<std::string, std::shared_ptr<FuzzyControlSystem>>>>(visitor.visit(tree));

    return fuzzy_rule_program;
  } catch (const std::exception &e) {
    throw utils::ExceptionHandler::AutoPasException(
        fmt::format("Error while parsing fuzzy rule file: \"{}\". Error: {}", fuzzyRuleFilename, e.what()));
  }
}

std::shared_ptr<FuzzyControlSettings> FuzzyTuning::getFuzzyControlSettings() const { return _fuzzyControlSettings; }

const std::map<std::string, std::shared_ptr<FuzzyControlSystem>> &FuzzyTuning::getFuzzyControlSystems() const {
  return _fuzzyControlSystems;
}

const std::map<std::string, std::shared_ptr<OutputMapper>> &FuzzyTuning::getOutputMappings() const {
  return _outputMappings;
}

#endif

}  // namespace autopas