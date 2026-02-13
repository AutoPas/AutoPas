#include <memory>

#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyLogReplayer.h"
#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedTuning.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/NumberSetFinite.h"

/**
 * The program analyzes the log of a tuning phase with respect to a given rule file. The following aspects are checked:
 *   - if the given log contradicts any rules from the file.
 *   - how many configs would have been thrown out by rules.
 *   - how much time could have been saved.
 *   - what the actual best configuration was.
 *
 * @param argc >=2
 * @param argv tuningRules.rule tuningLog.txt...
 * @return 0 on success
 */
int main(int argc, char **argv) {
  autopas::Logger logger{};

  if (argc <= 2) {
    std::cerr
        << "Usage: ruleChecker tuningRules.rule tuningLog.txt...\n"
           "Where:\n"
           " - tuningRules.rule is a file containing rules for the rule based tuning.\n"
           " - tuningLog.txt is one or many log files created via the `--use-tuning-log true` option to md-flexible."
        << std::endl;
    exit(-1);
  }

  if (not getenv("DISABLE_DEBUG_LOG")) {
    autopas::Logger::get()->set_level(spdlog::level::info);
  }

  // Map configuration to indices of files where they were best
  std::map<autopas::Configuration, std::vector<int>> bestConfigs;

  int numErrors = 0;
  constexpr double bigErrorThreshold = 1.15;
  int numBigErrors = 0;
  auto errorHandler = [&](const autopas::RuleSyntax::ConfigurationOrder &order,
                          const autopas::Configuration &actualBetterConfig, unsigned long betterRuntime,
                          const autopas::Configuration &shouldBeBetterConfig, unsigned long shouldBeBetterRuntime,
                          const autopas::LiveInfo &liveInfo) {
    numErrors++;

    auto factorDifference = static_cast<double>(shouldBeBetterRuntime) / static_cast<double>(betterRuntime);
    auto factorDifferenceRounded = std::round(factorDifference * 100) / 100;

    if (factorDifference >= bigErrorThreshold) {
      numBigErrors++;
    }

    AutoPasLog(ERROR,
               "\n"
               "\tError in ConfigurationOrder {}:\n"
               "\t\t{}ns for config\t{}\n"
               "\t\t{}ns for config\t{}\n"
               "\t\tx{} difference",
               order.toString(), betterRuntime, actualBetterConfig.toShortString(), shouldBeBetterRuntime,
               shouldBeBetterConfig.toShortString(), factorDifferenceRounded);
  };

  unsigned long tuningTimeSum = 0;
  unsigned long wouldHaveSkippedTuningTimeSum = 0;
  const std::string rulesfile{argv[1]};
  const autopas::NumberSetFinite<double> csfs({1., 2.});
  // @TODO: Have rules for triwise interactions
  const std::set<autopas::Configuration> searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      autopas::ContainerOption::getAllOptions(), autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), autopas::DataLayoutOption::getAllOptions(),
      autopas::Newton3Option::getAllOptions(), &csfs, autopas::InteractionTypeOption::pairwise);

  for (int i = 2; i < argc; i++) {
    const std::string filename{argv[i]};
    // -1 because the first file is the rules file
    AutoPasLog(INFO, "Checking file {}: {}", i - 1, filename);
    auto strategy =
        std::make_shared<autopas::RuleBasedTuning>(searchSpace, /*verification mode*/ true, rulesfile, errorHandler);
    autopas::TuningStrategyLogReplayer logReplayer{filename, strategy, searchSpace};
    auto optBestConfig = logReplayer.replay();
    AutoPasLog(INFO, "");

    if (optBestConfig.has_value()) {
      bestConfigs[optBestConfig.value()].push_back(i);
    }
    tuningTimeSum += strategy->getLifetimeTuningTime();
    wouldHaveSkippedTuningTimeSum += strategy->getLifetimeWouldHaveSkippedTuningTime();
  }

  std::stringstream str;
  for (const auto &[config, filenames] : bestConfigs) {
    str << "\t"
        << "Best in " << filenames.size() << " scenarios:\t" << config.toShortString() << " (file numbers: ";
    for (const auto &file : filenames) {
      str << file << ", ";
    }
    str << ")" << std::endl;
  }

  // -2 because 0 is exe 1 is rules.
  AutoPasLog(INFO, "Finished replaying {} scenarios!", argc - 2);
  AutoPasLog(INFO, "\nSummary of best configurations:\n{}", str.str());
  AutoPasLog(INFO, "In sum, found {} errors! Of these, {} errors where greater than {}", numErrors, numBigErrors,
             bigErrorThreshold);

  auto savedTuningTimeRatio =
      static_cast<double>(wouldHaveSkippedTuningTimeSum) / static_cast<double>(tuningTimeSum) * 100;
  auto savedTuningTimeRatioRounded = std::round(savedTuningTimeRatio * 100) / 100;
  AutoPasLog(INFO, "Overall, {}% of the tuning time would have been saved.", savedTuningTimeRatioRounded);
}