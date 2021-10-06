#include "autopas/selectors/tuningStrategy/TuningStrategyLoggerProxy.h"
#include "autopas/selectors/tuningStrategy/ruleBasedTuning/RuleBasedTuning.h"

int main(int argc, char** argv) {
  autopas::Logger::create();

  if(argc <= 1) {
    std::cerr << "Please provide the data files as arguments" << std::endl;
    exit(-1);
  }

  if(not getenv("DISABLE_DEBUG_LOG")) {
    autopas::Logger::get()->set_level(spdlog::level::debug);
  }

  auto containers = autopas::ContainerOption::getAllOptions();
  auto traversals = autopas::TraversalOption::getAllOptions();
  auto loadEstimators = autopas::LoadEstimatorOption::getAllOptions();
  auto dataLayouts = autopas::DataLayoutOption::getAllOptions();
  auto newton3Options = autopas::Newton3Option::getAllOptions();

  for(int i = 1; i < argc; i++) {
    AutoPasLog(info, "Checking with file {}", argv[i]);
    auto strategy = std::make_unique<autopas::RuleBasedTuning>(
        containers, std::set<double>({1.}), traversals, loadEstimators, dataLayouts, newton3Options,
        true);
    autopas::TuningStrategyLogReplayer logReplayer{argv[i], std::move(strategy)};
    logReplayer.replay();
  }


  autopas::Logger::unregister();
}