#include "autopas/selectors/tuningStrategy/TuningStrategyLoggerProxy.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyFactory.h"

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
    auto strategy =
        autopas::TuningStrategyFactory::generateTuningStrategy(
            autopas::TuningStrategyOption::fullSearch, containers,
            *std::make_unique<autopas::NumberSetFinite<double>>(std::set<double>({1.})),
            traversals, loadEstimators, dataLayouts, newton3Options, 10, 1.2,
            5, 0, 3,
            autopas::AcquisitionFunctionOption::upperConfidenceBound,
            autopas::ExtrapolationMethodOption::linearRegression, "", autopas::MPIStrategyOption::noMPI,
            autopas::AutoPas_MPI_Comm::AUTOPAS_MPI_COMM_NULL);
    autopas::TuningStrategyLogReplayer logReplayer{argv[i], std::move(strategy)};
    logReplayer.replay();
  }


  autopas::Logger::unregister();
}