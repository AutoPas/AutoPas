/**
 * @file MPIParallelizedStrategy.cpp
 * @author W. Thieme
 * @date 17.09.2020
 */

#include "MPIParallelizedStrategy.h"

#include <cstddef>

#include "autopas/tuning/Configuration.h"
#include "options/DataLayoutOption.h"
#include "options/Newton3Option.h"
#include "utils/AutoPasConfigurationCommunicator.h"
#include "utils/WrapMPI.h"

namespace autopas {

MPIParallelizedStrategy::MPIParallelizedStrategy(const Configuration &fallbackConfiguration,
                                                 const AutoPas_MPI_Comm &comm, double mpiTuningMaxDifferenceForBucket,
                                                 double mpiTuningWeightForMaxDensity)
    : _comm(comm),
      _fallbackConfiguration(fallbackConfiguration),
      _mpiTuningMaxDifferenceForBucket(mpiTuningMaxDifferenceForBucket),
      _mpiTuningWeightForMaxDensity(mpiTuningWeightForMaxDensity) {}

bool MPIParallelizedStrategy::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                  const EvidenceCollection &evidenceCollection) {
  // sanity check
  if (_bucket == AUTOPAS_MPI_COMM_NULL) {
    AutoPasLog(WARN, "_bucket was AUTOPAS_MPI_COMM_NULL");
    _bucket = _comm;
  }

  // All ranks should stay in tuning mode equally long so that none settles on an optimum
  // before the other's data is there.
  const auto [myBestConf, myBestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  const auto globallyBestConfig =
      utils::AutoPasConfigurationCommunicator::findGloballyBestConfiguration(_bucket, myBestConf, myBestEvidence.value);

  const auto myQueueSize = static_cast<unsigned int>(configQueue.size());
  unsigned int globallyLongestQueueSize{};
  AutoPas_MPI_Allreduce(&myQueueSize, &globallyLongestQueueSize, 1, AUTOPAS_MPI_UNSIGNED_INT, AUTOPAS_MPI_MAX, _bucket);
  // if our queue is about to be emptied but others still have things to test, stretch our queue by inserting the
  // current global optimum
  if (myQueueSize == 1 and globallyLongestQueueSize > myQueueSize) {
    // depending on if the global optimum is applicable to our subdomain insert this or the fallback.
    if (_rejectedConfigurations.count(globallyBestConfig) == 0) {
      // insert at the front so that the remaining config is at the back and picked next.
      configQueue.insert(configQueue.begin(), globallyBestConfig);
    } else {
      _usingFallbackConfig = true;
      configQueue.insert(configQueue.begin(), _fallbackConfiguration);
    }
  }
  // make sure not everyone is in fallback mode.
  bool everyoneUsingFallback{false};
  AutoPas_MPI_Reduce(&_usingFallbackConfig, &everyoneUsingFallback, 1, AUTOPAS_MPI_CXX_BOOL, AUTOPAS_MPI_LAND, 0,
                     _bucket);
  if (everyoneUsingFallback) {
    utils::ExceptionHandler::exception(
        "MPIParallelizedStrategy::optimizeSuggestions: All ranks defaulted to their fallback configuration. This means "
        "no rank could find an applicable configuration in their respective configuration sub queue");
  }
  // MPIParallelizedStrategy does no intentional config wipes to stop the tuning phase
  return false;
}

Configuration MPIParallelizedStrategy::createFallBackConfiguration(const std::set<Configuration> &searchSpace,
                                                                   const InteractionTypeOption &interactionType) {
  Configuration fallBackConfig{ContainerOption::linkedCells,
                               1.,
                               TraversalOption::lc_c08,
                               LoadEstimatorOption::none,
                               DataLayoutOption::aos,
                               Newton3Option::disabled,
                               interactionType};

  if (interactionType == InteractionTypeOption::triwise) {
    fallBackConfig.traversal = TraversalOption::lc_c01;
  }

  // Go through the search space and see if SoA or N3 are allowed.
  bool foundSoA{false};
  bool foundN3Enabled{false};
  for (const auto &conf : searchSpace) {
    if (not foundSoA and conf.dataLayout == DataLayoutOption::soa) {
      fallBackConfig.dataLayout = DataLayoutOption::soa;
      foundSoA = true;
    }
    if (not foundN3Enabled and conf.newton3 == Newton3Option::enabled) {
      fallBackConfig.newton3 = Newton3Option::enabled;
      foundN3Enabled = true;
    }
    if (foundSoA and foundN3Enabled) {
      break;
    }
  }
  return fallBackConfig;
}

void MPIParallelizedStrategy::rejectConfiguration(const Configuration &configuration, bool indefinitely) {
  _rejectedConfigurations.insert(configuration);
}

bool MPIParallelizedStrategy::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                    const EvidenceCollection &evidenceCollection) {
  // clear rejects since they now might be valid.
  _rejectedConfigurations.clear();

  // Rebuild the grouping (_bucket) of ranks according to subdomain similarity.
  autopas::utils::AutoPasConfigurationCommunicator::distributeRanksInBuckets(
      _comm, &_bucket, _smoothedPDBinStdDevDensity, _smoothedPDBinMaxDensity, _mpiTuningMaxDifferenceForBucket,
      _mpiTuningWeightForMaxDensity);

  // Grab our subset of configurations from the queue
  int myBucketRank{};
  AutoPas_MPI_Comm_rank(_bucket, &myBucketRank);
  int numRanksInBucket{};
  AutoPas_MPI_Comm_size(_bucket, &numRanksInBucket);
  std::vector<Configuration> myConfigQueue{};
  // Check if there are configs left for this rank.
  if (static_cast<size_t>(myBucketRank) < configQueue.size()) {
    // If there are enough configurations, distribute the configs in a round-robin fashion.
    // This way all ranks get a as similar distribution of containers as possible.

    // Equivalent to ceil(configQueue / numRanksInBucket) just purely in integers.
    const auto configsPerRankCeiled = (configQueue.size() + numRanksInBucket - 1) / numRanksInBucket;
    myConfigQueue.reserve(configsPerRankCeiled);
    for (size_t i = myBucketRank; i < configQueue.size(); i += numRanksInBucket) {
      myConfigQueue.push_back(configQueue[i]);
    }
  } else {
    // If there are not enough configurations pick one randomly.
    std::uniform_int_distribution<decltype(_rng)::result_type> distribution(0, configQueue.size() - 1);
    const auto randomConfigId = distribution(_rng);
    myConfigQueue.push_back(configQueue[randomConfigId]);
  }
  configQueue = myConfigQueue;

  // MPIParallelizedStrategy does no intentional config wipes to stop the tuning phase
  return false;
}

bool MPIParallelizedStrategy::needsDomainSimilarityStatistics() const { return true; }

void MPIParallelizedStrategy::receiveDomainSimilarityStatistics(double pdBinStdDevDensity, double pdBinMaxDensity) {
  _smoothedPDBinStdDevDensity = pdBinStdDevDensity;
  _smoothedPDBinMaxDensity = pdBinMaxDensity;
}

const AutoPas_MPI_Comm &MPIParallelizedStrategy::getBucket() const { return _bucket; }

TuningStrategyOption MPIParallelizedStrategy::getOptionType() const {
  return TuningStrategyOption::mpiDivideAndConquer;
}
}  // namespace autopas