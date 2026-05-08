/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 30.04.2020
 */

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "WrapMPI.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

/**
 * Provides several functions for handling configurations among mpi ranks.
 * This includes functionality for (de)serialization of configurations, splitting up search spaces based on ranks,
 * and finding the globally optimal configuration given time measurements.
 */

namespace autopas::utils::AutoPasConfigurationCommunicator {

/**
 * type definition for the serialization of configurations. A serialized config is an array of 16 bytes:
 * - container: 1 byte (<= less than 256 containers)
 * - traversal: 1 byte (<= less than 256 traversals)
 * - load estimator: 1 byte (<= less than 256 load estimators)
 * - data layout: 1 byte (<= less than 256 data layouts)
 * - newton3: 1 byte (<= less than 256 newton3 options)
 * - ompKind: 1 byte (<= less than 256 ompKind options)
 * - interactionType: 1 byte (<= less than interaction types)
 * Todo we probably don't need to send this as both sender and receiver should know interaction type
 * - cell size factor: 8 byte double
 * ToDo We probably don't need this to be so high precision. I wouldn't be surprised if a single byte suffices.
 * - omp chunk size: 2 byte unsigned short int (< 65536 => a reasonable assumption -> there's a good chance we only
 *   need 1 byte.)
 * */
using SerializedConfiguration = std::array<std::byte, 17>;

/**
 * Simply a shorter way of static_casting from Option to std::byte.
 * @tparam TOption
 * @param option
 * @return
 */
template <typename TOption>
inline std::byte castToByte(TOption option) {
  return static_cast<std::byte>(static_cast<typename TOption::Value>(option));
}

/**
 * Distribute ranks in buckets, which contain only ranks with similar scenarios.
 * Each bucket then has its own search space.
 *
 * If negative (=invalid) values for pdBinStdDevDensity or pdBinMaxDensity are passed, every rank gets its own bucket.
 *
 * @param comm MPI communicator
 * @param bucket new MPI communicator for its bucket
 * @param smoothedPDBinStdDevDensity standard deviation in particle-dependent bin density smoothed over last 10
 * iterations. See LiveInfo::gather() for more details.
 * @param smoothedPDBinMaxDensity maximum particle-dependent bin density smoothed over last 10 iterations. See
 * LiveInfo::gather() for more details.
 * @param MPITuningMaxDifferenceForBucket For MPI-tuning: Maximum of the relative difference in the comparison metric
 * for two ranks which exchange their tuning information.
 * @param MPITuningWeightForMaxDensity For MPI-tuning: Weight for maxDensity in the calculation for bucket distribution.
 */
void distributeRanksInBuckets(AutoPas_MPI_Comm comm, AutoPas_MPI_Comm *bucket, double smoothedPDBinStdDevDensity,
                              double smoothedPDBinMaxDensity, double MPITuningMaxDifferenceForBucket,
                              double MPITuningWeightForMaxDensity);

/**
 * Serializes a configuration object for communication via MPI.
 * @param configuration: the configuration to be sent.
 * @return The serialization
 */
SerializedConfiguration serializeConfiguration(Configuration configuration);

/**
 * Serialize a vector of configuration objects into a vector of bytes via serializeConfiguration().
 * @param configurations
 * @return Concatenated byte representations of configurations.
 */
std::vector<std::byte> serializeConfigurations(const std::vector<Configuration> &configurations);

/**
 * Recreates a Configuration object from the object obtained by _serializeConfiguration.
 * @param config: The SerializedConfiguration objects returned by _serializeConfiguration.
 * @return The deserialized Configuration object.
 */
Configuration deserializeConfiguration(SerializedConfiguration config);

/**
 * Deserialize a vector of bytes into a vector of configurations.
 * @param configurationsSerialized
 * @return
 */
std::vector<Configuration> deserializeConfigurations(const std::vector<std::byte> &configurationsSerialized);

/**
 * Handles communication to select the globally best configuration.
 * @param comm: The communicator used for sending and receiving the optimal configuration.
 * @param localOptimalConfig: The locally optimal configuration to be compared with others.
 * @param localOptimalTime: The time measured for localOptimalConfig.
 * @return The globally optimal configuration.
 */
Configuration findGloballyBestConfiguration(AutoPas_MPI_Comm comm, Configuration localOptimalConfig,
                                            long localOptimalTime);

/**
 * Gather the configuration vectors of all ranks in the communicator into one vector at the root process.
 * @param comm
 * @param localConfigurations
 * @param root
 * @return Combined vector of configurations
 */
std::vector<Configuration> gatherConfigurations(AutoPas_MPI_Comm comm,
                                                const std::vector<Configuration> &localConfigurations, int root);

}  // namespace autopas::utils::AutoPasConfigurationCommunicator
