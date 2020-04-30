/**
 * @file AutoPasMessagePasser.h
 * @author W. Thieme
 * @date 04/30/20
 */

#include "WrapMPI.h"
#include "autopas/selectors/Configuration.h"

namespace autopas {

Configuration globalOptimalConfiguration(Configuration config, size_t traversalTime) {
  int worldRank;
  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &worldRank);

  struct {
    size_t val;
    int rank;
  }in[1], out[1];
  in->val = traversalTime;
  in->rank = worldRank;
  AutoPas_MPI_Allreduce(in, out, 1, AUTOPAS_MPI_LONG_INT, AUTOPAS_MPI_MINLOC, AUTOPAS_MPI_COMM_WORLD);

  // The rank with the best configuration sends it to all other ranks
  struct {
    int container, traversal, dataLayout, newton3;
    double cellSizeFactor;
  } resultStruct;
  if (out->rank == worldRank) {
    resultStruct.container = static_cast<ContainerOption::Value>(config.container);
    resultStruct.cellSizeFactor = config.cellSizeFactor;
    resultStruct.traversal = static_cast<TraversalOption::Value>(config.traversal);
    resultStruct.dataLayout = static_cast<DataLayoutOption::Value>(config.dataLayout);
    resultStruct.newton3 = static_cast<Newton3Option::Value>(config.newton3);
  }
  AutoPas_MPI_Bcast(&resultStruct, sizeof(resultStruct), AUTOPAS_MPI_BYTE, out->rank, AUTOPAS_MPI_COMM_WORLD);

  Configuration result = Configuration(
          static_cast<ContainerOption::Value>(resultStruct.container),
          config.cellSizeFactor,
          static_cast<TraversalOption::Value>(resultStruct.traversal),
          static_cast<DataLayoutOption::Value>(resultStruct.dataLayout),
          static_cast<Newton3Option::Value>(resultStruct.newton3));

  AutoPasLog(debug, "Global optimal configuration: {}; took {} nanoseconds", result.toString(), out->val);
  return result;
}

}