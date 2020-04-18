/**
 * @file WrapMPI.h
 * @author W. Thieme
 * @date 4/17/20
 */

#pragma once

/**
 * Provide non-MPI versions of the needed MPI function calls.
 *
 * Extend when necessary.
 */

// @todo find a more elegant way of dealing with the MPI enums
enum Autopas_MPI_Comm {
  AUTOPAS_MPI_COMM_WORLD,
};

enum Autopas_MPI_Error {
  AUTOPAS_MPI_SUCCESS,
  AUTOPAS_MPI_ERR_COMM,
  AUTOPAS_MPI_ERR_ARG,
  AUTOPAS_UNKNOWN_MPI_ERR,
};


#if defined(AUTOPAS_MPI)
#include <mpi.h>
#endif

namespace autopas{

#if defined(AUTOPAS_MPI)

/**
 * Translator for the return types of MPI functions
 * @param err: a value return by an MPI function call
 * @return: the appropriate Autopas equivalent
 */
int translate(int err) {
  switch (err) {
    case MPI_SUCCESS:
      return AUTOPAS_MPI_SUCCESS;
    case MPI_ERR_COMM:
      return AUTOPAS_MPI_ERR_COMM;
    case MPI_ERR_ARG:
      return AUTOPAS_MPI_ERR_ARG;
    default:
      return AUTOPAS_UNKNOWN_MPI_ERR;
  }
}

/**
 * Wrapper for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: outputs number of processes in the group of comm
 * @return: Autopas_MPI error value
 */
inline int Autopas_MPI_Comm_size(Autopas_MPI_Comm comm, int *size) { 
  switch(comm) {
    case AUTOPAS_MPI_COMM_WORLD:
      return translate(MPI_Comm_size(MPI_COMM_WORLD, size));
    default:
      return AUTOPAS_MPI_ERR_COMM;
  }
}

/**
 * Wrapper for MPI_Comm_rank
 * @param comm: communicator (handle)
 * @param rank: outputs rank of the process
 * @return: MPI error value
 */
inline int Autopas_MPI_Comm_rank(Autopas_MPI_Comm comm, int *rank) { 
  switch(comm) {
    case AUTOPAS_MPI_COMM_WORLD:
      return translate(MPI_Comm_rank(MPI_COMM_WORLD, rank));
    default:
      return AUTOPAS_MPI_ERR_COMM;
  }
}

#else

/**
 * Dummy for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: always outputs 1
 * @return: always returns AUTOPAS_MPI_SUCCESS
 */
inline int Autopas_MPI_Comm_size(Autopas_MPI_Comm comm, int *size) {
  if (nullptr == size) {
    return AUTOPAS_MPI_ERR_ARG;
  }
  *size = 1;
  return AUTOPAS_MPI_SUCCESS;
}

/**
 * Dummy for MPI_Comm_rank
 * @param comm: communicator (handle)
 * @param rank: always outputs 0
 * @return: always returns AUTOPAS_MPI_SUCCESS
 */
inline int Autopas_MPI_Comm_rank(Autopas_MPI_Comm comm, int *rank) {
  if (nullptr == rank) {
    return AUTOPAS_MPI_ERR_ARG;
  }
  *rank = 0;
  return AUTOPAS_MPI_SUCCESS;
}

#endif
}
