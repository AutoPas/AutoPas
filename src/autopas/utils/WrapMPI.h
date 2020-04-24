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
enum AutoPas_MPI_Comm {
  AUTOPAS_MPI_COMM_WORLD,
};

enum AutopPas_MPI_Error {
  AUTOPAS_MPI_SUCCESS,
  AUTOPAS_MPI_ERR_COMM,
  AUTOPAS_MPI_ERR_ARG,
  AUTOPAS_MPI_ERR_TYPE,
  AUTOPAS_UNKNOWN_MPI_ERR, // NOT the AutoPas version of MPI_ERR_UNKNOWN
};

enum AutoPas_MPI_Datatype {
  AUTOPAS_CONFIG,
};

enum AutoPas_MPI_Status {
  AUTOPAS_MPI_STATUS_IGNORE,
};

#if defined(AUTOPAS_MPI)
#include <mpi.h>
#include <stddef.h>
#endif

namespace autopas{

/**
 * Extends MPI with a Datatype for AutoPas Configarations
 */
#if defined(AUTOPAS_MPI)
struct Config{
  short container, traversal, dataLayout, newton3;
  double cellSizeFactor;
};

MPI_Datatype _init_config_type() {
  MPI_Datatype config;
  const int array_of_blocklengths[] = {4, 1};
  const MPI_Aint array_of_displacement[] = { offsetof(Config, container), offsetof(Config, cellSizeFactor) };
  const MPI_Datatype array_of_datatypes[] = { MPI_SHORT, MPI_DOUBLE };
  MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacement, array_of_datatypes, &config);
  return config;
}

struct {
  MPI_Datatype AUTOPAS_MPI_CONFIG = _init_config_type();
}_AutoPas_MPI_Datatype_Handles;

/**
 * Translator for the return types of MPI functions
 * @param err: a value return by an MPI function call
 * @return: the appropriate AutoPas equivalent
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
 * @return: AutoPas_MPI error value
 */
inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size) {
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
inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank) {
  switch(comm) {
    case AUTOPAS_MPI_COMM_WORLD:
      return translate(MPI_Comm_rank(MPI_COMM_WORLD, rank));
    default:
      return AUTOPAS_MPI_ERR_COMM;
  }
}

/**
 * Helper for AutoPas_MPI_Send
 * Not to be used externally
 * @param all: same as AutoPas_MPI_Send
 */
inline int _AutoPas_MPI_Send_Helper(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, AutoPas_MPI_Comm comm) {
  switch(comm) {
    case AUTOPAS_MPI_COMM_WORLD:
      return translate(MPI_Send(buf, count, datatype, dest, tag, MPI_COMM_WORLD));
    default:
      return AUTOPAS_MPI_ERR_COMM;
  }
}

/**
 * Wrapper for MPI_Send
 * @param buf: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param dest: rank of destination process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @return AutoPas_MPI error value
 */
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag, AutoPas_MPI_Comm comm) {
  switch(datatype) {
    case AUTOPAS_CONFIG:
      return _AutoPas_MPI_Send_Helper(buf, count, _AutoPas_MPI_Datatype_Handles.AUTOPAS_MPI_CONFIG, dest, tag, comm);
    default:
      return AUTOPAS_MPI_ERR_TYPE;
  }
}


/**
 * Helper for AutoPas_MPI_Recv
 * Not to be used externally
 * @param all: same as AutoPas_MPI_Recv
 */
inline int _AutoPas_MPI_Recv_helper(void *buf, int count, MPI_Datatype datatype, int source, int tag,
        AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  switch(comm) {
    case AUTOPAS_MPI_COMM_WORLD:
      return translate(MPI_Recv(buf, count, datatype, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE/* should be changed*/));
    default:
      return AUTOPAS_MPI_ERR_COMM;
  }
}
/**
 * Wrapper for MPI_Recv
 * @param buf: outputs receive buffer
 * @param count: maximum number of elements in receive buffer
 * @param datatype: type of elements in receive buffer
 * @param source: rank of source process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @param status: outputs status object
 * @return AutoPas_MPI error value
 */
inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_MPI_Datatype datatype, int source, int tag,
             AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  switch(datatype) {
    case AUTOPAS_CONFIG:
      return _AutoPas_MPI_Recv_helper(buf, count, _AutoPas_MPI_Datatype_Handles.AUTOPAS_MPI_CONFIG, source, tag,
              comm, status);
    default:
      return AUTOPAS_MPI_ERR_TYPE;
  }
}

#else

/**
 * Dummy for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: always outputs 1
 * @return: always returns AUTOPAS_MPI_SUCCESS
 */
inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size) {
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
inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank) {
  if (nullptr == rank) {
    return AUTOPAS_MPI_ERR_ARG;
  }
  *rank = 0;
  return AUTOPAS_MPI_SUCCESS;
}

/**
 * Dummy for MPI_Send
 * @param buf: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param dest: rank of destination process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @return always returns AUTOPAS_MPI_SUCCESS
 */
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag, AutoPas_MPI_Comm comm) {
  return AUTOPAS_MPI_SUCCESS;
}

/**
 * Dummy for MPI_Recv
 * @param buf: outputs nullptr
 * @param count: maximum number of elements in receive buffer
 * @param datatype: type of elements in receive buffer
 * @param source: rank of source process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @param status: outputs the input
 * @return always return AUTOPAS_MPI_SUCCESS
 */
inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_MPI_Datatype datatype, int source, int tag,
                            AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  buf = nullptr;
  return AUTOPAS_MPI_SUCCESS;
}

#endif
}
