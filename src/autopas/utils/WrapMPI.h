/**
 * @file WrapMPI.h
 * @author W. Thieme
 * @date 4/17/20
 */

#pragma once

/**
 * Provide non-MPI versions of the needed MPI function calls.
 * Extend MPI functionality by AutoPas-specifics (e.g. the Config datatype)
 *
 * May be extended when necessary.
 */

// initialize to size  of the type in bytes
enum AutoPas_Datatype {
  AUTOPAS_CONFIG = 16,
};


struct Config_struct {
  char container, traversal, dataLayout, newton3;
  double cellSizeFactor;
};


#if defined(AUTOPAS_MPI)
#include <mpi.h>
#endif

namespace autopas {

/**
 * Extends MPI with a Datatype for AutoPas Configarations
 */
#if defined(AUTOPAS_MPI)

// MPI_Comm
#define AUTOPAS_MPI_COMM_WORLD MPI_COMM_WORLD

// MPI_Datatype
#define AUTOPAS_MPI_CXX_BOOL MPI_CXX_BOOL
#define AUTOPAS_MPI_LONG_INT MPI_LONG_INT

// MPI_Op
#define AUTOPAS_MPI_LAND MPI_LAND
#define AUTOPAS_MPI_MINLOC MPI_MINLOC

// MPI_Status
#define AUTOPAS_MPI_STATUS_IGNORE MPI_STATUS_IGNORE

using AutoPas_MPI_Comm = MPI_Comm;
using AutoPas_MPI_Datatype = MPI_Datatype;
using AutoPas_MPI_Op = MPI_Op;
using AutoPas_MPI_Status = MPI_Status;
using AutoPas_MPI_Request = MPI_Request;

/**
 * Function used internally to define AUTOPAS_MPI_CONFIG as an MPI_Datatype
 * @return Handle for AUTOPAS_MPI_CONFIG
 */
MPI_Datatype _init_config_type() {
  MPI_Datatype config;
  const int array_of_blocklengths[] = {4, 1};
  const MPI_Aint array_of_displacement[] = { offsetof(Config_struct, container), offsetof(Config_struct, cellSizeFactor) };
  const MPI_Datatype array_of_datatypes[] = { MPI_CHAR, MPI_DOUBLE };
  MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacement, array_of_datatypes, &config);
  return config;
}

struct {
  MPI_Datatype AUTOPAS_MPI_CONFIG = nullptr;
}_AutoPas_Datatype_Handles;

/**
 * Gives the MPI handle for a given AutoPas_Datatype
 * @param datatype: the AutoPas_Datatype to convert
 * @return the MPI handle for datatype or MPI_DATATYPE_NULL for undefined inputs
 */
inline MPI_Datatype AutoPas_to_MPI_datatype(AutoPas_Datatype datatype) {
  MPI_Datatype result;
  switch (datatype) {
    case AUTOPAS_CONFIG:
      if (_AutoPas_Datatype_Handles.AUTOPAS_MPI_CONFIG == nullptr) {
        _AutoPas_Datatype_Handles.AUTOPAS_MPI_CONFIG = _init_config_type();
        MPI_Type_commit(&_AutoPas_Datatype_Handles.AUTOPAS_MPI_CONFIG);
      }
      result = _AutoPas_Datatype_Handles.AUTOPAS_MPI_CONFIG; break;
    default:
      return MPI_DATATYPE_NULL;
  }
  return result;
}

/**
 * Wrapper for MPI_Error_string
 * @param errorcode
 * @param string: output string
 * @param resultlen: length of output
 * @return MPI error value
 */
inline int AutoPas_MPI_Error_string(int errorcode, char *string, int *resultlen) {
  return MPI_Error_string(errorcode, string, resultlen);
}

/**
 * Wrapper for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: outputs number of processes in the group of comm
 * @return: MPI error value
 */
inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size) { return MPI_Comm_size(comm, size); }

/**
 * Wrapper for MPI_Comm_rank
 * @param comm: communicator (handle)
 * @param rank: outputs rank of the process
 * @return: MPI error value
 */
inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank) { return MPI_Comm_rank(comm, rank); }

/**
 * Wrapper for MPI_Send
 * @param buf: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param dest: rank of destination process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype,
        int dest, int tag, AutoPas_MPI_Comm comm) {
  return MPI_Send(buf, count, datatype, dest, tag, comm);
}
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_Datatype datatype,
        int dest, int tag, AutoPas_MPI_Comm comm) {
  return MPI_Send(buf, count, AutoPas_to_MPI_datatype(datatype), dest, tag, comm);
}

/**
 * Wrapper for MPI_Recv
 * @param buf: outputs receive buffer
 * @param count: maximum number of elements in receive buffer
 * @param datatype: type of elements in receive buffer
 * @param source: rank of source process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @param status: currently ignored
 * @return MPI error value
 */
inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_MPI_Datatype datatype, int source, int tag,
        AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  return MPI_Recv(buf, count, datatype, source, tag, comm, status);
}
inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_Datatype datatype, int source, int tag,
        AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  return MPI_Recv(buf, count, AutoPas_to_MPI_datatype(datatype), source, tag, comm, status);
}

/**
 * Wrapper for MPI_Bcast
 * @param buffer: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  return MPI_Bcast(buffer, count, datatype, root, comm);
}
inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  return MPI_Bcast(buffer, count, AutoPas_to_MPI_datatype(datatype), root, comm);
}

/**
 * Wrapper for MPI_Allreduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs receive buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
        AutoPas_MPI_Datatype datatype, AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}
inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
        AutoPas_Datatype datatype, AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count, AutoPas_to_MPI_datatype(datatype), op, comm);
}

#else

enum AutoPas_MPI_Comm {
  AUTOPAS_MPI_COMM_WORLD,
};

enum AutoPas_MPI_Error {
  AUTOPAS_MPI_SUCCESS = 0,
  AUTOPAS_MPI_ERR_ARG,
  AUTOPAS_MPI_ERR_COMM,
  AUTOPAS_MPI_ERR_TYPE,
};

// initialize values to the size of the respective type in bytes
enum AutoPas_MPI_Datatype {
  AUTOPAS_MPI_LONG_INT = 12,
};

enum AutoPas_MPI_Op {
  AUTOPAS_MPI_MINLOC,
};

struct _AutoPas_MPI_Status{
  int count, cancelled, AUTOPAS_MPI_SOURCE, AUTOPAS_MPI_TAG, AUTOPAS_MPI_ERROR;
};
using AutoPas_MPI_Status = struct _AutoPas_MPI_Status;

#define AUTOPAS_MPI_STATUS_IGNORE nullptr

/**
 * Dummy for MPI_Error_string
 * @param errorcode
 * @param string: output string
 * @param resultlen: length of output
 * @return AUTOPAS_MPI_SUCCESS or AUTOPAS_MPI_ERR_ARG in case of undefined error code
 */
inline int AutoPas_MPI_Error_string(int errorcode, char *string, int *resultlen) {
  switch (errorcode) {
    // @todo implement error strings
    default:
      return AUTOPAS_MPI_ERR_ARG;
  }
}

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
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_Datatype datatype, int dest, int tag, AutoPas_MPI_Comm comm) {
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
inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_Datatype datatype, int source, int tag,
        AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  buf = nullptr;
  return AUTOPAS_MPI_SUCCESS;
}

/**
 * Dummy for MPI_Bcast
 * @param buffer: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @return AUTOPAS_MPI_SUCCESS
 */
inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  return AUTOPAS_MPI_SUCCESS;
}
inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  return AUTOPAS_MPI_SUCCESS;
}

/**
 * Dummy for MPI_Allreduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs sendbuf
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param comm: communicator (handle)
 * @return AUTOPAS_MPI_SUCCESS
 */
inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                                 AutoPas_MPI_Datatype datatype, AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  memcpy(recvbuf, sendbuf, datatype);
  return AUTOPAS_MPI_SUCCESS;
}
inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                                 AutoPas_Datatype datatype, AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  memcpy(recvbuf, sendbuf, datatype);
  return AUTOPAS_MPI_SUCCESS;
}

#endif
}
