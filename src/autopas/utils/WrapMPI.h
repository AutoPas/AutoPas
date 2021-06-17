/**
 * @file WrapMPI.h
 * @author W. Thieme
 * @date 17.04.2020
 */

#pragma once

/**
 * Provide non-MPI versions of the needed MPI function calls.
 *
 * May be extended when necessary.
 */

#if defined(AUTOPAS_INTERNODE_TUNING)
#include <mpi.h>
#else
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#endif

namespace autopas {

#if defined(AUTOPAS_INTERNODE_TUNING)

// MPI_Comm
/** Wrapper for MPI_COMM_NULL */
#define AUTOPAS_MPI_COMM_NULL MPI_COMM_NULL
/** Wrapper for MPI_COMM_WORLD */
#define AUTOPAS_MPI_COMM_WORLD MPI_COMM_WORLD

// MPI_Datatype
/** Wrapper for MPI_BYTE */
#define AUTOPAS_MPI_BYTE MPI_BYTE
/** Wrapper for MPI_CXX_BOOL */
#define AUTOPAS_MPI_CXX_BOOL MPI_CXX_BOOL
/** Wrapper for MPI_CHAR */
#define AUTOPAS_MPI_CHAR MPI_CHAR
/** Wrapper for MPI_INT */
#define AUTOPAS_MPI_INT MPI_INT
/** Wrapper for MPI_UNSIGNED LONG */
#define AUTOPAS_MPI_UNSIGNED_LONG MPI_UNSIGNED_LONG

// MPI_Op
/** Wrapper for MPI_LAND */
#define AUTOPAS_MPI_LAND MPI_LAND
/** Wrapper for MPI_MIN */
#define AUTOPAS_MPI_MIN MPI_MIN
/** Wrapper for MPI_MINLOC */
#define AUTOPAS_MPI_MINLOC MPI_MINLOC

// MPI_Status
/** Wrapper for MPI_STATUS IGNORE */
#define AUTOPAS_MPI_STATUS_IGNORE MPI_STATUS_IGNORE

// MPI_Request
/** Wrapper for MPI_REQUEST_NULL */
#define AUTOPAS_MPI_REQUEST_NULL MPI_REQUEST_NULL

/** Wrapper for MPI_MAX_ERROR_STRING */
#define AUTOPAS_MPI_MAX_ERROR_STRING MPI_MAX_ERROR_STRING

using AutoPas_MPI_Comm = MPI_Comm;
using AutoPas_MPI_Datatype = MPI_Datatype;
using AutoPas_MPI_Op = MPI_Op;
using AutoPas_MPI_Status = MPI_Status;
using AutoPas_MPI_Request = MPI_Request;

#else

/**
 * Dummy for MPI_Comm.
 */
enum AutoPas_MPI_Comm {
  AUTOPAS_MPI_COMM_NULL = 0,
  AUTOPAS_MPI_COMM_WORLD,
};

/**
 * Dummy for MPI_Datatype.
 * initialize values to the size of the respective type in bytes.
 */
enum AutoPas_MPI_Datatype {
  AUTOPAS_MPI_BYTE = 1,
  AUTOPAS_MPI_CXX_BOOL = sizeof(bool),
  AUTOPAS_MPI_CHAR = sizeof(char),
  AUTOPAS_MPI_INT = sizeof(int),
  AUTOPAS_MPI_UNSIGNED_LONG = sizeof(unsigned long),
};

/**
 * Dummy for MPI_Op.
 */
enum AutoPas_MPI_Op {
  AUTOPAS_MPI_LAND,
  AUTOPAS_MPI_MIN,
  AUTOPAS_MPI_MINLOC,
};

/**
 * @struct AutoPas_MPI_Status
 * Dummy for MPI_Status
 * @var AutoPas_MPI_Status::count
 * additional field that the MPI standard does not necessitate, but that is often used in implementations of MPI_Status
 * @var AutoPas_MPI_Status::cancelled
 * additional field that the MPI standard does not necessitate, but that is often used in implementations of MPI_Status
 * @var AutoPas_MPI_Status::AUTOPAS_MPI_SOURCE
 * Dummy for MPI_Status::MPI_SOURCE
 * @var AutoPas_MPI_Status::AUTOPAS_MPI_TAG
 * Dummy for MPI_Status::MPI_TAG
 * @var AutoPas_MPI_Status::AUTOPAS_MPI_ERROR
 * Dummy for MPI_Status::MPI_ERROR
 */
struct AutoPas_MPI_Status {
  int count, cancelled, AUTOPAS_MPI_SOURCE, AUTOPAS_MPI_TAG, AUTOPAS_MPI_ERROR;
};
/** Dummy for MPI_STATUS_IGNORE */
#define AUTOPAS_MPI_STATUS_IGNORE nullptr

/**
 * Dummy for MPI_Request.
 */
enum AutoPas_MPI_Request {
  AUTOPAS_MPI_REQUEST_NULL,
  _AUTOPAS_MPI_COMPLETED_REQUEST,
  _AUTOPAS_MPI_INCOMPLETE_REQUEST,
};

/**
 * Dummy for MPI_Error
 */
enum AutoPas_MPI_Error {
  AUTOPAS_MPI_SUCCESS = 0,
  AUTOPAS_MPI_ERR_ARG,
  AUTOPAS_MPI_ERR_COMM,
  AUTOPAS_MPI_ERR_RANK,
  AUTOPAS_MPI_ERR_REQUEST,
  AUTOPAS_MPI_ERR_TYPE,
};
/** Dummy for MPI_MAX_ERROR_STRING */
#define AUTOPAS_MPI_MAX_ERROR_STRING 256

#endif

/**
 * Wrapper for MPI_Init
 * Also defines the AutoPas communicator, so it needs to be called when AutoPas should use MPI
 * @param argc: Pointer to the number of arguments
 * @param argv: Pointer to the argument vector
 * @return MPI error value
 */
inline int AutoPas_MPI_Init(int *argc, char ***argv);

/**
 * Wrapper for MPI_Finalize
 * Also frees the AutoPas communicator, so it needs to be called if AutoPas_MPI_Init was called
 * @return MPI error value
 */
inline int AutoPas_MPI_Finalize();

/**
 * Wrapper for MPI_Finalized
 * @param flag: returns true if (AutoPas_)MPI_Finalize has been called
 * @return MPI error value
 */
inline int AutoPas_MPI_Finalized(int *flag);

/**
 * Wrapper for MPI_Error_string
 * @param errorcode: MPI error value
 * @param string: output string
 * @param resultlen: length of output
 * @return MPI error value
 */
inline int AutoPas_MPI_Error_string(int errorcode, char *string, int *resultlen);

/**
 * Wrapper for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: outputs number of processes in the group of comm
 * @return: MPI error value
 */
inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size);

/**
 * Wrapper for MPI_Comm_rank
 * @param comm: communicator (handle)
 * @param rank: outputs rank of the process
 * @return: MPI error value
 */
inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank);

/**
 * Wrapper for MPI_Comm_dup
 * @param comm: Communicator to be duplicated (handle)
 * @param newComm: outputs new communicator over the same group as comm
 * @return MPI error value
 */
inline int AutoPas_MPI_Comm_dup(AutoPas_MPI_Comm comm, AutoPas_MPI_Comm *newComm);

/**
 * Wrapper for MPI_Comm_free
 * @param comm: communicator to be freed (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Comm_free(AutoPas_MPI_Comm *comm);

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
inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                            AutoPas_MPI_Comm comm);

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
                            AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status);

/**
 * Wrapper for MPI_Bcast
 * @param buffer: send buffer for root. receive buffer for others
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm);

/**
 * Wrapper for MPI_Ibcast
 * @param buffer: send buffer for root. receive buffer for others
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return
 */
inline int AutoPas_MPI_Ibcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm,
                              AutoPas_MPI_Request *request);

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
inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                 AutoPas_MPI_Op op, AutoPas_MPI_Comm comm);

/**
 * Wrapper for MPI_Iallreduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs receive buffer
 * @param count: number of elements in send send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return
 */
inline int AutoPas_MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                  AutoPas_MPI_Op op, AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request);

/**
 * Wrapper for MPI_Barrier
 * @param comm communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Barrier(AutoPas_MPI_Comm comm);

/**
 * Wrapper for MPI_Ibarrier
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Ibarrier(AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request);

/**
 * Wrapper for MPI_Test
 * @param request: request to be tested.
 * @param flag: outputs true if operation complete
 * @param status: outputs status object. May be AUTOPAS_MPI_STATUS_IGNORE
 * @return MPI error value
 */
inline int AutoPas_MPI_Test(AutoPas_MPI_Request *request, int *flag, AutoPas_MPI_Status *status);

/**
 * Wrapper for MPI_Wait
 * @param request: request to be waited for.
 * @param status: outputs status object. May be AUTOPAS_MPI_STATUS_IGNORE
 * @return MPI error value
 */
inline int AutoPas_MPI_Wait(AutoPas_MPI_Request *request, AutoPas_MPI_Status *status);

/**
 * Wrapper for MPI_Request_free
 * @param request: request to be freed (handle). Will be set to AUTOPAS_MPI_REQUEST_NULL
 * @return MPI error value
 */
inline int AutoPas_MPI_Request_free(AutoPas_MPI_Request *request);

/**
 * Wrapper for MPI_Cart_create
 * @param comm: The AutoPas_MPI_Communicator for which to generate the cartesian grid.
 * @param nDims: The number of dimensions in the resulting cartesian grid.
 * @param dims: The size of the cartesian grid in each dimension.
 * @param periods: An array defining for each dimension if it has periodic boundaries.
 * @param reorder: Defines if ranking may be reordered (true) or not (false).
 * @param cartCommunicator: The resulting MPI communicator.
 */
inline int AutoPas_MPI_Cart_create(AutoPas_MPI_Comm comm, int nDims, const int *dims, const int *periods, int reorder,
                                   AutoPas_MPI_Comm *comm_cart);

/**
 * Wrapper for MPI_Cart_get.
 * @param comm: Communicator with Cartesian structure (handle).
 * @param maxdims: Length of vectors dims, periods, and coords in the calling program (integer).
 * @param dims: Number of processes for each Cartesian dimension (array of integers).
 * @param periods: Periodicity (true/false) for each Cartesian dimension (array of logicals).
 * @param coords: Coordinates of calling process in Cartesian structure (array of integers).
 */
inline int AutoPas_MPI_Cart_get(AutoPas_MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]);

/**
 * Wrapper for MPI_Isend.
 * @param buf: Initial address of send buffer (choice).
 * @param count: Number of elements in send buffer (integer).
 * @param datatype: Datatype of each send buffer element (handle).
 * @param dest: Rank of destination (integer).
 * @param tag: Message tag (integer).
 * @param comm: Communicator (handle).
 */
inline int AutoPas_MPI_Isend(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                             AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request);

/**
 * Wrapper for MPI_Probe.
 * @param source: Source rank (integer).
 * @param tag: Tag value (integer).
 * @param comm: Communicator (handle).
 */
inline int AutoPas_MPI_Probe(int source, int tag, AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status);

/**
 * Wrapper for MPI_Get_count.
 * @param status: Return status of receive operation (status).
 * @param datatype: Datatype of each receive buffer element (handle).
 * @param count: Number of received elements (integer).
 */
inline int AutoPas_MPI_Get_count(const AutoPas_MPI_Status *status, AutoPas_MPI_Datatype datatype, int *count);

/**
 * Wrapper for MPI_Waitall.
 * @param count: Lists length (integer).
 * @param array_of_requests: Array of requests (array of handles).
 * @param array_of_statuses: Array of status objects (array of status).
 */
inline int AutoPas_MPI_Waitall(int count, AutoPas_MPI_Request array_of_requests[],
                               AutoPas_MPI_Status *array_of_statuses);

#if defined(AUTOPAS_INTERNODE_TUNING)

inline int AutoPas_MPI_Init(int *argc, char ***argv) { return MPI_Init(argc, argv); }

inline int AutoPas_MPI_Finalize() { return MPI_Finalize(); }

inline int AutoPas_MPI_Finalized(int *flag) { return MPI_Finalized(flag); }

inline int AutoPas_MPI_Error_string(int errorcode, char *string, int *resultlen) {
  return MPI_Error_string(errorcode, string, resultlen);
}

inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size) { return MPI_Comm_size(comm, size); }

inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank) { return MPI_Comm_rank(comm, rank); }

inline int AutoPas_MPI_Comm_dup(AutoPas_MPI_Comm comm, AutoPas_MPI_Comm *newComm) {
  return MPI_Comm_dup(comm, newComm);
}

inline int AutoPas_MPI_Comm_free(AutoPas_MPI_Comm *comm) { return MPI_Comm_free(comm); }

inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                            AutoPas_MPI_Comm comm) {
  return MPI_Send(buf, count, datatype, dest, tag, comm);
}

inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_MPI_Datatype datatype, int source, int tag,
                            AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  return MPI_Recv(buf, count, datatype, source, tag, comm, status);
}

inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  return MPI_Bcast(buffer, count, datatype, root, comm);
}

inline int AutoPas_MPI_Ibcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm,
                              AutoPas_MPI_Request *request) {
  return MPI_Ibcast(buffer, count, datatype, root, comm, request);
}

inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                 AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

inline int AutoPas_MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                  AutoPas_MPI_Op op, AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  return MPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
}

inline int AutoPas_MPI_Barrier(AutoPas_MPI_Comm comm) { return MPI_Barrier(comm); }

inline int AutoPas_MPI_Ibarrier(AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  return MPI_Ibarrier(comm, request);
}

inline int AutoPas_MPI_Test(AutoPas_MPI_Request *request, int *flag, AutoPas_MPI_Status *status) {
  return MPI_Test(request, flag, status);
}

inline int AutoPas_MPI_Wait(AutoPas_MPI_Request *request, AutoPas_MPI_Status *status) {
  return MPI_Wait(request, status);
}

inline int AutoPas_MPI_Request_free(AutoPas_MPI_Request *request) { return MPI_Request_free(request); }

inline int AutoPas_MPI_Cart_create(AutoPas_MPI_Comm comm, int nDims, const int *dims, const int *periods, int reorder,
                                   AutoPas_MPI_Comm *comm_cart) {
  return MPI_Cart_create(comm, nDims, dims, periods, reorder, comm_cart);
}

inline int AutoPas_MPI_Cart_get(AutoPas_MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) {
  return MPI_Cart_get(comm, maxdims, dims, periods, coords);
}

inline int AutoPas_MPI_Isend(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                             AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}

inline int AutoPas_MPI_Probe(int source, int tag, AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  return MPI_Probe(source, tag, comm, status);
}

inline int AutoPas_MPI_Get_count(const AutoPas_MPI_Status *status, AutoPas_MPI_Datatype datatype, int *count) {
  return MPI_Get_count(status, datatype, count);
}

inline int AutoPas_MPI_Waitall(int count, AutoPas_MPI_Request array_of_requests[],
                               AutoPas_MPI_Status *array_of_statuses) {
  return MPI_Waitall(count, array_of_requests, array_of_statuses);
}

#else

inline int AutoPas_MPI_Init(int *argc, char ***argv) { return AUTOPAS_MPI_SUCCESS; }

inline int AutoPas_MPI_Finalize() { return AUTOPAS_MPI_SUCCESS; }

inline int AutoPas_MPI_Finalized(int *flag) {
  *flag = 1;
  return AUTOPAS_MPI_SUCCESS;
}

int AutoPas_MPI_Error_string(int errorcode, char *string, int *resultlen) {
  static const std::map<int, const char *> errorStrings = {
      {AUTOPAS_MPI_SUCCESS, "MPI_SUCCESS: no errors"},
      {AUTOPAS_MPI_ERR_ARG, "MPI_ERR_ARG: invalid argument of some other kind"},
      {AUTOPAS_MPI_ERR_COMM, "MPI_ERR_COMM: invalid communicator"},
      {AUTOPAS_MPI_ERR_RANK, "MPI_ERR_RANK: invalid rank"},
      {AUTOPAS_MPI_ERR_REQUEST, "MPI_ERR_REQUEST: invalid Request"},
      {AUTOPAS_MPI_ERR_TYPE, "MPI_ERR_TYPE: invalid datatype"},
  };
  snprintf(string, AUTOPAS_MPI_MAX_ERROR_STRING, "%s", errorStrings.at(errorcode));
  *resultlen = strnlen(string, AUTOPAS_MPI_MAX_ERROR_STRING);
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Comm_size(AutoPas_MPI_Comm comm, int *size) {
  if (nullptr == size) {
    return AUTOPAS_MPI_ERR_ARG;
  }
  *size = 1;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Comm_rank(AutoPas_MPI_Comm comm, int *rank) {
  if (nullptr == rank) {
    return AUTOPAS_MPI_ERR_ARG;
  }
  *rank = 0;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Comm_dup(AutoPas_MPI_Comm comm, AutoPas_MPI_Comm *newComm) {
  *newComm = comm;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Comm_free(AutoPas_MPI_Comm *comm) {
  *comm = AUTOPAS_MPI_COMM_NULL;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Send(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                            AutoPas_MPI_Comm comm) {
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Recv(void *buf, int count, AutoPas_MPI_Datatype datatype, int source, int tag,
                            AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  // if source==0 there is an error, because you cannot receive from yourself
  // if source >0 there is an error, because there only exists one process
  return AUTOPAS_MPI_ERR_RANK;
}

inline int AutoPas_MPI_Bcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm) {
  if (root > 0) {
    return AUTOPAS_MPI_ERR_RANK;
  } else {
    return AUTOPAS_MPI_SUCCESS;
  }
}

inline int AutoPas_MPI_Ibcast(void *buffer, int count, AutoPas_MPI_Datatype datatype, int root, AutoPas_MPI_Comm comm,
                              AutoPas_MPI_Request *request) {
  if (root > 0) {
    *request = AUTOPAS_MPI_REQUEST_NULL;
    return AUTOPAS_MPI_ERR_RANK;
  } else {
    *request = _AUTOPAS_MPI_COMPLETED_REQUEST;
    return AUTOPAS_MPI_SUCCESS;
  }
}

inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                 AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  memcpy(recvbuf, sendbuf, datatype * static_cast<size_t>(count));
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                  AutoPas_MPI_Op op, AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  *request = _AUTOPAS_MPI_COMPLETED_REQUEST;
  return AutoPas_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

inline int AutoPas_MPI_Barrier(AutoPas_MPI_Comm comm) { return AUTOPAS_MPI_SUCCESS; }

inline int AutoPas_MPI_Ibarrier(AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  *request = _AUTOPAS_MPI_COMPLETED_REQUEST;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Test(AutoPas_MPI_Request *request, int *flag, AutoPas_MPI_Status *status) {
  *request = AUTOPAS_MPI_REQUEST_NULL;
  *flag = 1;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Wait(AutoPas_MPI_Request *request, AutoPas_MPI_Status *status) {
  *request = AUTOPAS_MPI_REQUEST_NULL;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Request_free(AutoPas_MPI_Request *request) {
  if (*request != _AUTOPAS_MPI_COMPLETED_REQUEST) {
    return AUTOPAS_MPI_ERR_REQUEST;
  } else {
    *request = AUTOPAS_MPI_REQUEST_NULL;
    return AUTOPAS_MPI_SUCCESS;
  }
}

inline int AutoPas_MPI_Cart_create(AutoPas_MPI_Comm comm, int nDims, const int *dims, const int *periods, int reorder,
                                   AutoPas_MPI_Comm *comm_cart) {
  *comm_cart = AUTOPAS_MPI_COMM_WORLD;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Cart_get(AutoPas_MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) {
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Isend(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                             AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Probe(int source, int tag, AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status) {
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Get_count(const AutoPas_MPI_Status *status, AutoPas_MPI_Datatype datatype, int *count) {
  *count = 0;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Waitall(int count, AutoPas_MPI_Request array_of_requests[],
                               AutoPas_MPI_Status *array_of_statuses) {
  return AUTOPAS_MPI_SUCCESS;
}

#endif
}  // namespace autopas
