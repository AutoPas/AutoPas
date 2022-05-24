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

#include <limits.h>
#include <stdint.h>

#if defined(AUTOPAS_INCLUDE_MPI)
#include <mpi.h>
#else
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#endif

namespace autopas {

#if defined(AUTOPAS_INCLUDE_MPI)

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
/** Wrapper for MPI_UNSIGNED */
#define AUTOPAS_MPI_UNSIGNED_INT MPI_UNSIGNED
/** Wrapper for MPI_LONG */
#define AUTOPAS_MPI_LONG MPI_LONG
/** Wrapper for MPI_UNSIGNED_LONG */
#define AUTOPAS_MPI_UNSIGNED_LONG MPI_UNSIGNED_LONG
/** Wrapper for MPI_DOUBLE */
#define AUTOPAS_MPI_DOUBLE MPI_DOUBLE

// MPI_Op
/** Wrapper for MPI_LAND */
#define AUTOPAS_MPI_LAND MPI_LAND
/** Wrapper for MPI_MIN */
#define AUTOPAS_MPI_MIN MPI_MIN
/** Wrapper for MPI_MINLOC */
#define AUTOPAS_MPI_MINLOC MPI_MINLOC
/** Wrapper for MPI_SUM */
#define AUTOPAS_MPI_SUM MPI_SUM

// MPI Constants
/** Wrapper for MPI_IN_PLACE  */
#define AUTOPAS_MPI_IN_PLACE MPI_IN_PLACE

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
  COMM_NULL = 0,
  COMM_WORLD,
};
/** Wrapper for MPI_COMM_NULL */
#define AUTOPAS_MPI_COMM_NULL autopas::AutoPas_MPI_Comm::COMM_NULL
/** Wrapper for MPI_COMM_WORLD */
#define AUTOPAS_MPI_COMM_WORLD autopas::AutoPas_MPI_Comm::COMM_WORLD

/**
 * Dummy for MPI_Datatype.
 * initialize values to the size of the respective type in bytes.
 */
enum AutoPas_MPI_Datatype {
  BYTE = 1,
  CXX_BOOL = sizeof(bool),
  CHAR = sizeof(char),
  UNSIGNED_CHAR = sizeof(unsigned char),
  UNSIGNED_SHORT = sizeof(unsigned short),
  INT = sizeof(int),
  UNSIGNED_INT = sizeof(unsigned int),
  UNSIGNED_LONG = sizeof(unsigned long),
  UNSIGNED_LONG_LONG = sizeof(unsigned long long),
  LONG = sizeof(double),
  DOUBLE = sizeof(double)
};
// MPI_Datatype
/** Wrapper for MPI_BYTE */
#define AUTOPAS_MPI_BYTE autopas::AutoPas_MPI_Datatype::BYTE
/** Wrapper for MPI_CXX_BOOL */
#define AUTOPAS_MPI_CXX_BOOL autopas::AutoPas_MPI_Datatype::CXX_BOOL
/** Wrapper for MPI_CHAR */
#define AUTOPAS_MPI_CHAR autopas::AutoPas_MPI_Datatype::CHAR
/** Wrapper for MPI_INT */
#define AUTOPAS_MPI_INT autopas::AutoPas_MPI_Datatype::INT
/** Wrapper for MPI_UNSIGNED */
#define AUTOPAS_MPI_UNSIGNED_INT autopas::AutoPas_MPI_Datatype::UNSIGNED_INT
/** Wrapper for MPI_LONG */
#define AUTOPAS_MPI_LONG autopas::AutoPas_MPI_Datatype::LONG
/** Wrapper for MPI_UNSIGNED LONG */
#define AUTOPAS_MPI_UNSIGNED_LONG autopas::AutoPas_MPI_Datatype::UNSIGNED_LONG
/** Wrapper for MPI_DOUBLE */
#define AUTOPAS_MPI_DOUBLE autopas::AutoPas_MPI_Datatype::DOUBLE

/**
 * Dummy for MPI_Op.
 */
enum AutoPas_MPI_Op { LAND, MIN, MINLOC, SUM };
// MPI_Op
/** Wrapper for MPI_LAND */
#define AUTOPAS_MPI_LAND autopas::AutoPas_MPI_Op::LAND
/** Wrapper for MPI_MIN */
#define AUTOPAS_MPI_MIN autopas::AutoPas_MPI_Op::MIN
/** Wrapper for MPI_MINLOC */
#define AUTOPAS_MPI_MINLOC autopas::AutoPas_MPI_Op::MINLOC
/** Wrapper for MPI_SUM */
#define AUTOPAS_MPI_SUM autopas::AutoPas_MPI_Op::SUM

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
  REQUEST_NULL,
  COMPLETED_REQUEST,
  INCOMPLETE_REQUEST,
};
/** Dummy for MPI_REQUEST_NULL */
#define AUTOPAS_MPI_REQUEST_NULL autopas::AutoPas_MPI_Request::REQUEST_NULL

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
/** Indicator for Collectives to happe in-place */
#define AUTOPAS_MPI_IN_PLACE ((void *)1)

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
 * Wrapper for MPI_Reduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs receive buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param root: the rank of the root process
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int AutoPas_MPI_Reduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                              AutoPas_MPI_Op op, int root, AutoPas_MPI_Comm comm);

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
 * @param comm_cart: The resulting MPI communicator.
 * @return MPI error value
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
 * @return MPI error value
 */
inline int AutoPas_MPI_Cart_get(AutoPas_MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]);

/**
 * Wrapper for MPI_Cart_coords
 * Determines process coords in cartesian topology given rank in group.
 * @param comm Communicator with cartesian structure (handle).
 * @param rank Rank of a process within group of comm (integer).
 * @param maxdims Length of vector coords in the calling program (integer).
 * @param coords Integer array (of size ndims) containing the Cartesian coordinates of specified process (integer).
 * @return MPI error value
 */
inline int AutoPas_MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
/**
 * Wrapper for MPI_Dims_create.
 * Creates a division of processors in a cartesian grid. Typically this is a factorization of nnodes in ndims factors.
 * @param nnodes Number of ranks to divide.
 * @param ndims Number of dimension over which to spread the ranks.
 * @param dims Output parameter. Should be an array of size ndims.
 * @return MPI error value
 */
inline int AutoPas_MPI_Dims_create(int nnodes, int ndims, int dims[]);

/**
 * Wrapper for MPI_Isend.
 * @param buf: Initial address of send buffer (choice).
 * @param count: Number of elements in send buffer (integer).
 * @param datatype: Datatype of each send buffer element (handle).
 * @param dest: Rank of destination (integer).
 * @param tag: Message tag (integer).
 * @param comm: Communicator (handle).
 * @param request: A pointer to the created send request.
 * @return MPI error value
 */
inline int AutoPas_MPI_Isend(const void *buf, int count, AutoPas_MPI_Datatype datatype, int dest, int tag,
                             AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request);

/**
 * Wrapper for MPI_Probe.
 * @param source: Source rank (integer).
 * @param tag: Tag value (integer).
 * @param comm: Communicator (handle).
 * @param status: The status of the probed request.
 * @return MPI error value
 */
inline int AutoPas_MPI_Probe(int source, int tag, AutoPas_MPI_Comm comm, AutoPas_MPI_Status *status);

/**
 * Wrapper for MPI_Get_count.
 * @param status: Return status of receive operation (status).
 * @param datatype: Datatype of each receive buffer element (handle).
 * @param count: Number of received elements (integer).
 * @return MPI error value
 */
inline int AutoPas_MPI_Get_count(const AutoPas_MPI_Status *status, AutoPas_MPI_Datatype datatype, int *count);

/**
 * Wrapper for MPI_Waitall.
 * @param count: Lists length (integer).
 * @param array_of_requests: Array of requests (array of handles).
 * @param array_of_statuses: Array of status objects (array of status).
 * @return MPI error value
 */
inline int AutoPas_MPI_Waitall(int count, AutoPas_MPI_Request array_of_requests[],
                               AutoPas_MPI_Status *array_of_statuses);

/**
 * Wrapper for MPI_Allgather
 * @param buffer_send: send buffer
 * @param count_send: number of elements in send buffer
 * @param datatype_send: type of elements in send buffer
 * @param buffer_recv: receive buffer
 * @param count_recv: number of elements received from each rank
 * @param datatype_recv: type of elements in receive buffer
 * @param comm: communicator (handle)
 * @return
 */
inline int AutoPas_MPI_Allgather(void *buffer_send, int count_send, AutoPas_MPI_Datatype datatype_send,
                                 void *buffer_recv, int count_recv, AutoPas_MPI_Datatype datatype_recv,
                                 AutoPas_MPI_Comm comm);

/**
 * Wrapper for MPI_Comm_split
 * @param old_communicator: old communicator (handle)
 * @param color: determines which ranks ar in the same bucket
 * @param key: determines rank order in new communicator
 * @param new_communicator: pointer to new communicator
 * @return
 */
inline int AutoPas_MPI_Comm_split(AutoPas_MPI_Comm old_communicator, int color, int key,
                                  AutoPas_MPI_Comm *new_communicator);

#if defined(AUTOPAS_INCLUDE_MPI)

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

inline int AutoPas_MPI_Comm_split(AutoPas_MPI_Comm comm, int color, int key, AutoPas_MPI_Comm *newcomm) {
  return MPI_Comm_split(comm, color, key, newcomm);
}

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

inline int AutoPas_MPI_Reduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                              AutoPas_MPI_Op op, int root, AutoPas_MPI_Comm comm) {
  return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
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

inline int AutoPas_MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) {
  return MPI_Cart_coords(comm, rank, maxdims, coords);
}

inline int AutoPas_MPI_Dims_create(int nnodes, int ndims, int dims[]) { return MPI_Dims_create(nnodes, ndims, dims); }

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

inline int AutoPas_MPI_Allgather(void *buffer_send, int count_send, AutoPas_MPI_Datatype datatype_send,
                                 void *buffer_recv, int count_recv, AutoPas_MPI_Datatype datatype_recv,
                                 AutoPas_MPI_Comm comm) {
  return MPI_Allgather(buffer_send, count_send, datatype_send, buffer_recv, count_recv, datatype_recv, comm);
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
    *request = REQUEST_NULL;
    return AUTOPAS_MPI_ERR_RANK;
  } else {
    *request = COMPLETED_REQUEST;
    return AUTOPAS_MPI_SUCCESS;
  }
}

inline int AutoPas_MPI_Reduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                              AutoPas_MPI_Op op, int root, AutoPas_MPI_Comm comm) {
  if (sendbuf != AUTOPAS_MPI_IN_PLACE) {
    memcpy(recvbuf, sendbuf, datatype * static_cast<size_t>(count));
  }
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                 AutoPas_MPI_Op op, AutoPas_MPI_Comm comm) {
  return AutoPas_MPI_Reduce(sendbuf, recvbuf, count, datatype, op, 0, comm);
}

inline int AutoPas_MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, AutoPas_MPI_Datatype datatype,
                                  AutoPas_MPI_Op op, AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  *request = COMPLETED_REQUEST;
  return AutoPas_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

inline int AutoPas_MPI_Barrier(AutoPas_MPI_Comm comm) { return AUTOPAS_MPI_SUCCESS; }

inline int AutoPas_MPI_Ibarrier(AutoPas_MPI_Comm comm, AutoPas_MPI_Request *request) {
  *request = COMPLETED_REQUEST;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Test(AutoPas_MPI_Request *request, int *flag, AutoPas_MPI_Status *status) {
  *request = REQUEST_NULL;
  *flag = 1;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Wait(AutoPas_MPI_Request *request, AutoPas_MPI_Status *status) {
  *request = REQUEST_NULL;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Request_free(AutoPas_MPI_Request *request) {
  if (*request != COMPLETED_REQUEST) {
    return AUTOPAS_MPI_ERR_REQUEST;
  } else {
    *request = REQUEST_NULL;
    return AUTOPAS_MPI_SUCCESS;
  }
}

inline int AutoPas_MPI_Allgather(void *buffer_send, int count_send, AutoPas_MPI_Datatype datatype_send,
                                 void *buffer_recv, int count_recv, AutoPas_MPI_Datatype datatype_recv,
                                 AutoPas_MPI_Comm comm) {
  if (buffer_send != AUTOPAS_MPI_IN_PLACE) {
    for (long i = 0; i < (count_recv / count_send); i++)
      // offsets from pointers are of type ptrdiff_t which is an alias for long. Hence, i should be long.
      memcpy(static_cast<char *>(buffer_recv) + (i * count_send * sizeof(datatype_send)), buffer_send,
             count_send * sizeof(datatype_send));
  }
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Comm_split(AutoPas_MPI_Comm old_communicator, int color, int key,
                                  AutoPas_MPI_Comm *new_communicator) {
  new_communicator = &old_communicator;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Cart_create(AutoPas_MPI_Comm comm, int nDims, const int *dims, const int *periods, int reorder,
                                   AutoPas_MPI_Comm *comm_cart) {
  *comm_cart = AUTOPAS_MPI_COMM_WORLD;
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Cart_get(AutoPas_MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) {
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) {
  for (int i = 0; i < maxdims; ++i) {
    coords[i] = 0;
  }
  return AUTOPAS_MPI_SUCCESS;
}

inline int AutoPas_MPI_Dims_create(int nnodes, int ndims, int dims[]) {
  // in non-MPI case nnodes should always be 1
  dims[0] = 1;
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
