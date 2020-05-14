/**
 * @file WrapMPI.h
 * @author W. Thieme
 * @date 4/17/20
 */

#pragma once

/**
 * Provide non-MPI versions of the needed MPI function calls.
 *
 * May be extended when necessary.
 */

#if defined(AUTOPAS_MPI)
#include <mpi.h>
namespace autopas {}
#else
#include <set>

namespace autopas {

enum MPI_Comm {
  MPI_COMM_NULL = 0,
  MPI_COMM_WORLD,
};

// initialize values to the size of the respective type in bytes
enum MPI_Datatype {
  MPI_BYTE = 1,
  MPI_CXX_BOOL = sizeof(bool),
  MPI_INT = sizeof(int),
  MPI_LONG_INT = sizeof(long) + sizeof(int),
  MPI_UNSIGNED_LONG = sizeof(unsigned long),
};

enum MPI_Op {
  MPI_LAND,
  MPI_MIN,
  MPI_MINLOC,
};

struct MPI_Status{
  int count, cancelled, MPI_SOURCE, AUTOPAS_MPI_TAG, AUTOPAS_MPI_ERROR;
} MPI_STATUS_IGNORE;

enum MPI_Request {
  MPI_REQUEST_NULL,
  _MPI_COMPLETED_REQUEST,
  _MPI_INCOMPLETE_REQUEST,
};

enum MPI_Error {
  MPI_SUCCESS = 0,
  MPI_ERR_ARG,
  MPI_ERR_COMM,
  MPI_ERR_RANK,
  MPI_ERR_REQUEST,
  MPI_ERR_TYPE,
};
#define MPI_MAX_ERROR_STRING 256

/**
 * Dummy for MPI_Init
 * Also defines the AutoPas communicator, so it needs to be called when AutoPas should use MPI
 * @param argc: Pointer to the number of arguments
 * @param argv: Pointer to the argument vector
 * @return MPI error value
 */
inline int MPI_Init(int *argc, char ***argv) {
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Finalize
 * Also frees the AutoPas communicator, so it needs to be called if MPI_Init was called
 * @return MPI error value
 */
inline int MPI_Finalize() {
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Finalized
 * @param flag: returns true if (AutoPas_)MPI_Finalize has been called
 * @return MPI error value
 */
inline int MPI_Finalized(int *flag) {
  *flag = 1;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Error_string
 * @param errorcode: MPI error value
 * @param string: output string
 * @param resultlen: length of output
 * @return MPI error value
 */
inline int MPI_Error_string(int errorcode, char *string, int *resultlen) {
  static const std::map<int, const char *> errorStrings = {
          {MPI_SUCCESS,     "MPI_SUCCESS: no errors"},
          {MPI_ERR_ARG,     "MPI_ERR_ARG: invalid argument of some other kind"},
          {MPI_ERR_COMM,    "MPI_ERR_COMM: invalid communicator"},
          {MPI_ERR_RANK,    "MPI_ERR_RANK: invalid rank"},
          {MPI_ERR_REQUEST, "MPI_ERR_REQUEST: invalid Request"},
          {MPI_ERR_TYPE,    "MPI_ERR_TYPE: invalid datatype"},
  };
  strcpy(string, errorStrings.at(errorcode));
  *resultlen = strlen(string);
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Comm_size
 * @param comm: communicator (handle)
 * @param size: outputs number of processes in the group of comm
 * @return: MPI error value
 */
inline int MPI_Comm_size(MPI_Comm comm, int *size) {
  if (nullptr == size) {
    return MPI_ERR_ARG;
  }
  *size = 1;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Comm_rank
 * @param comm: communicator (handle)
 * @param rank: outputs rank of the process
 * @return: MPI error value
 */
inline int MPI_Comm_rank(MPI_Comm comm, int *rank) {
  if (nullptr == rank) {
    return MPI_ERR_ARG;
  }
  *rank = 0;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Comm_dup
 * @param comm: Communicator to be duplicated (handle)
 * @param newComm: outputs new communicator over the same group as comm
 * @return MPI error value
 */
inline int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newComm) {
  *newComm = comm;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Comm_free
 * @param comm: communicator to be freed (handle)
 * @return MPI error value
 */
inline int MPI_Comm_free(MPI_Comm *comm) {
  *comm = MPI_COMM_NULL;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Send
 * @param buf: send buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param dest: rank of destination process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Recv
 * @param buf: outputs receive buffer
 * @param count: maximum number of elements in receive buffer
 * @param datatype: type of elements in receive buffer
 * @param source: rank of source process
 * @param tag: message tag
 * @param comm: communicator (handle)
 * @param status: currently ignored
 * @return MPI error value
 */
inline int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                    MPI_Comm comm, MPI_Status *status) {
  // if source==0 there is an error, because you cannot receive from yourself
  // if source >0 there is an error, because there only exists one process
  return MPI_ERR_RANK;
}

/**
 * Dummy for MPI_Bcast
 * @param buffer: send buffer for root. receive buffer for others
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  if (root > 0) {
    return MPI_ERR_RANK;
  } else {
    return MPI_SUCCESS;
  }
}

/**
 * Dummy for MPI_Ibcast
 * @param buffer: send buffer for root. receive buffer for others
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param root: rank of the process sending the broadcast
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return
 */
inline int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                      MPI_Request *request) {
  if (root > 0) {
    *request = MPI_REQUEST_NULL;
    return MPI_ERR_RANK;
  } else {
    *request = _MPI_COMPLETED_REQUEST;
    return MPI_SUCCESS;
  }
}

/**
 * Dummy for MPI_Allreduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs receive buffer
 * @param count: number of elements in send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param comm: communicator (handle)
 * @return MPI error value
 */
inline int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  memcpy(recvbuf, sendbuf, datatype * count);
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Iallreduce
 * @param sendbuf: send buffer
 * @param recvbuf: outputs receive buffer
 * @param count: number of elements in send send buffer
 * @param datatype: type of elements in send buffer
 * @param op: reduction operation (handle)
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return
 */
inline int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
                          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                          MPI_Request *request) {
  *request = _MPI_COMPLETED_REQUEST;
  return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
}

/**
 * Dummy for MPI_Ibarrier
 * @param comm: communicator (handle)
 * @param request: outputs communication request (handle)
 * @return MPI error value
 */
inline int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request) {
  *request = _MPI_COMPLETED_REQUEST;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Test
 * @param request: request to be tested.
 * @param flag: outputs true if operation complete
 * @param status: outputs status object. May be MPI_STATUS_IGNORE
 * @return MPI error value
 */
inline int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) {
  *request = MPI_REQUEST_NULL;
  *flag = 1;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Wait
 * @param request: request to be waited for.
 * @param status: outputs status object. May be MPI_STATUS_IGNORE
 * @return MPI error value
 */
inline int MPI_Wait(MPI_Request *request, MPI_Status *status) {
  *request = MPI_REQUEST_NULL;
  return MPI_SUCCESS;
}

/**
 * Dummy for MPI_Request_free
 * @param request: request to be freed (handle). Will be set to MPI_REQUEST_NULL
 * @return MPI error value
 */
inline int MPI_Request_free(MPI_Request *request) {
  if (*request != _MPI_COMPLETED_REQUEST) {
    return MPI_ERR_REQUEST;
  } else {
    *request = MPI_REQUEST_NULL;
    return MPI_SUCCESS;
  }
}

}
#endif
