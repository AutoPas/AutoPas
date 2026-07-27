#ifndef ALL_GENERAL_FUNCTION_HEADER_INCLUDED
#define ALL_GENERAL_FUNCTION_HEADER_INCLUDED

#include <mpi.h>
#include <cmath>
#include <vector>

namespace ALL{
namespace Functions {

/// function to compute the one-dimensional shift of the border between the
/// local process and the process indicated by neighbor_rank
/// @tparam T data type for vertices and related data
/// @tparam W data type for work and related data
/// @param[in] remote_rank the MPI rank of the neighbors the border is
/// shared with
/// @param[in] local_coord the cartesian coordinate of the local domain in the
/// process grid in the direction of the border shift
/// @param[in] global_dim the dimension of the process grid in the direction of
/// the border shift
/// @param[in] local_work the work on the local process
/// @param[in] remote_work the work on the neighboring process
/// @param[in] local_size the size of the local domain in the dimension of the
/// border shift
/// @param[in] remote_size the size of the neighbor domain in the dimension of
/// the border shift
/// @param[in] gamma the correction value for the shift, needed to regulate the
/// width of the shift, e.g. to avoid borders shifts where next-neighbor borders
/// are crossed by one another
/// @param[in] minSize optional minimum size of the domain in the dimension of
/// the border shift
/// @result the shift of the border in relation to the local domain
template <typename T, typename W>
T borderShift1d(const int remote_rank, const int local_coord,
                const int global_dim, const W local_work, const W remote_work,
                const T local_size, const T remote_size, const T gamma,
                const T minSize) {
  // calculate shift of borders:
  // s = 0.5 * gamma * (W_r - W_l) / (W_r + W_l) * (d_r + d_l)
  // *_r = * of neighbor (remote)
  // *_l = * of local domain

  T shift;
  // check if a sensible shift can be made
  if (remote_rank != MPI_PROC_NULL && local_coord != global_dim - 1 &&
      !(remote_work == (T)0 && local_work == (T)0)) {
    shift = 1.0 / gamma * 0.5 * (remote_work - local_work) /
            (remote_work + local_work) * (local_size + remote_size);

    // determine a maximum move in order to avoid overlapping borders and to
    // guarantee a stable shift of the domain borders (avoid the loss of a
    // domain due to too large shifts of its borders [no crossing of borders])
    T maxmove;
    if (shift > 0.0)
      maxmove = 0.49 * (remote_size - minSize);
    else
      maxmove = 0.49 * (local_size - minSize);

    if (std::abs(shift) > maxmove) {
      shift = (shift < 0.0) ? -maxmove : maxmove;
    }
  } else {
    shift = (T)0;
  }

  return shift;
}

}}// namespace ALL::Functions
#endif
